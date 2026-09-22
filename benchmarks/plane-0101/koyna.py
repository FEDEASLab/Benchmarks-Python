#
# Koyna Dam 
#

# from tqdm import tqdm
tqdm = lambda x: x

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

import xara
from xara.helpers import find_node, find_nodes

import veux

cwd = Path.cwd()
if (cwd / 'data').exists():
    KOYNA_DIR = cwd

OUT = KOYNA_DIR / 'out' 
OUT.mkdir(parents=True, exist_ok=True)



def support_nodes(model):
    return [tag for tag in model.getNodeTags() if abs(model.nodeCoord(tag)[1]) < 1.0e-9]


def upstream_face_nodes(model):
    return sorted(
        [(tag, model.nodeCoord(tag)) for tag in model.getNodeTags() if abs(model.nodeCoord(tag)[0]) < 1.0e-9],
        key=lambda item: item[1][1],
    )


def hydrodynamic_loading(model, h_water, rho_water, thickness, gravity):
    face_nodes = upstream_face_nodes(model)
    added_masses = {}
    hydro_loads = {}

    for index, (tag, (_, y)) in enumerate(face_nodes):
        if y > h_water:
            break

        if index == 0:
            tributary = (face_nodes[index + 1][1][1] - y) / 2.0
        elif index == len(face_nodes) - 1:
            tributary = (y - face_nodes[index - 1][1][1]) / 2.0
        else:
            tributary = (face_nodes[index + 1][1][1] - face_nodes[index - 1][1][1]) / 2.0

        mass_per_area = 0.875 * rho_water * np.sqrt(h_water * (h_water - y))
        node_mass = mass_per_area * tributary * thickness
        pressure = rho_water * gravity * (h_water - y)

        added_masses[tag] = (node_mass, 0.0)
        hydro_loads[tag] = (pressure * tributary * thickness, 0.0)

    return added_masses, hydro_loads


def model_edges(model):
    edges = set()
    for ele in model.getEleTags():
        nodes = model.eleNodes(ele)
        for i, j in ((0, 1), (1, 2), (2, 3), (3, 0)):
            edges.add(tuple(sorted((nodes[i], nodes[j]))))
    return sorted(edges)

#
# Model
#

def create_model(koyna, material, element="Q4", units=None, mesh=None):
    if mesh is None:
        mesh = {}
    
    nx      = mesh.get('nx',      20)
    ny_base = mesh.get('ny_base', 24)
    ny_top  = mesh.get('ny_top',  14)

    m = units.meter

    L_base = koyna.L_base
    H_base = koyna.H_base
    H_top  = koyna.H_top
    thickness = koyna.thickness
    h_water = koyna.h_water

    rho = material.asdict().get("density")



    model = xara.Model(ndm=2, ndf=2)

    model.material(material)

    section = xara.PlaneSection("PlaneStress", material=material, thickness=thickness)

    model.section(section)


    plane_args = {"section": section, "b": [0.0, -rho*units.gravity]}

    base_surface = model.surface(
        (nx, ny_base),
        element=element,
        args=plane_args,
        points={
            1: (0.0, 0.0),
            2: (L_base, 0.0),
            3: (19.25 * m, H_base),
            4: (0.0, H_base),
        },
    )

    crest_surface = model.surface(
        (nx, ny_top),
        element=element,
        args=plane_args,
        points={
            1: (    0.0, H_base),
            2: (19.25*m, H_base),
            3: (14.80*m, H_base + H_top),
            4: (    0.0, H_base + H_top),
        },
    )

    for tag in support_nodes(model):
        model.fix(tag, 1, 1)

    added_masses, hydro_loads = hydrodynamic_loading(model, 
                                                     koyna.h_water, 
                                                     koyna.rho_water, 
                                                     koyna.thickness, 
                                                     units.gravity)
    for tag, (mx, my) in added_masses.items():
        model.mass(tag, mx, my)
    return model 

#
# Gravity, Eigenvalues, and Damping
#
def static_analysis(model, added_masses, hydro_loads):

    print(f'Upstream water-loaded nodes: {len(hydro_loads)}')

    model.pattern('Plain', 1,'Constant')
    for tag, (fx, fy) in hydro_loads.items():
        model.load(tag, (fx, fy), pattern=1)


    model.constraints('Plain')
    model.numberer('RCM')
    model.system('BandGeneral')
    model.test('EnergyIncr', 1.0e-12, 50)
    model.algorithm('Newton')
    model.integrator('LoadControl', 1.0)
    model.analysis('Static')

    static_ok = model.analyze(1)
    if static_ok != 0:
        raise RuntimeError(f'Gravity analysis failed with code {static_ok}')

    model.loadConst(time=0.0)


#
# Transient Analysis
#
def dynamic_analysis(koyna, model, accel_x, accel_y, dt_data, duration, units):
    crest = koyna.crest_node(model)
    model.wipeAnalysis()


    model.timeSeries('Path', 11, dt=dt_data, values=accel_x, factor=g)
    model.timeSeries('Path', 12, dt=dt_data, values=accel_y, factor=g)

    model.pattern("UniformExcitation", 11, 1, accel=11)
    model.pattern("UniformExcitation", 12, 2, accel=12)

    model.constraints('Plain')
    model.numberer('RCM')
    model.system('UmfPack')
    model.test('EnergyIncr', 1.0e-12, 5, 1)
    model.algorithm('Newton')
    model.integrator('Newmark', 0.5, 0.25)
    model.analysis('Transient')

    dt = 0.005

    crest_disp=[model.nodeDisp(crest,1)]
    steps = int(duration/dt)
    for i in tqdm(range(steps)):
        status = model.analyze(1, dt)
        if status != xara.successful:
            break
            # raise RuntimeError(f"analysis failed at time {model.getTime()}")
        crest_disp.append(model.nodeDisp(crest,1))


    M = model.getTangent(m=1, c=0, k=0)
    print(f"Mass = {np.sum(M)}")
    print(np.diag(model.getTangent(m=1,c=0,k=0)))


class KoynaDam:
    def __init__(self, units, mesh=None):
        self.units = units
        m  = units.meter
        kg = units.kilogram

        self.L_base = 70.0 * m
        self.H_base = 66.5 * m
        self.H_top  = 36.5 * m
        self.thickness = 1.0 # TODO
        self.h_water = 91.75 * m
        self.rho_water = 1000.0 * kg / (m ** 3)

        self.mesh = mesh

    def crest_node(self, model):
        return find_node(model, x=0, y=self.H_base + self.H_top)

    def create_model(self, material, element="Q4"):
        return create_model(self, 
                            material, 
                            element=element, 
                            units=self.units, 
                            mesh=self.mesh)

    def static_analysis(self, model):
        added_masses, hydro_loads = hydrodynamic_loading(model, 
                                                         self.h_water, 
                                                         self.rho_water, 
                                                         self.thickness, 
                                                         self.units.gravity)
        return static_analysis(model, added_masses, hydro_loads)


    def apply_damping(self, model, zeta=0.03):
        eigenvalues = model.eigen(5)
        periods = [2.0 * np.pi / np.sqrt(value) for value in eigenvalues]

        omega_1 = np.sqrt(eigenvalues[0])
        beta_k = 2.0 * zeta / omega_1
        model.rayleigh(0.0, 0.0, 0.0, beta_k)

        print('Gravity analysis completed.')
        print('First five periods (s):')
        for i, period in enumerate(periods, start=1):
            print(f'  Mode {i}: {period:.6e}')
        print(f'Stiffness-proportional damping beta_k = {beta_k:.6e}')


#!/usr/bin/env python
"""
Cook's membrane

Tapered quadrilateral panel clamped at the left edge with a uniform
in-plane shear traction applied along the right edge. The reported
quantity is the vertical displacement at the midpoint of the tip
edge, point (48, 52). 
The widely cited converged value for the
standard parameters (E=1, nu=1/3, P=1, t=1) is 23.96.
"""
import xara
from xara.load import SurfaceLoad, Line
from xara.helpers import find_nodes, find_node


class CooksMembrane:
    name = "Cook's Membrane"
    reference = 23.96
    metric = "v(48,52)"

    def __init__(self):
        self.points = {
            1: [ 0.0,  0.0, 0.0],
            2: [48.0, 44.0, 0.0],
            3: [48.0, 60.0, 0.0],
            4: [ 0.0, 44.0, 0.0],
        }
        self.E = 1.0
        self.nu = 1.0 / 3.0
        self.thickness = 1.0
        self.shear_total = 1.0
        self.edge_length = 16.0
        self.tip_xy = (48.0, 52.0)

    def solve(self, element, mesh):
        nx, ny = mesh

        model = xara.Model(ndm=3, ndf=6)

        material = xara.MultiaxialMaterial("ElasticIsotropic", E=self.E, nu=self.nu)
        model.material(material)

        section = xara.ShellSection("Elastic", material, self.thickness)
        model.section(section)

        model.surface((nx, ny),
                      element=element,
                      args={"section": section},
                      order=1,
                      points=self.points)

        # Clamp the left edge
        for node in find_nodes(model, x=0.0):
            model.fix(node, (1, 1, 1, 1, 1, 1))

        # Uniform shear traction on the right edge
        q = self.shear_total / self.edge_length
        traction = lambda s: [0.0, q, 0.0]

        right_nodes = sorted(find_nodes(model, x=48.0),
                             key=lambda n: model.nodeCoord(n)[1])
        load = SurfaceLoad(Line(model, right_nodes), traction)

        model.pattern(xara.StaticPattern([load]))

        analysis = xara.StaticAnalysis(model)
        analysis.analyze()

        tip = find_node(model, x=self.tip_xy[0], y=self.tip_xy[1])
        return model.nodeDisp(tip)[1]


def run_validation(problem, elements, meshes):
    title = "=== {} ===".format(problem.name)
    print(title)
    print("Reference {} = {}".format(problem.metric, problem.reference))
    print()

    col_elem, col_mesh, col_disp = 14, 12, 16
    header = ("{:<%d}{:<%d}{:<%d}" % (col_elem, col_mesh, col_disp)).format(
        "Element", "Mesh", "Displacement")
    print(header)
    print("-" * (col_elem + col_mesh + col_disp))

    for element in elements:
        if element is None:
            print()
            continue
        for i, mesh in enumerate(meshes):
            try:
                value = problem.solve(element, mesh)
                value_str = "{:.6f}".format(value)
            except Exception as exc:
                raise exc
                value_str = "FAILED ({})".format(type(exc).__name__)

            mesh_str = "{}x{}".format(mesh[0], mesh[1])
            elem_str = element if i == 0 else ""
            print(("{:<%d}{:<%d}{:<%d}" % (col_elem, col_mesh, col_disp)).format(
                elem_str, mesh_str, value_str))
    print()


if __name__ == "__main__":

    elements = [
        "ShellMITC4",
        "ASDShellQ4",
        "ShellQ4/ASD",
        "ShellQ4/T01",
        "ShellQ4/E5",
        "ShellQ4/F",

        "PlateQ4/L01",
        "PlateQ4/L02",
        "PlateQ4/U",
        "PlateQ4/E5",
    ]
    meshes = [
        (2, 2), 
        # (8, 8), 
        # (16, 16), 
        # (32, 32)
    ]

    problems = [
        (CooksMembrane(), meshes),
    ]

    for problem, meshes in problems:
        run_validation(problem, elements, meshes)

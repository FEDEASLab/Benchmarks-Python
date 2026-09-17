import xara
from xara.helpers import find_node
from xara.units.mks import cm, MPa, N, m
import matplotlib.pyplot as plt
import veux, veux.motion


class Span15:
    # Encastred beam with distributed load
    def __init__(self, w, L, E, I, k=1, A=0, G=0):
        self.w = w
        self.L = L
        self.E = E
        self.I = I
        self.k = k

        self.uy = lambda x: \
            w*x**2/(24*E*I)*(L - x)**2 - \
                ((w/(k*G*A)*L**2/24*(1 - 12*x/L + 12*(x/L)**2) + w*L**2/(24*G*A*k)) if G and A else 0)


class Span1:
    # Simple span with distributed load
    def __init__(self, w, L, E, I, k=1, A=0, G=0):
        self.k = k
        self.A = A
        self.G = G
        self.w = w
        self.L = L
        self.E = E
        self.I = I

    def uy(self, x):
        L = self.L
        E = self.E
        I = self.I
        w = self.w
        k = self.k
        G = self.G
        A = self.A
        xi = x/L
        return w/(24*E*I)*(L**3*x - 2*L*x**3 + x**4) \
           + (
               w*L**2/(2*k*G*A)*(xi - xi**2) if G*A else 0
           )

def create_section(model, type, tag, shape, material):

    model.material('ElasticIsotropic', 1, **material)

    if type == "Elastic":
        cmm = shape.cmm()
        cnn = shape.cnn()
        cnv = shape.cnv()
        cnm = shape.cnm()
        cmw = shape.cmw()
        A = cnn[0,0]
        model.section("ElasticFrame", tag,
                        E=E,
                        G=G,
                        A=A,
                        Ay=1*A,
                        Az=1*A,
                        Qy=cnm[0,1],
                        Qz=cnm[2,0],
                        Iy=cmm[1,1],
                        Iz=cmm[2,2],
                        J =shape._analysis.torsion_constant(),
                        Ry= cnv[1,0],
                        Rz=-cnv[2,0],
                        Sy= cmw[1,0],
                        Sz=-cmw[2,0]
        )
    else:

        model.section("ShearFiber", tag, GJ=0)

        for fiber in shape.create_fibers():

            model.fiber(**fiber, material=1, section=1)


def create_prism(length:    float,
                 element:   str,
                 shape:     dict,
                 material:  dict,
                 boundary:  tuple,
                 section:   str = "ShearFiber",
                 vertical:  int = 3,
                 shear:     int = 0,
                 order:     int = 0,
             #   orient = (0,  -1, 0)
                 geometry:  str = None,
                 transform: str = None,
                 divisions: int = 1,
                 rotation = None):

    L  = length

    # Number of elements discretizing the column
    ne = divisions

    nen = order + 2
    nn = ne*(nen-1)+1

    model = xara.Model(ndm=3, ndf=6)

    for i in range(1, nn+1):
        x = (i-1)*L/float(nn-1)
        location = (x, 0, 0)

        if rotation is not None:
            location = tuple(rotation@location)

        model.node(i, location)


    # Define boundary conditions
    model.fix( 1, boundary[0])
    model.fix(nn, boundary[1])

    #
    # Define cross-section 
    #
    create_section(model, section, 1, shape, material)

    # Define geometric transformation
    if vertical == 3:
        orient = (0, 0, 1)
    if rotation is not None:
        orient = tuple(map(float, rotation@orient))
    model.geomTransf(transform, 1, *orient)


    # Define elements
    for i in range(ne):
        start = i * (nen - 1) + 1
        nodes = list(range(start, start + nen))
        if geometry is None or geometry == "Linear" or "Exact" in element:
            model.element(element, i+1, nodes,
                        section=1,
                        shear=shear,
                        transform=1)
        else:
            model.element(element, i, nodes,
                        section=1,
                        order={"Linear": 0, "delta": 1}[geometry],
                        transform=1)

    return model


if __name__ == "__main__":
    from xsection.library import Rectangle 

    length = 50*cm
    E = 220e3 * MPa
    G = 0.5*E/(1 + 0.3)
    w = 2.0e4 * N/m

    shape = Rectangle(d=1*cm, b=10*cm)
    I = shape.cmm()[1,1]
    A = shape.cnn()[0,0]

    divisions = 8
    assert divisions % 2 == 0, "Number of divisions must be even"

    model = create_prism(
        length=length,
        element="ExactFrame",
        shape=shape,
        material = {"E": E, "G": G},
        vertical = 3,
        boundary = ((1, 1, 1, 1, 0, 0),
                    (1, 1, 1, 1, 0, 0)),
        section="Elastic",
        transform="Linear",
        shear=1,
        order=1,
        divisions = divisions
    )

    model.pattern("Plain", 1, "Linear")

    if False:
        for i in range(divisions):
            model.eleLoad("-type", "BeamUniform", (0.0, w), ele=i+1, pattern=1)
    else:
        model.eleLoad("Frame", "Heaviside",
                    basis = "local",
                    force = [0, 0, w],
                    pattern=1,
                    elements=list(range(1, divisions+1))
        )


    steps = 20
    model.system('Umfpack')
    model.integrator("LoadControl", 1/steps)
    model.test("Energy", 1e-16, 10, 1)
    model.algorithm("Newton")
    model.analysis("Static")

    artist = veux.create_artist(model, vertical=3,
                                model_config=dict(extrude_outline=shape))
    motion = veux.motion.Motion(artist)

    u = [0]
    P = [0]
    for i in range(steps):

        if model.analyze(1) != 0:
            print(f"Failed at time = {model.getTime()}")
            break

        u.append(model.nodeDisp(find_node(model, x=length/2), 3))
        P.append(model.getTime())


    uy = Span1(w=w, L=length, E=E, I=I, k=1, A=A, G=G).uy(length/2)
    fig, ax = plt.subplots()
    ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.axvline(0, color='black', linestyle='-', linewidth=1)
    ax.plot(u, P)
    ax.plot([0,uy], [0,1.0], color='red', linestyle='--', linewidth=1, label="Linear")
    plt.show()

    model.reactions()
    print(
        w*length / sum(model.nodeResponse(i, 3, "reactionForce") for i in model.getNodeTags()),
    )

    motion.add_to(artist.canvas)
    veux.serve(artist)

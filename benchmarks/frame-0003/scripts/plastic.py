from xara.benchmarks import Prism
from xara.helpers import find_node
from xara.units.mks import cm, MPa, N, m
from xsection.library import Rectangle 
import matplotlib.pyplot as plt
import veux, veux.motion

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


def analyze(model, prism, w):
    divisions = prism.divisions
    model.pattern("Plain", 1, "Linear")

    model.eleLoad("Frame", "Heaviside",
                basis = "local",
                force = [0, 0, w],
                pattern=1,
                elements=list(range(1, divisions))
    )


    steps = 40
    model.system('Umfpack')
    model.integrator("LoadControl", 1/steps)
    model.test("Energy", 1e-14, 40, 1)
    model.algorithm("Newton")
    model.analysis("Static")


    u = [0]
    P = [0]
    for i in range(steps):

        if model.analyze(1) != 0:
            print(f"Failed at time = {model.getTime()}")
            break

        u.append(model.nodeDisp(find_node(model, x=length/2), 3))
        P.append(model.getTime())

    return u,P

if __name__ == "__main__":

    length = 50*cm
    w  = 2e6 * N/m

    #
    # Material
    #
    E  = 220e3*MPa
    G  = E/(1 + 0.27)/2
    Fy = 300*MPa

    #
    # Section
    #
    shape = Rectangle(d=5*cm, b=10*cm, mesh_scale=1/800)
    I = shape.cmm()[1,1]
    A = shape.cnn()[0,0]

    divisions = 4
    assert divisions % 2 == 0, "Number of divisions must be even"

    #
    # Reference solution
    #
    uy = Span1(w=w, L=length, E=E, I=I, k=1, A=A, G=G).uy(length/2)
    fig, ax = plt.subplots()
    ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.axvline(0, color='black', linestyle='-', linewidth=1)
    ax.plot([0,uy], [0,1.0], color='red', linestyle='--', linewidth=1, label="Linear")
    markers = iter(["s",  "+", "x", ".", "^", "<", ">"])
    

    #
    # Finite Element Analysis
    #
    for element in ["ForceFrame", "ExactFrame"]:
        for section in ["Elastic", "ShearFiber"]:
            for shear in [0,1]:

                if "Exact" in element and shear == 0:
                    continue

                prism = Prism(
                    length=length,
                    element=element,
                    section=section,
                    shape=shape,
                    material = {
                        "E": E, 
                        "G": G,
                        "name": "J2", # "J2Simplified"
                        "Fy": Fy, 
                        "Hiso": 0.02*E,
                        "Fs": Fy, "Hsat": 0
                    },
                    vertical = 3,
                    boundary = ((1, 1, 1, 1, 0, 0),
                                (0, 1, 1, 1, 0, 0)),
                    transform="Corotational02",
                    shear=shear,
                    order=1 if "Exact" not in element else 3,
                    divisions = divisions
                )

                model = prism.create_model()
                u, P = analyze(model, prism, w)
                ax.plot(u, P, next(markers), label=f"{element}, {section}, shear: {bool(shear)}")

    ax.legend()
    ax.set_ylabel("$\\lambda$")
    ax.set_xlabel("$u$ [m]")
    ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    plt.show()


    model.reactions()
    print(
      w*length / sum(model.nodeResponse(i, 3, "reactionForce") for i in model.getNodeTags()),
    )
    print(f"{E = }, {G = }")

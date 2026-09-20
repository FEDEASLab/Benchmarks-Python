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



def analyze(model, prism, w, basis):
    divisions = prism.divisions
    model.pattern("Plain", 1, "Linear")

    if basis is None:
        for i in range(divisions):
            model.eleLoad("-type", "BeamUniform", (0.0, w), 
                          ele=i+1, pattern=1)
    
    else:
        model.eleLoad("Frame", "Heaviside",
                    basis = basis,
                    force = [0, 0, w],
                    pattern=1,
                    elements=list(range(1, divisions))
        )


    steps = 40
    model.system('Umfpack')
    model.integrator("LoadControl", 1/steps)
    model.test("Energy", 1e-16, 40, 0)
    model.algorithm("Newton")
    model.analysis("Static")

    # artist = veux.create_artist(model, vertical=3,
    #                             model_config=dict(extrude_outline=shape))
    # motion = veux.motion.Motion(artist)

    u = [0]
    P = [0]
    for i in range(steps):

        if model.analyze(1) != 0:
            print(f"Failed at time = {model.getTime()}")
            raise RuntimeError(f"Analysis failed at step {i+1}")
            break

        # motion.draw_sections(rotation=model.nodeRotation,
        #                      position=model.nodeDisp)
        # motion.advance(time=model.getTime()*10)

        u.append(model.nodeDisp(find_node(model, x=length/2), 3))
        P.append(model.getTime())

    return u, P


if __name__ == "__main__":

    length = 50*cm
    E = 220e3 * MPa
    G = 0.5*E/(1 + 0.27)
    w = 1.0e6 * N/m
    shear = 1

    shape = Rectangle(d=1*cm, b=100*cm)

    I = shape.cmm()[1,1]
    A = shape.cnn()[0,0]


    uy = Span1(w=w, L=length, E=E, I=I, k=1, A=A, G=G*shear).uy(length/2)
    fig, ax = plt.subplots()
    ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.axvline(0, color='black', linestyle='-', linewidth=1)
    ax.plot([0,uy], [0,1.0], color='red', linestyle='--', linewidth=1, label="Linear")
    markers = iter(["s", "^", "v", "o", "+", "x", ".", "^", "<", ">"])


    c = 1
    for basis in ["director", "reference"]:
        for (element, n_e) in [("ExactFrame",  2),  ("ForceFrame",  2), 
                               ("ExactFrame", 16)]: #,  ("ForceFrame", 16)]:

            assert n_e % 2 == 0, "Number of divisions must be even"


            if "Exact" in element and basis is None:
                continue 

            prism = Prism(
                length=length,
                element=element,
                shape=shape,
                material = {"E": E, "G": G},
                vertical = 3,
                boundary = ((1, 1, 1, 1, 0, 0),
                            (c, 1, 1, 1, 0, 0)),
                section="Elastic",
                transform="Corotational02",
                shear=shear,
                order=1 if "Exact" not in element else 3,
                divisions = n_e
            )

            model = prism.create_model()
            u, P = analyze(model, prism, w, basis)
            
            ax.plot(u, P, 
                    ":" if n_e > 10 else next(markers), 
                    label=f"{element} (${n_e = }$), {basis} load")

    ax.legend()
    ax.grid()
    ax.set_ylabel("$\\lambda$")
    ax.set_xlabel("$u$ [m]")
    ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    plt.show()


    # motion.add_to(artist.canvas)
    # veux.serve(artist)

import veux
from veux.config import LineStyle
import numpy as np
import xara.units.iks as units
from xara.units.iks import foot, inch
from xsection.library import from_aisc
from veux.utility.alpha_shape import alpha_shape
import matplotlib.pyplot as plt


if __name__ == "__main__":

    shape = from_aisc("W14x48", units=units)
    print(shape.summary())


    points = np.array([
        fiber.location for fiber in shape.create_fibers(origin="centroid", mesh_scale=1/2.5)
    ])

    fig, ax = plt.subplots()
    ax.scatter(*zip(*points), 0.2, marker=".")
    ax.set_aspect("equal")
    alpha = alpha_shape(points, bound_ratio=0.01) #radius=0.05)
    ax.scatter(*zip(*alpha), 0.2, marker=".", color="red")
    plt.show()


    mesh = shape.mesh
    artist = veux.create_artist((mesh.nodes, mesh.cells()), ndf=1)
    artist.draw_surfaces()
    d = 8*inch
#   R = np.eye(3)
    R = np.array([[0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 0]])
    for i,v in enumerate(d/5*R.T[1:]):
        artist.canvas.plot_vectors([[0,0,0]], [v],
                                   line_style=LineStyle(color=("green", "blue")[i]),
                                   extrude=True)
#   artist.draw_origin(extrude=True)
    artist.draw_outlines()
    veux.serve(artist)

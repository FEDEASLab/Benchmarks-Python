#
# Bathe's curved cantilever analysis with Crisfield-Jelenic load paths
#
# [1] Jelenić G, Crisfield MA (1999) "Geometrically exact 3D beam theory:
#     implementation of a strain-invariant finite element for statics and
#     dynamics."
#     Computer Methods in Applied Mechanics and Engineering,  171(1–2):141–171. 
#     https://doi.org/10/dj37b3
#
# [2] Simo, J. C., and L. Vu-Quoc. 
#     "A Three-Dimensional Finite-Strain Rod Model. Part II: Computational Aspects."
#     Computer Methods in Applied Mechanics and Engineering 58, no. 1 (1986): 79–116. 
#     https://doi.org/10.1016/0045-7825(86)90079-4.
#

from xara.examples.perez2024nonlinear import Ex_4_4

def test_perez():
    Ex_4_4()


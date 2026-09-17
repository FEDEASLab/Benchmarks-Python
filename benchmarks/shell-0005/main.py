#!/usr/bin/env python
"""
Plate with rectangular cut-out, line load on the inner top edge.

Geometry:
  Plate     500 x 300 mm
  Cut-out   250 x 180 mm, centred (corners at (125, 60) and (375, 240))
  Thickness 20 mm
Material:
  E = 200 000 N/mm^2,  nu = 0.25
BCs:
  Left edge  (x = 0)   clamped, all 6 DOFs
  Top edge   (y = 300) clamped, all 6 DOFs
  Right and bottom edges free
Load:
  1 kN/mm directed along -z on the inner top edge of the cut-out (y = 240,
  125 <= x <= 375).
Quantity of interest:
  Transverse deflection w at point A = (375, 60), the bottom-right hole corner.
"""

from __future__ import annotations

import argparse

import xara
from xara.helpers import find_node
from xara.load import Line, SurfaceLoad


# Geometry (mm)
PLATE_WIDTH_MM = 500.0
PLATE_HEIGHT_MM = 300.0

# Rectangular hole: axis-aligned, centred in the plate.
HOLE_X_MIN_MM = 125.0
HOLE_X_MAX_MM = 375.0
HOLE_Y_MIN_MM = 60.0
HOLE_Y_MAX_MM = 240.0

THICKNESS_MM = 20.0

# Benchmark point A: bottom-right corner of the hole (see module docstring).
POINT_A_X_MM = HOLE_X_MAX_MM
POINT_A_Y_MM = HOLE_Y_MIN_MM

# --- Material ---------------------------------------------------------------
E_MODULUS_MPA = 200_000.0  # N/mm^2
POISSON_RATIO = 0.25

# --- Loading ---------------------------------------------------------------
# Uniform line load on the inner top edge of the hole, global +z component (N/mm).
LINE_LOAD_Z_N_PER_MM = -1000.0

# Nodes closer than this (mm) to a boundary line are treated as on that line.
COORD_MATCH_TOL_MM = 1e-6


# Eight quadrilateral patches cover the plate; the hole is left empty.
# Keys: row (B=bottom, M=middle, T=top) + column (L=left, M=middle, R=right).
# Each value is four (x, y) corners in counter-clockwise order on the mid-surface.
_PATCH_CORNERS_MM: dict[str, list[tuple[float, float]]] = {
    "BL": [
        (0.0, 0.0),
        (HOLE_X_MIN_MM, 0.0),
        (HOLE_X_MIN_MM, HOLE_Y_MIN_MM),
        (0.0, HOLE_Y_MIN_MM),
    ],
    "BM": [
        (HOLE_X_MIN_MM, 0.0),
        (HOLE_X_MAX_MM, 0.0),
        (HOLE_X_MAX_MM, HOLE_Y_MIN_MM),
        (HOLE_X_MIN_MM, HOLE_Y_MIN_MM),
    ],
    "BR": [
        (HOLE_X_MAX_MM, 0.0),
        (PLATE_WIDTH_MM, 0.0),
        (PLATE_WIDTH_MM, HOLE_Y_MIN_MM),
        (HOLE_X_MAX_MM, HOLE_Y_MIN_MM),
    ],
    "ML": [
        (0.0, HOLE_Y_MIN_MM),
        (HOLE_X_MIN_MM, HOLE_Y_MIN_MM),
        (HOLE_X_MIN_MM, HOLE_Y_MAX_MM),
        (0.0, HOLE_Y_MAX_MM),
    ],
    "MR": [
        (HOLE_X_MAX_MM, HOLE_Y_MIN_MM),
        (PLATE_WIDTH_MM, HOLE_Y_MIN_MM),
        (PLATE_WIDTH_MM, HOLE_Y_MAX_MM),
        (HOLE_X_MAX_MM, HOLE_Y_MAX_MM),
    ],
    "TL": [
        (0.0, HOLE_Y_MAX_MM),
        (HOLE_X_MIN_MM, HOLE_Y_MAX_MM),
        (HOLE_X_MIN_MM, PLATE_HEIGHT_MM),
        (0.0, PLATE_HEIGHT_MM),
    ],
    "TM": [
        (HOLE_X_MIN_MM, HOLE_Y_MAX_MM),
        (HOLE_X_MAX_MM, HOLE_Y_MAX_MM),
        (HOLE_X_MAX_MM, PLATE_HEIGHT_MM),
        (HOLE_X_MIN_MM, PLATE_HEIGHT_MM),
    ],
    "TR": [
        (HOLE_X_MAX_MM, HOLE_Y_MAX_MM),
        (PLATE_WIDTH_MM, HOLE_Y_MAX_MM),
        (PLATE_WIDTH_MM, PLATE_HEIGHT_MM),
        (HOLE_X_MAX_MM, PLATE_HEIGHT_MM),
    ],
}

# Physical width/height of each patch column (L/M/R) and row (B/M/T), in mm.
_PATCH_COLUMN_WIDTH_MM = {
    "L": HOLE_X_MIN_MM,
    "M": HOLE_X_MAX_MM - HOLE_X_MIN_MM,
    "R": PLATE_WIDTH_MM - HOLE_X_MAX_MM,
}
_PATCH_ROW_HEIGHT_MM = {
    "B": HOLE_Y_MIN_MM,
    "M": HOLE_Y_MAX_MM - HOLE_Y_MIN_MM,
    "T": PLATE_HEIGHT_MM - HOLE_Y_MAX_MM,
}


def _surface_subdivisions(
    patch_key: str, target_edge_length_mm: float
) -> tuple[int, int]:
    """
    Number of finite elements along each parent direction for one patch.

    ``patch_key`` is two letters: row (B/M/T) then column (L/M/R).  The counts
    are chosen so that the shorter patch dimension is near ``target_edge_length_mm``.
    """
    row_letter, col_letter = patch_key[0], patch_key[1]
    num_elem_x = max(
        1, round(_PATCH_COLUMN_WIDTH_MM[col_letter] / target_edge_length_mm)
    )
    num_elem_y = max(
        1, round(_PATCH_ROW_HEIGHT_MM[row_letter] / target_edge_length_mm)
    )
    return num_elem_x, num_elem_y


def _corner_points_for_xara_surface(
    corners_xy: list[tuple[float, float]],
) -> dict[int, list[float]]:
    """Map 2-D patch corners to Xara's 1-based point dict with z = 0."""
    return {
        index + 1: [float(x), float(y), 0.0]
        for index, (x, y) in enumerate(corners_xy)
    }


def build_model(element: str, order: int, target_edge_length_mm: float) -> xara.Model:
    """
    Assemble the mesh, boundary conditions, and line load for this benchmark.

    Parameters
    ----------
    element
        Shell element type passed to ``Model.surface``.
    order
        Parent-element order (1 for bilinear, 2 for biquadratic).
    target_edge_length_mm
        Target physical edge length (mm); actual spacing is rounded per patch.

    Returns
    -------
    xara.Model
        Fully constrained and loaded 3-D shell model (6 DOFs per node).
    """
    model = xara.Model(ndm=3, ndf=6)

    material = xara.TriaxialMaterial(
        "ElasticIsotropic", E=E_MODULUS_MPA, nu=POISSON_RATIO
    )
    model.material(material)

    section = xara.ShellSection("Elastic", material, THICKNESS_MM)
    model.section(section)

    # Coincident patch corners share nodes automatically.
    for patch_key, corners_xy in _PATCH_CORNERS_MM.items():
        model.surface(
            _surface_subdivisions(patch_key, target_edge_length_mm),
            element=element,
            args={"section": section},
            order=order,
            points=_corner_points_for_xara_surface(corners_xy),
        )

    # Clamped outer edges: left (x = 0) and top (y = plate height).
    for node_tag in model.getNodeTags():
        x_mm, y_mm, _z_mm = model.nodeCoord(node_tag)
        if abs(x_mm) < COORD_MATCH_TOL_MM or abs(
            y_mm - PLATE_HEIGHT_MM
        ) < COORD_MATCH_TOL_MM:
            model.fix(node_tag, (1, 1, 1, 1, 1, 1))

    # Mindlin plate elements carry a drilling DOF; pin it at every node.
    if "plate" in element.lower():
        for node_tag in model.getNodeTags():
            try:
                model.fix(node_tag, (0, 0, 0, 0, 0, 1))
            except Exception:
                # Node may already have a conflicting fixity from the edge clamps.
                pass

    # Inner top edge of the hole: y = HOLE_Y_MAX_MM, x between hole x-bounds.
    cutout_top_edge_node_tags = [
        tag
        for tag in model.getNodeTags()
        if abs(model.nodeCoord(tag)[1] - HOLE_Y_MAX_MM) < COORD_MATCH_TOL_MM
        and HOLE_X_MIN_MM - COORD_MATCH_TOL_MM
        <= model.nodeCoord(tag)[0]
        <= HOLE_X_MAX_MM + COORD_MATCH_TOL_MM
    ]
    cutout_top_edge_node_tags.sort(key=lambda t: model.nodeCoord(t)[0])

    def uniform_line_load_global(*_coords: object) -> list[float]:
        """Line-load intensity callback: constant (0, 0, q_z) in global axes."""
        return [0.0, 0.0, LINE_LOAD_Z_N_PER_MM]

    distributed_load = SurfaceLoad(
        Line(model, cutout_top_edge_node_tags), uniform_line_load_global
    )
    model.pattern(xara.StaticPattern([distributed_load]))
    return model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "-e",
        "--element",
        default="ASDShellQ4",
        choices=["ASDShellQ4", "ShellMITC4", "HeterosisPlate", "ShellMITC9"],
        help="shell element name (default: ASDShellQ4)",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=1,
        choices=(1, 2),
        help="mesh order, 1 for Q4 or 2 for Q9 (default: 1)",
    )
    parser.add_argument(
        "--h",
        type=float,
        default=25.0,
        help="target element edge length in mm (default: 25)",
    )
    args = parser.parse_args()

    model = build_model(args.element, args.order, args.h)
    print(f"  Element: {args.element}")
    print(f"  {len(model.getNodeTags())} nodes")

    analysis = xara.StaticAnalysis(model, system="Umfpack")
    print(analysis)
    analysis.analyze()

    node_a = find_node(model, x=POINT_A_X_MM, y=POINT_A_Y_MM)
    w_vertical_mm = model.nodeDisp(node_a)[2]

    print(f"  w at point A     = {w_vertical_mm:.6e} mm")

    import veux

    veux.serve(veux.render(model))


if __name__ == "__main__":
    main()

"""
Clamped circular plate under a center point load.

The benchmark compares nodal transverse displacement profiles against the
closed-form Reissner-Mindlin solution for a clamped circular plate under a
center point load.

The HeterosisPlate case is handled separately because its external OpenSees
connectivity uses 9 shell nodes, while the physical plate interpolation uses
Q8 transverse displacement and Q9 rotations.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest
import xara
from opensees.errors import XaraError
from xara.helpers import find_node


RADIUS = 5.0
THICKNESS = 2.0
E_MOD = 10.92e5
NU = 0.3
KAPPA = 5.0 / 6.0

P_TOTAL = 1.0
P_QUARTER = -P_TOTAL / 4.0

SECTION_TAG = 1
LOAD_PATTERN_TAG = 1
MESH_SUBDIVISIONS_PER_DIRECTION = 12

MAX_POINTWISE_REL_ERROR = 5.0e-2
MAX_POINTWISE_ABS_ERROR = 1.0e-2
REFERENCE_VALUE_FLOOR = 5.0e-2

ELEMENT_ORDER = {
    "ASDShellQ4": 1,
    "ShellMITC4": 1,
    "HeterosisPlate": 2,
    "ShellMITC9": 2,
}


def bending_rigidity() -> float:
    """Return plate flexural rigidity D = E t^3 / (12 * (1 - nu^2))."""
    return E_MOD * THICKNESS**3 / (12.0 * (1.0 - NU**2))


def shear_rigidity() -> float:
    """Return Mindlin shear stiffness kappa * G * t."""
    shear_modulus = E_MOD / (2.0 * (1.0 + NU))
    return KAPPA * shear_modulus * THICKNESS


def normalized_displacement(w_down: np.ndarray) -> np.ndarray:
    """Return Hughes normalized displacement, with w positive downward."""
    return 16.0 * math.pi * bending_rigidity() * w_down / (
        P_TOTAL * RADIUS**2
    )


def reissner_clamped_center_load_w(
    radius_from_center: np.ndarray,
) -> np.ndarray:
    """Return Reissner-Mindlin displacement for a clamped circular plate."""
    rho = np.asarray(radius_from_center, dtype=float) / RADIUS
    rho = np.clip(rho, 1.0e-14, 1.0)

    flexural = bending_rigidity()
    shear = shear_rigidity()

    w_kirchhoff = (
        P_TOTAL
        * RADIUS**2
        / (16.0 * math.pi * flexural)
        * (1.0 - rho**2 + 2.0 * rho**2 * np.log(rho))
    )
    w_shear = P_TOTAL / (2.0 * math.pi * shear) * np.log(1.0 / rho)

    return w_kirchhoff + w_shear


def quarter_disk_points() -> dict[int, list[float]]:
    """Return Q9 control points for the mapped quarter disk."""
    r = RADIUS
    a = math.pi / 8.0

    return {
        1: [0.0, 0.0, 0.0],
        2: [r, 0.0, 0.0],
        3: [r / math.sqrt(2.0), r / math.sqrt(2.0), 0.0],
        4: [0.0, r, 0.0],
        5: [0.5 * r, 0.0, 0.0],
        6: [r * math.cos(a), r * math.sin(a), 0.0],
        7: [r * math.cos(3.0 * a), r * math.sin(3.0 * a), 0.0],
        8: [0.0, 0.5 * r, 0.0],
        9: [0.5 * r / math.sqrt(2.0), 0.5 * r / math.sqrt(2.0), 0.0],
    }


def add_to_fixity(fixity: list[int], dofs_1based: tuple[int, ...]) -> None:
    """Set selected 1-based degrees of freedom in a 6-entry fixity vector."""
    for dof in dofs_1based:
        fixity[dof - 1] = 1


def nodes_from_edges(edge_segments: list[tuple[int, int]]) -> list[int]:
    """Return sorted unique node tags from a sequence of edge node pairs."""
    return sorted({node for edge in edge_segments for node in edge})


def is_heterosis_plate(element_name: str) -> bool:
    """Return true for the HeterosisPlate element."""
    return element_name.lower() == "heterosisplate"


def heterosis_rotation_only_nodes(model: xara.Model) -> set[int]:
    """
    Return nodes that belong to the Q9 rotation field but not to the Q8 w field.

    HeterosisPlate uses the first 8 connectivity nodes for transverse
    displacement interpolation and all 9 connectivity nodes for rotation
    interpolation. Therefore, the 9th local node is rotation-only unless the
    same global node appears in the first 8 positions of another element.
    """
    w_nodes: set[int] = set()
    rotation_center_nodes: set[int] = set()

    for element_tag in model.getEleTags():
        connectivity = [int(tag) for tag in model.eleNodes(element_tag)]

        if len(connectivity) == 9:
            w_nodes.update(connectivity[:8])
            rotation_center_nodes.add(connectivity[8])

    return rotation_center_nodes - w_nodes


def apply_boundary_conditions(
    model: xara.Model,
    surface: Any,
    element_name: str,
    n_el: int,
) -> None:
    """Apply symmetry, clamp, and Heterosis inactive-degree constraints."""
    edge_node_pairs = list(surface.walk_edge())

    x_symmetry_nodes = nodes_from_edges(edge_node_pairs[:n_el])
    clamped_outer_arc_nodes = nodes_from_edges(edge_node_pairs[n_el : 3 * n_el])
    y_symmetry_nodes = nodes_from_edges(edge_node_pairs[3 * n_el :])

    node_fixities = {
        node_tag: [0, 0, 0, 0, 0, 0]
        for node_tag in model.getNodeTags()
    }

    for node_tag in x_symmetry_nodes:
        add_to_fixity(node_fixities[node_tag], (2, 4))

    for node_tag in y_symmetry_nodes:
        add_to_fixity(node_fixities[node_tag], (1, 5))

    for node_tag in clamped_outer_arc_nodes:
        add_to_fixity(node_fixities[node_tag], (3, 4, 5))

    if is_heterosis_plate(element_name):
        for node_tag in model.getNodeTags():
            add_to_fixity(node_fixities[node_tag], (1, 2, 6))

        for node_tag in heterosis_rotation_only_nodes(model):
            add_to_fixity(node_fixities[node_tag], (3,))

    elif "plate" in element_name.lower():
        for node_tag in model.getNodeTags():
            add_to_fixity(node_fixities[node_tag], (6,))

    for node_tag, fixity in node_fixities.items():
        if any(fixity):
            model.fix(node_tag, tuple(fixity))


def build_model(element_name: str, order: int, n_el: int) -> tuple[xara.Model, int]:
    """Build the shell-0001 model and return the model and center node tag."""
    model = xara.Model(ndm=3, ndf=6)

    model.section("ElasticShell", SECTION_TAG, E_MOD, NU, THICKNESS)

    surface = model.surface(
        (n_el, n_el),
        element=element_name,
        args={"section": SECTION_TAG},
        order=order,
        points=quarter_disk_points(),
    )

    apply_boundary_conditions(model, surface, element_name, n_el)

    center = find_node(model, x=0.0, y=0.0)

    model.pattern("Plain", LOAD_PATTERN_TAG, "Linear")
    model.load(
        center,
        (0.0, 0.0, P_QUARTER, 0.0, 0.0, 0.0),
        pattern=LOAD_PATTERN_TAG,
    )

    return model, center


def nodal_profile(
    model: xara.Model,
    excluded_node_tags: set[int] | None = None,
) -> np.ndarray:
    """Return sorted nodal samples (rho, radius, w_down)."""
    excluded_node_tags = excluded_node_tags or set()
    samples = []

    for tag in model.getNodeTags():
        if tag in excluded_node_tags:
            continue

        x, y, _ = model.nodeCoord(tag)
        radius = math.hypot(x, y)
        rho = radius / RADIUS

        if 0.05 < rho < 0.95:
            w_down = -model.nodeDisp(tag)[2]
            samples.append((rho, radius, w_down))

    if not samples:
        raise RuntimeError("No nodal profile points found in 0.05 < r/R < 0.95.")

    return np.asarray(sorted(samples, key=lambda row: row[0]), dtype=float)


def assert_pointwise_match(
    model: xara.Model,
    element_name: str,
    mesh_order: int,
    max_relative_error: float,
    max_absolute_error: float,
    excluded_node_tags: set[int] | None = None,
) -> None:
    """Check the normalized displacement profile against the reference solution."""
    profile = nodal_profile(model, excluded_node_tags=excluded_node_tags)
    radius = profile[:, 1]

    y_fem = normalized_displacement(profile[:, 2])
    y_ref = normalized_displacement(reissner_clamped_center_load_w(radius))

    pointwise_abs_error = np.abs(y_fem - y_ref)
    pointwise_rel_error = pointwise_abs_error / np.maximum(
        np.abs(y_ref),
        REFERENCE_VALUE_FLOOR,
    )

    passing_mask = (
        (pointwise_rel_error <= max_relative_error)
        | (pointwise_abs_error <= max_absolute_error)
    )
    failing_indices = np.where(~passing_mask)[0]

    assert failing_indices.size == 0, (
        f"{element_name} (order={mesh_order}, "
        f"n_el={MESH_SUBDIVISIONS_PER_DIRECTION}) has "
        f"{failing_indices.size}/{len(pointwise_abs_error)} nodes violating "
        f"rel<={max_relative_error:.3e} OR abs<={max_absolute_error:.3e} "
        f"(reference floor={REFERENCE_VALUE_FLOOR:.3e}); "
        f"max rel={pointwise_rel_error.max():.6e}, "
        f"max abs={pointwise_abs_error.max():.6e}"
    )


@pytest.mark.parametrize(
    "element_name,mesh_order",
    [
        (name, order)
        for name, order in sorted(ELEMENT_ORDER.items())
        if name != "HeterosisPlate"
    ],
)
def test_hughes_clamped_plate_pointwise_error(
    element_name: str,
    mesh_order: int,
) -> None:
    """Verify shell elements node-by-node against the Reissner reference."""
    try:
        model, _ = build_model(
            element_name,
            mesh_order,
            MESH_SUBDIVISIONS_PER_DIRECTION,
        )
    except XaraError as exc:
        if "unknown element type" in str(exc).lower():
            pytest.skip(f"{element_name} is not available in this OpenSees build.")
        raise

    analysis = xara.StaticAnalysis(model, system="Umfpack")
    analysis.analyze()

    assert_pointwise_match(
        model=model,
        element_name=element_name,
        mesh_order=mesh_order,
        max_relative_error=MAX_POINTWISE_REL_ERROR,
        max_absolute_error=MAX_POINTWISE_ABS_ERROR,
    )


def test_hughes_clamped_plate_pointwise_error_heterosis() -> None:
    """Verify HeterosisPlate using connectivity-aware displacement sampling."""
    element_name = "HeterosisPlate"
    mesh_order = ELEMENT_ORDER[element_name]

    try:
        model, _ = build_model(
            element_name,
            mesh_order,
            MESH_SUBDIVISIONS_PER_DIRECTION,
        )
    except XaraError as exc:
        if "unknown element type" in str(exc).lower():
            pytest.skip(f"{element_name} is not available in this OpenSees build.")
        raise

    analysis = xara.StaticAnalysis(model, system="Umfpack")
    analysis.analyze()

    excluded_nodes = heterosis_rotation_only_nodes(model)

    assert_pointwise_match(
        model=model,
        element_name=element_name,
        mesh_order=mesh_order,
        max_relative_error=MAX_POINTWISE_REL_ERROR,
        max_absolute_error=MAX_POINTWISE_ABS_ERROR,
        excluded_node_tags=excluded_nodes,
    )
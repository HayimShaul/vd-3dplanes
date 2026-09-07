"""Types for the zone of a query line in a 2D arrangement."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.arrangement2d.types import EdgePiece, Face, Vertex
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar


class ZoneError(ValueError):
    """Base for zone degeneracies."""


class QueryOverlapsArrangement(ZoneError):
    """Query line is coincident with an arrangement line.

    Phase 10 handles this inside ``compute_zone`` / ``compute_crossings``.
    The exception remains for callers that want to reject the case.
    """


class QueryThroughVertex(ZoneError):
    """Query line passes through an arrangement vertex.

    Phase 10 handles this inside ``compute_zone``. The exception remains
    for callers that want to reject the case.
    """


class PointNotInOpenFace(ZoneError):
    """The sample point lies on an arrangement feature, or in no open face."""


@dataclass(frozen=True, slots=True)
class Crossing:
    """Intersection of the query line with one arrangement edge.

    Sorted along the query by ``t = parameter_on_line(query, point)``.
    ``vertex_id`` is set iff the intersection is an arrangement vertex.
    """

    point: Point2D
    t: Scalar
    edge_id: int
    vertex_id: int | None


@dataclass(frozen=True, slots=True)
class Zone:
    """Ordered features crossed by a query line that is not an arrangement line.

    ``faces[0]`` is the face containing a point of ``L`` before the first
    crossing. After that, for a transverse crossing ``faces[i+1]`` is the
    face on the other side of ``edges[i]``. Through a vertex, ``faces[i+1]``
    is the open face of ``L`` just after that vertex. On an overlapping
    query, ``faces`` are the incident faces of the coincident line in
    order of first appearance along ``L``. ``vertices`` lists vertex
    features in walk order.
    """

    query: Line2D
    faces: tuple[Face, ...]
    edges: tuple[EdgePiece, ...]
    vertices: tuple[Vertex, ...]
    crossings: tuple[Crossing, ...]

    @property
    def cells(self) -> tuple[Face, ...]:
        """Design name for ``faces``."""
        return self.faces


@dataclass(frozen=True, slots=True)
class SupportingLineZone:
    """Vertices of every face incident to an arrangement line ``L``.

    Split into vertices that lie on ``L`` (triples, later) versus opposite
    vertices (alignment candidates). This is not ``COMPUTE_ZONE``: a
    supporting line overlaps its own edges.
    """

    line_index: int
    faces: tuple[Face, ...]
    vertices_on_line: tuple[Vertex, ...]
    opposite_vertices: tuple[Vertex, ...]

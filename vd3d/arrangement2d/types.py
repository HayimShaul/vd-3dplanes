"""DCEL types for an unbounded 2D line arrangement."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar

from vd3d.arrangement2d.predicates import Direction

EdgeKind = str  # "segment" | "ray" | "line"


@dataclass(slots=True)
class Vertex:
    id: int
    point: Point2D
    line_indices: tuple[int, ...]
    outgoing: tuple[int, ...] = ()


@dataclass(slots=True)
class EdgePiece:
    """Undirected piece of a supporting line: a segment, a ray, or a whole line."""

    id: int
    kind: EdgeKind
    line_index: int
    start_vertex: int | None
    end_vertex: int | None
    t_min: Scalar | None
    t_max: Scalar | None


@dataclass(slots=True)
class HalfEdge:
    id: int
    origin_id: int | None
    twin_id: int
    edge_id: int
    line_index: int
    direction: Direction
    next_id: int | None = None
    prev_id: int | None = None
    face_id: int | None = None


@dataclass(slots=True)
class Face:
    id: int
    half_edge_id: int | None
    unbounded: bool
    representative: Point2D


@dataclass(slots=True)
class Arrangement2D:
    lines: tuple[Line2D, ...]
    vertices: tuple[Vertex, ...]
    edges: tuple[EdgePiece, ...]
    half_edges: tuple[HalfEdge, ...]
    faces: tuple[Face, ...]

    def cycle(self, face: Face) -> tuple[HalfEdge, ...]:
        """Boundary half-edges of ``face``, in walk order (face on the left)."""
        if face.half_edge_id is None:
            return ()
        start = self.half_edges[face.half_edge_id]
        out = [start]
        current = self.half_edges[start.next_id]  # type: ignore[index]
        guard = len(self.half_edges) + 1
        while current.id != start.id:
            out.append(current)
            if current.next_id is None:
                raise RuntimeError(f"half-edge {current.id} has no next")
            current = self.half_edges[current.next_id]
            if len(out) > guard:
                raise RuntimeError(f"face {face.id} boundary did not close")
        return tuple(out)

    @property
    def bounded_faces(self) -> tuple[Face, ...]:
        return tuple(face for face in self.faces if not face.unbounded)

    @property
    def unbounded_faces(self) -> tuple[Face, ...]:
        return tuple(face for face in self.faces if face.unbounded)

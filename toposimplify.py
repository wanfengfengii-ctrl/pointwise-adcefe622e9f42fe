"""Topology-preserving simplification of shared region boundaries.

Implements simplify(topology, tolerance) from TOPOLOGY.md: vertices are
removed from shared arcs built by Archive.build_topology; every removal
applies to all users of the arc simultaneously and is accepted only when
the cumulative tolerance and every topology invariant still hold.  All
arithmetic is exact Fraction arithmetic; each call restarts from the
build-time topology and never mutates its inputs.
"""
from fractions import Fraction
from functools import cmp_to_key

from polyarchive import (Region, Topology, _SEAL, _segment_kind, _on,
                         _interiors_overlap)


def _segments(points, closed):
    segs = list(zip(points, points[1:]))
    if closed:
        segs.append((points[-1], points[0]))
    return segs


def _within_tolerance(p, a, b, tol2):
    """Squared distance from p to the finite segment a-b is at most tol2."""
    dx, dy = b[0] - a[0], b[1] - a[1]
    px, py = p[0] - a[0], p[1] - a[1]
    length2 = dx * dx + dy * dy
    if length2 == 0:
        dist2 = px * px + py * py
    else:
        dot = px * dx + py * dy
        if dot <= 0:
            dist2 = px * px + py * py
        elif dot >= length2:
            qx, qy = p[0] - b[0], p[1] - b[1]
            dist2 = qx * qx + qy * qy
        else:
            dist2 = px * px + py * py - dot * dot / length2
    return dist2 <= tol2


def _arc_simple(points, closed):
    """The polyline does not cross, touch or backtrack over itself."""
    if len(points) < 2 or any(a == b for a, b in zip(points, points[1:])):
        return False
    segs = _segments(points, closed)
    count = len(segs)
    loop = not closed and points[0] == points[-1]
    for i in range(count):
        for j in range(i + 1, count):
            adjacent = (j == i + 1
                        or (i == 0 and j == count - 1 and (closed or loop)))
            kind = _segment_kind(segs[i][0], segs[i][1],
                                 segs[j][0], segs[j][1])
            if kind != ("point" if adjacent else "none"):
                return False
    return True


def _arcs_meet_ok(points_a, closed_a, points_b, closed_b):
    """Two arcs may only meet at shared build-time endpoint junctions."""
    ends_a = set() if closed_a else {points_a[0], points_a[-1]}
    ends_b = set() if closed_b else {points_b[0], points_b[-1]}
    common = ends_a & ends_b
    for s in _segments(points_a, closed_a):
        for t in _segments(points_b, closed_b):
            kind = _segment_kind(s[0], s[1], t[0], t[1])
            if kind == "none":
                continue
            if kind != "point":
                return False
            contacts = {p for p in s + t
                        if _on(p, s[0], s[1]) and _on(p, t[0], t[1])}
            if not contacts <= common:
                return False
    return True


def _angle_cmp(u, v):
    half_u = 0 if (u[1] > 0 or (u[1] == 0 and u[0] > 0)) else 1
    half_v = 0 if (v[1] > 0 or (v[1] == 0 and v[0] > 0)) else 1
    if half_u != half_v:
        return -1 if half_u < half_v else 1
    cross = u[0] * v[1] - u[1] * v[0]
    if cross == 0:
        return 0
    return -1 if cross > 0 else 1


def _junction_order(ends, points_of):
    """Arc-end labels around one junction in counterclockwise order."""
    vecs = []
    for arc_index, end in ends:
        pts = points_of(arc_index)
        base, tip = (pts[0], pts[1]) if end == 0 else (pts[-1], pts[-2])
        vecs.append(((arc_index, end), (tip[0] - base[0], tip[1] - base[1])))
    vecs.sort(key=cmp_to_key(lambda x, y: _angle_cmp(x[1], y[1])))
    for i in range(len(vecs)):
        ux, uy = vecs[i][1]
        vx, vy = vecs[(i + 1) % len(vecs)][1]
        if ux * vy - uy * vx == 0 and ux * vx + uy * vy > 0:
            return None
    return tuple(label for label, _ in vecs)


def _same_cycle(first, second):
    if first is None or second is None or len(first) != len(second):
        return False
    count = len(first)
    return any(all(first[(s + i) % count] == second[i] for i in range(count))
               for s in range(count))


def simplify(topology, tolerance):
    """Remove vertices from shared arcs; see TOPOLOGY.md for the contract."""
    if (not isinstance(topology, Topology)
            or getattr(topology, "_seal", None) is not _SEAL):
        raise ValueError("topology must come from Archive.build_topology")
    if not isinstance(tolerance, Fraction) or tolerance < 0:
        raise ValueError("tolerance must be a nonnegative Fraction")
    tol2 = tolerance * tolerance
    arcs = topology.arcs
    closed = [arc.closed for arc in arcs]
    kept = [list(range(len(arc.points))) for arc in arcs]

    def current(index):
        arc = arcs[index]
        return tuple(arc.points[i] for i in kept[index])

    ends_at = {}
    for index, arc in enumerate(arcs):
        if not arc.closed:
            ends_at.setdefault(arc.points[0], []).append((index, 0))
            ends_at.setdefault(arc.points[-1], []).append((index, 1))
    base_order = {junction: _junction_order(ends, current)
                  for junction, ends in ends_at.items()}

    def expand(refs, points_of):
        if len(refs) == 1:
            index, direction = refs[0]
            pts = points_of(index)
            if direction == -1:
                pts = pts[::-1]
            return pts if closed[index] else pts[:-1]
        out = []
        for index, direction in refs:
            pts = points_of(index)
            if direction == -1:
                pts = pts[::-1]
            out.extend(pts[1:] if out else pts)
        return tuple(out[:-1])

    def build_region(region_id, points_of):
        entries = topology.rings[region_id]
        return Region(expand(entries[0], points_of),
                      tuple(expand(refs, points_of) for refs in entries[1:]))

    regions = {rid: build_region(rid, current) for rid in topology.rings}

    def acceptable(index, pos):
        chain = kept[index]
        original = arcs[index].points
        lo, hi = chain[pos - 1], chain[(pos + 1) % len(chain)]
        a, b = original[lo], original[hi]
        cursor = lo
        while True:  # cumulative: every original vertex of the subchain
            if not _within_tolerance(original[cursor], a, b, tol2):
                return False
            if cursor == hi:
                break
            cursor = (cursor + 1) % len(original)
        trial = tuple(p for i, p in enumerate(current(index)) if i != pos)
        if not _arc_simple(trial, closed[index]):
            return False

        def points_of(other):
            return trial if other == index else current(other)

        for other in range(len(arcs)):
            if other != index and not _arcs_meet_ok(trial, closed[index],
                                                    current(other),
                                                    closed[other]):
                return False
        if not closed[index]:
            for junction in {trial[0], trial[-1]}:
                if not _same_cycle(
                        _junction_order(ends_at[junction], points_of),
                        base_order[junction]):
                    return False
        rebuilt = {}
        for rid in sorted({user[0] for user in arcs[index].users}):
            try:
                rebuilt[rid] = build_region(rid, points_of)
            except ValueError:
                return False
        for rid, region in rebuilt.items():
            for sid in topology.rings:
                if sid != rid and _interiors_overlap(
                        region, rebuilt.get(sid, regions[sid])):
                    return False
        return True

    while True:
        candidates = []
        for index, arc in enumerate(arcs):
            chain = kept[index]
            positions = (range(len(chain)) if arc.closed
                         else range(1, len(chain) - 1))
            candidates.extend((arc.points[chain[pos]], index, pos)
                              for pos in positions)
        candidates.sort(key=lambda item: item[0])
        for _, index, pos in candidates:
            if acceptable(index, pos):
                del kept[index][pos]
                for rid in {user[0] for user in arcs[index].users}:
                    regions[rid] = build_region(rid, current)
                break
        else:
            break
    return {rid: regions[rid] for rid in topology.rings}

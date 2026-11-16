"""Exact independent polygon storage; Python 3.11 standard library."""
from dataclasses import dataclass
from fractions import Fraction
import json

Point = tuple[Fraction, Fraction]
Ring = tuple[Point, ...]


def _sequence(value):
    if not isinstance(value, (list, tuple)):
        raise ValueError("expected list or tuple")
    return value


def _point(value):
    _sequence(value)
    if len(value) != 2 or not all(isinstance(x, Fraction) for x in value):
        raise ValueError("coordinates must be Fraction pairs")
    return tuple(value)


def _cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on(p, a, b):
    return (_cross(a, b, p) == 0
            and min(a[0], b[0]) <= p[0] <= max(a[0], b[0])
            and min(a[1], b[1]) <= p[1] <= max(a[1], b[1]))


def _edges(ring):
    return zip(ring, ring[1:] + ring[:1])


def _segment_kind(a, b, c, d):
    ab_c, ab_d = _cross(a, b, c), _cross(a, b, d)
    cd_a, cd_b = _cross(c, d, a), _cross(c, d, b)
    if ab_c == ab_d == cd_a == cd_b == 0:
        common = {p for p in (a, b, c, d) if _on(p, a, b) and _on(p, c, d)}
        return "overlap" if len(common) > 1 else "point" if common else "none"
    if ab_c * ab_d < 0 and cd_a * cd_b < 0:
        return "cross"
    if (_on(c, a, b) or _on(d, a, b)
            or _on(a, c, d) or _on(b, c, d)):
        return "point"
    return "none"


def _twice_area(ring):
    return sum((a[0] * b[1] - b[0] * a[1] for a, b in _edges(ring)), Fraction(0))


def _location(ring, p):
    inside = False
    for a, b in _edges(ring):
        if _on(p, a, b):
            return "boundary"
        if (a[1] > p[1]) != (b[1] > p[1]):
            x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if p[0] < x:
                inside = not inside
    return "interior" if inside else "exterior"


def _ring(value):
    ring = tuple(_point(p) for p in _sequence(value))
    if len(ring) < 3 or any(a == b for a, b in _edges(ring)):
        raise ValueError("ring needs at least three distinct consecutive vertices")
    if _twice_area(ring) == 0:
        raise ValueError("zero-area ring")
    edges = list(_edges(ring))
    for i, (a, b) in enumerate(edges):
        for j in range(i + 1, len(edges)):
            c, d = edges[j]
            kind = _segment_kind(a, b, c, d)
            adjacent = j == i + 1 or (i == 0 and j == len(edges) - 1)
            if (adjacent and kind != "point") or (not adjacent and kind != "none"):
                raise ValueError("ring is not simple")
    return ring


def _rings_meet(a, b):
    return any(_segment_kind(p, q, r, s) != "none"
               for p, q in _edges(a) for r, s in _edges(b))


@dataclass(frozen=True)
class Region:
    outer: Ring
    holes: tuple[Ring, ...] = ()

    def __post_init__(self):
        outer = _ring(self.outer)
        holes = tuple(_ring(h) for h in _sequence(self.holes))
        for h in holes:
            if _rings_meet(outer, h) or _location(outer, h[0]) != "interior":
                raise ValueError("hole must be strictly inside outer ring")
        for i, a in enumerate(holes):
            for b in holes[i + 1:]:
                if (_rings_meet(a, b) or _location(a, b[0]) != "exterior"
                        or _location(b, a[0]) != "exterior"):
                    raise ValueError("holes must have disjoint closures")
        object.__setattr__(self, "outer", outer)
        object.__setattr__(self, "holes", holes)

    @property
    def area(self):
        return (abs(_twice_area(self.outer))
                - sum((abs(_twice_area(h)) for h in self.holes), Fraction(0))) / 2

    def locate(self, point):
        point = _point(point)
        result = _location(self.outer, point)
        if result != "interior":
            return result
        for hole in self.holes:
            result = _location(hole, point)
            if result == "boundary":
                return result
            if result == "interior":
                return "exterior"
        return "interior"


def _identifier(value):
    if not isinstance(value, str) or not value:
        raise ValueError("region ID must be a nonempty string")
    return value


def _keys(value, expected):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError("unexpected JSON object fields")


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError("nonfinite JSON number: " + value)


def _decode_ring(value):
    if not isinstance(value, list):
        raise ValueError("wire ring must be an array")
    points = []
    for pair in value:
        if not isinstance(pair, list) or len(pair) != 2:
            raise ValueError("wire point must be a two-element array")
        coords = []
        for text in pair:
            if not isinstance(text, str):
                raise ValueError("wire coordinate must be a canonical rational string")
            try:
                number = Fraction(text)
            except (ValueError, ZeroDivisionError) as exc:
                raise ValueError("invalid rational") from exc
            if str(number) != text:
                raise ValueError("noncanonical rational")
            coords.append(number)
        points.append(tuple(coords))
    return tuple(points)


def _encode_ring(ring):
    return [[str(x), str(y)] for x, y in ring]


class Archive:
    """Independent region IDs; put validates before replacing any stored value."""

    def __init__(self):
        self._regions = {}

    @property
    def ids(self):
        return tuple(sorted(self._regions))

    def put(self, region_id, outer, holes=()):
        region_id = _identifier(region_id)
        candidate = Region(outer, holes)
        self._regions[region_id] = candidate

    def read(self, region_id):
        return self._regions[_identifier(region_id)]

    def export(self, region_ids=None):
        if region_ids is None:
            selected = self.ids
        else:
            selected = tuple(_identifier(x) for x in _sequence(region_ids))
            if len(set(selected)) != len(selected):
                raise ValueError("duplicate selected region ID")
            selected = tuple(sorted(selected))
        records = []
        for region_id in selected:
            region = self.read(region_id)
            records.append({"id": region_id, "outer": _encode_ring(region.outer),
                            "holes": [_encode_ring(h) for h in region.holes]})
        return {"version": 1, "regions": records}

    def dumps(self, region_ids=None):
        return json.dumps(self.export(region_ids), ensure_ascii=True, sort_keys=True,
                          separators=(",", ":")) + "\n"

    def build_topology(self, region_ids=None):
        """Unified shared-boundary representation of the selected regions.

        Selection follows export(); the archive is never modified.  See
        TOPOLOGY.md for the contract.
        """
        if region_ids is None:
            selected = self.ids
        else:
            selected = tuple(_identifier(x) for x in _sequence(region_ids))
            if len(set(selected)) != len(selected):
                raise ValueError("duplicate selected region ID")
            selected = tuple(sorted(selected))
        regions = [(region_id, self.read(region_id)) for region_id in selected]
        edges = sum(len(region.outer) + sum(map(len, region.holes))
                    for _, region in regions)
        if edges > _MAX_TOPOLOGY_EDGES:
            raise ValueError("too many input edges for topology")
        for i, (_, first) in enumerate(regions):
            for _, second in regions[i + 1:]:
                if _interiors_overlap(first, second):
                    raise ValueError("region interiors overlap")
        arcs, rings = _topology_data(regions)
        return Topology(arcs, rings, _seal=_SEAL)

    @classmethod
    def loads(cls, text):
        if not isinstance(text, str):
            raise ValueError("archive input must be text")
        try:
            data = json.loads(text, object_pairs_hook=_object, parse_constant=_reject_constant)
        except (ValueError, RecursionError) as exc:
            raise ValueError("invalid archive JSON") from exc
        _keys(data, ("version", "regions"))
        if type(data["version"]) is not int or data["version"] != 1:
            raise ValueError("unsupported archive version")
        if not isinstance(data["regions"], list):
            raise ValueError("regions must be an array")
        archive = cls()
        for record in data["regions"]:
            _keys(record, ("id", "outer", "holes"))
            region_id = _identifier(record["id"])
            if region_id in archive._regions:
                raise ValueError("duplicate region ID")
            if not isinstance(record["holes"], list):
                raise ValueError("holes must be an array")
            archive.put(region_id, _decode_ring(record["outer"]),
                        tuple(_decode_ring(h) for h in record["holes"]))
        return archive


# --- Shared-boundary topology (TOPOLOGY.md) ---------------------------------

_MAX_TOPOLOGY_EDGES = 120
_SEAL = object()


@dataclass(frozen=True)
class Arc:
    """Maximal boundary chain between junctions; shared geometry stored once.

    points: Fraction coordinate pairs in storage direction.  closed is True
    only for endpoint-less independent closed arcs (points does not repeat
    the first point; last-to-first is also an edge).  users is a tuple of
    (region_id, ring_index, direction) with direction +1/-1, one per user.
    """
    points: tuple
    closed: bool
    users: tuple


class Topology:
    """Arc table plus per-region ring references built by build_topology."""

    __slots__ = ("arcs", "rings", "_seal")

    def __init__(self, arcs, rings, _seal=None):
        if _seal is not _SEAL:
            raise ValueError("Topology must come from Archive.build_topology")
        self.arcs = arcs
        self.rings = rings
        self._seal = _seal


def _seg_dist2(p, a, b):
    """Exact squared distance from p to the finite segment a-b."""
    dx, dy = b[0] - a[0], b[1] - a[1]
    px, py = p[0] - a[0], p[1] - a[1]
    length2 = dx * dx + dy * dy
    if length2 == 0:
        return px * px + py * py
    dot = px * dx + py * dy
    if dot <= 0:
        return px * px + py * py
    if dot >= length2:
        qx, qy = p[0] - b[0], p[1] - b[1]
        return qx * qx + qy * qy
    return px * px + py * py - dot * dot / length2


def _interior_point(region):
    """An exact point strictly inside the region's material."""
    ring = region.outer
    a, b = ring[0], ring[1]
    dx, dy = b[0] - a[0], b[1] - a[1]
    nx, ny = (-dy, dx) if _twice_area(ring) > 0 else (dy, -dx)
    m = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
    best = None
    for boundary in (ring,) + region.holes:
        for c, d in _edges(boundary):
            if (c, d) == (a, b):
                continue
            dist = _seg_dist2(m, c, d)
            if best is None or dist < best:
                best = dist
    length2 = nx * nx + ny * ny
    step = min(best, length2) / (2 * length2)
    return (m[0] + nx * step, m[1] + ny * step)


def _interiors_overlap(first, second):
    """True iff the two regions share any positive-area interior."""
    for ring_a in (first.outer,) + first.holes:
        for ring_b in (second.outer,) + second.holes:
            for a, b in _edges(ring_a):
                for c, d in _edges(ring_b):
                    if _segment_kind(a, b, c, d) == "cross":
                        return True
    return (second.locate(_interior_point(first)) == "interior"
            or first.locate(_interior_point(second)) == "interior")


def _topology_data(regions):
    """Build (arcs, rings) for selected (region_id, Region) pairs.

    Edges are split at every input vertex so equal geometric atoms exist
    exactly once; atoms are chained into maximal arcs at junctions (degree
    != 2 or changing user ring-label sets).
    """
    rings = []
    for region_id, region in regions:
        rings.append(((region_id, 0), region.outer))
        rings.extend(((region_id, i + 1), hole)
                     for i, hole in enumerate(region.holes))
    points = set()
    for _, ring in rings:
        points.update(ring)
    atoms = {}
    sequences = {}
    for label, ring in rings:
        sequence = []
        for a, b in _edges(ring):
            dx, dy = b[0] - a[0], b[1] - a[1]
            on = sorted((p for p in points if _on(p, a, b)),
                        key=lambda p: (p[0] - a[0]) * dx + (p[1] - a[1]) * dy)
            for u, w in zip(on, on[1:]):
                key = (u, w) if u <= w else (w, u)
                atoms.setdefault(key, {})[label] = 1 if u <= w else -1
                sequence.append(key)
        sequences[label] = sequence
    adjacency = {}
    for u, w in atoms:
        adjacency.setdefault(u, []).append((u, w))
        adjacency.setdefault(w, []).append((u, w))
    for keys in adjacency.values():
        keys.sort()
    junctions = {p for p, keys in adjacency.items()
                 if len(keys) != 2
                 or set(atoms[keys[0]]) != set(atoms[keys[1]])}
    arcs = []
    atom_arc = {}
    visited = set()

    def register(pts, closed, dirs):
        index = len(arcs)
        arcs.append((pts, closed, dirs))
        for key, sign in dirs.items():
            atom_arc[key] = (index, sign)

    for junction in sorted(junctions):
        for key in adjacency[junction]:
            if key in visited:
                continue
            other = key[1] if key[0] == junction else key[0]
            pts = [junction, other]
            dirs = {key: 1 if key == (junction, other) else -1}
            visited.add(key)
            cur, prev = other, key
            while cur not in junctions:
                first, second = adjacency[cur]
                nxt = first if second == prev else second
                other = nxt[1] if nxt[0] == cur else nxt[0]
                dirs[nxt] = 1 if nxt == (cur, other) else -1
                visited.add(nxt)
                pts.append(other)
                prev, cur = nxt, other
            register(tuple(pts), False, dirs)
    for start in sorted(atoms):
        if start in visited:
            continue
        u, w = start
        pts = [u, w]
        dirs = {start: 1}
        visited.add(start)
        cur, prev = w, start
        while cur != u:
            first, second = adjacency[cur]
            nxt = first if second == prev else second
            other = nxt[1] if nxt[0] == cur else nxt[0]
            dirs[nxt] = 1 if nxt == (cur, other) else -1
            visited.add(nxt)
            pts.append(other)
            prev, cur = nxt, other
        pts.pop()
        register(tuple(pts), True, dirs)
    arc_objs = []
    for pts, closed, dirs in arcs:
        key0 = next(iter(dirs))
        forward = dirs[key0]
        users = tuple(sorted(
            (label[0], label[1], sign * forward)
            for label, sign in atoms[key0].items()))
        arc_objs.append(Arc(pts, closed, users))
    for arc, (_, _, dirs) in zip(arc_objs, arcs):
        expect = {(u[0], u[1]): u[2] for u in arc.users}
        for key, forward in dirs.items():
            assert {label: sign * forward
                    for label, sign in atoms[key].items()} == expect
    rings_map = {}
    for label, _ring in rings:
        region_id, _ = label
        refs = []
        for key in sequences[label]:
            index, forward = atom_arc[key]
            direction = atoms[key][label] * forward
            if refs and refs[-1][0] == index:
                assert refs[-1][1] == direction
            else:
                refs.append((index, direction))
        if len(refs) > 1 and refs[0][0] == refs[-1][0]:
            assert refs[0][1] == refs[-1][1]
            refs = refs[1:-1] + [refs[-1]]
        rings_map.setdefault(region_id, []).append(tuple(refs))
    rings_map = {region_id: tuple(entries)
                 for region_id, entries in rings_map.items()}
    for index, arc in enumerate(arc_objs):
        expect = {(u[0], u[1]): u[2] for u in arc.users}
        got = {}
        for region_id, entries in rings_map.items():
            for ring_index, refs in enumerate(entries):
                for arc_index, direction in refs:
                    if arc_index == index:
                        assert (region_id, ring_index) not in got
                        got[(region_id, ring_index)] = direction
        assert got == expect
    return tuple(arc_objs), rings_map

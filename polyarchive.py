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

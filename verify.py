"""Independent exact reference verification for the topology pipeline.

Runs the frozen baseline suite, then re-checks Archive.build_topology and
toposimplify.simplify against independently recomputed oracles:

  * public representation expansion (arcs + ring references reproduce the
    input rings up to cyclic start and inserted collinear split points),
  * shared geometry identity (stored once; users and references consistent;
    every user of a shared arc sees the identical simplified chain),
  * cumulative tolerance (every replacement segment stays within tolerance
    of every original vertex of its build-time subchain, recomputed here),
  * cross-region contacts (no crossings, no new contacts, junctions and
    their cyclic arc-end orders preserved, interiors still disjoint).

Nothing below reuses the simplifier's own accept/reject decisions; all
geometry is recomputed from the public topology attributes and the
returned regions with exact Fraction arithmetic.
"""
import subprocess
import sys
from fractions import Fraction as F
from functools import cmp_to_key
from pathlib import Path

from polyarchive import Archive, Region
from toposimplify import simplify


# --- independent exact geometry ---------------------------------------------

def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_seg(p, a, b):
    return (cross(a, b, p) == 0
            and min(a[0], b[0]) <= p[0] <= max(a[0], b[0])
            and min(a[1], b[1]) <= p[1] <= max(a[1], b[1]))


def seg_kind(a, b, c, d):
    abc, abd = cross(a, b, c), cross(a, b, d)
    cda, cdb = cross(c, d, a), cross(c, d, b)
    if abc == abd == cda == cdb == 0:
        common = {p for p in (a, b, c, d) if on_seg(p, a, b) and on_seg(p, c, d)}
        return "overlap" if len(common) > 1 else "point" if common else "none"
    if abc * abd < 0 and cda * cdb < 0:
        return "cross"
    if (on_seg(c, a, b) or on_seg(d, a, b)
            or on_seg(a, c, d) or on_seg(b, c, d)):
        return "point"
    return "none"


def edges(ring):
    return list(zip(ring, ring[1:] + ring[:1]))


def twice_area(ring):
    return sum((a[0] * b[1] - b[0] * a[1] for a, b in edges(ring)), F(0))


def dist2_point_seg(p, a, b):
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


def locate(ring, p):
    inside = False
    for a, b in edges(ring):
        if on_seg(p, a, b):
            return "boundary"
        if (a[1] > p[1]) != (b[1] > p[1]):
            x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if p[0] < x:
                inside = not inside
    return "interior" if inside else "exterior"


def region_locate(region, p):
    result = locate(region.outer, p)
    if result != "interior":
        return result
    for hole in region.holes:
        result = locate(hole, p)
        if result == "boundary":
            return result
        if result == "interior":
            return "exterior"
    return "interior"


def interior_point(region):
    """Independently recomputed strict interior point of a region."""
    ring = region.outer
    a, b = ring[0], ring[1]
    dx, dy = b[0] - a[0], b[1] - a[1]
    nx, ny = (-dy, dx) if twice_area(ring) > 0 else (dy, -dx)
    m = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
    best = None
    for boundary in (ring,) + region.holes:
        for c, d in edges(boundary):
            if (c, d) == (a, b):
                continue
            dist = dist2_point_seg(m, c, d)
            if best is None or dist < best:
                best = dist
    length2 = nx * nx + ny * ny
    step = min(best, length2) / (2 * length2)
    return (m[0] + nx * step, m[1] + ny * step)


def interiors_overlap(first, second):
    for ring_a in (first.outer,) + first.holes:
        for ring_b in (second.outer,) + second.holes:
            for e1 in edges(ring_a):
                for e2 in edges(ring_b):
                    if seg_kind(*e1, *e2) == "cross":
                        return True
    return (region_locate(second, interior_point(first)) == "interior"
            or region_locate(first, interior_point(second)) == "interior")


def angle_cmp(u, v):
    half_u = 0 if (u[1] > 0 or (u[1] == 0 and u[0] > 0)) else 1
    half_v = 0 if (v[1] > 0 or (v[1] == 0 and v[0] > 0)) else 1
    if half_u != half_v:
        return -1 if half_u < half_v else 1
    c = u[0] * v[1] - u[1] * v[0]
    if c == 0:
        return 0
    return -1 if c > 0 else 1


def cyclic_order(ends, points_of):
    vecs = []
    for index, end in ends:
        pts = points_of(index)
        a, b = (pts[0], pts[1]) if end == 0 else (pts[-1], pts[-2])
        vecs.append(((index, end), (b[0] - a[0], b[1] - a[1])))
    vecs.sort(key=cmp_to_key(lambda x, y: angle_cmp(x[1], y[1])))
    return tuple(label for label, _ in vecs)


def same_cycle(first, second):
    n = len(first)
    return len(first) == len(second) and any(
        all(first[(s + i) % n] == second[i] for i in range(n))
        for s in range(n))


# --- independent expansion of the public topology ----------------------------

def expand_ring(arcs, refs):
    if len(refs) == 1:
        index, direction = refs[0]
        pts = arcs[index].points
        if direction == -1:
            pts = pts[::-1]
        return tuple(pts) if arcs[index].closed else tuple(pts[:-1])
    out = []
    for index, direction in refs:
        pts = arcs[index].points
        if direction == -1:
            pts = pts[::-1]
        out.extend(pts[1:] if out else pts)
    return tuple(out[:-1])


def expand_topology(t):
    return {rid: [expand_ring(t.arcs, refs) for refs in entries]
            for rid, entries in t.rings.items()}


def matches_input(expanded, source):
    """expanded == source ring plus collinear split points, same direction."""
    n = len(expanded)
    for start in (i for i, p in enumerate(expanded) if p == source[0]):
        pos, ok = start, True
        for k in range(len(source)):
            a, b = source[k], source[(k + 1) % len(source)]
            if expanded[pos % n] != a:
                ok = False
                break
            dx, dy = b[0] - a[0], b[1] - a[1]
            limit = dx * dx + dy * dy
            prev = F(0)
            pos += 1
            while expanded[pos % n] != b:
                p = expanded[pos % n]
                proj = (p[0] - a[0]) * dx + (p[1] - a[1]) * dy
                if not (prev < proj < limit) or cross(a, b, p) != 0:
                    ok = False
                    break
                prev = proj
                pos += 1
                if pos - start > n:
                    ok = False
                    break
            if not ok:
                break
        if ok and pos % n == start:
            return True
    return False


def arc_segments(arc):
    segs = list(zip(arc.points, arc.points[1:]))
    if arc.closed:
        segs.append((arc.points[-1], arc.points[0]))
    return segs


def canon_cycle(seq):
    best = min(range(len(seq)), key=lambda i: seq[i:] + seq[:i])
    return seq[best:] + seq[:best]


def arc_chain(t, out, index, user):
    """Current chain of one arc, in storage direction, read off a user ring."""
    arc = t.arcs[index]
    rid, ring_index, direction = user
    region = out[rid]
    ring = region.outer if ring_index == 0 else region.holes[ring_index - 1]
    if arc.closed:
        cyc = ring if direction == 1 else tuple(reversed(ring))
        return canon_cycle(tuple(cyc))
    j1, j2 = arc.points[0], arc.points[-1]
    start, stop = (j1, j2) if direction == 1 else (j2, j1)
    n = len(ring)
    i = ring.index(start)
    if start == stop:
        # loop arc: the user ring is exactly this arc
        chain = [ring[(i + k) % n] for k in range(n)] + [start]
    else:
        chain = []
        while True:
            chain.append(ring[i])
            if ring[i] == stop:
                break
            i = (i + 1) % n
    if direction == -1:
        chain.reverse()
    return tuple(chain)


# --- generic pipeline verification -------------------------------------------

def check_pipeline(archive, ids, tolerance, label):
    before = archive.dumps()
    t = archive.build_topology(ids)
    snapshot = t.arcs
    expanded = expand_topology(t)
    for rid in ids:
        source = archive.read(rid)
        assert matches_input(expanded[rid][0], source.outer), (label, rid, "outer")
        for i, hole in enumerate(source.holes):
            assert matches_input(expanded[rid][i + 1], hole), (label, rid, i)
    # shared geometry stored once: arcs meet only at common endpoints
    for i, arc in enumerate(t.arcs):
        ends_i = set() if arc.closed else {arc.points[0], arc.points[-1]}
        for j in range(i + 1, len(t.arcs)):
            other = t.arcs[j]
            ends_j = (set() if other.closed
                      else {other.points[0], other.points[-1]})
            for s in arc_segments(arc):
                for u in arc_segments(other):
                    kind = seg_kind(*s, *u)
                    assert kind != "overlap", (label, "duplicated geometry", i, j)
                    assert kind != "cross", (label, "arc crossing", i, j)
                    if kind == "point":
                        contacts = {p for p in s + u
                                    if on_seg(p, *s) and on_seg(p, *u)}
                        assert contacts <= (ends_i & ends_j), (label, i, j)
    # users <-> ring references consistency
    for index, arc in enumerate(t.arcs):
        expect = {(u[0], u[1]): u[2] for u in arc.users}
        got = {}
        for rid, entries in t.rings.items():
            for ring_index, refs in enumerate(entries):
                for arc_index, direction in refs:
                    if arc_index == index:
                        assert (rid, ring_index) not in got
                        got[(rid, ring_index)] = direction
        assert got == expect, (label, "users/refs", index)
    out = simplify(t, tolerance)
    assert archive.dumps() == before, (label, "archive mutated")
    assert t.arcs == snapshot, (label, "topology mutated")
    assert set(out) == set(ids), (label, "output ids")
    assert all(isinstance(region, Region) for region in out.values())
    check_output(t, archive, ids, out, tolerance, label)
    return t, out


def check_output(t, archive, ids, out, tolerance, label):
    tol2 = tolerance * tolerance
    for rid in ids:
        source, result = archive.read(rid), out[rid]
        assert len(result.holes) == len(source.holes), (label, rid, "holes")
        for got, want in zip((result.outer,) + result.holes,
                             (source.outer,) + source.holes):
            assert twice_area(got) * twice_area(want) > 0, (label, rid, "direction")
    vertices = set()
    for arc in t.arcs:
        vertices.update(arc.points)
    for rid in ids:
        for ring in (out[rid].outer,) + out[rid].holes:
            for p in ring:
                assert p in vertices, (label, rid, "foreign vertex", p)
    # every user of a shared arc sees the identical current chain
    chains = {}
    for index, arc in enumerate(t.arcs):
        seen = [arc_chain(t, out, index, user) for user in arc.users]
        assert all(chain == seen[0] for chain in seen), (label, "diverged", index)
        chains[index] = seen[0]
        if not arc.closed:
            assert seen[0][0] == arc.points[0], (label, "start moved", index)
            assert seen[0][-1] == arc.points[-1], (label, "end moved", index)
    # cumulative tolerance over build-time subchains, recomputed
    for index, arc in enumerate(t.arcs):
        original = arc.points
        chain = chains[index]
        if arc.closed:
            start = original.index(chain[0])
            original = original[start:] + original[:start]
            segs = list(zip(chain, chain[1:])) + [(chain[-1], chain[0])]
        else:
            segs = list(zip(chain, chain[1:]))
        for u, w in segs:
            i = original.index(u)
            stop = original.index(w)
            while True:
                p = original[i]
                assert dist2_point_seg(p, u, w) <= tol2, (label, "tolerance", index, p)
                if i == stop:
                    break
                i = (i + 1) % len(original)
    # junctions and their cyclic arc-end orders preserved
    junctions = set()
    ends_at = {}
    for index, arc in enumerate(t.arcs):
        if arc.closed:
            continue
        junctions.add(arc.points[0])
        junctions.add(arc.points[-1])
        ends_at.setdefault(arc.points[0], []).append((index, 0))
        ends_at.setdefault(arc.points[-1], []).append((index, 1))
    for junction, ends in ends_at.items():
        before_order = cyclic_order(ends, lambda i: t.arcs[i].points)
        after_order = cyclic_order(ends, lambda i: chains[i])
        assert same_cycle(before_order, after_order), (label, "order", junction)
    # cross-region contacts: no crossings, no new contacts, no unshared overlap
    shared = [(index, chains[index])
              for index, arc in enumerate(t.arcs) if len(arc.users) > 1]

    def on_shared_chain(rid1, rid2, p):
        for index, chain in shared:
            if not {rid1, rid2} <= {u[0] for u in t.arcs[index].users}:
                continue
            segs = list(zip(chain, chain[1:]))
            if t.arcs[index].closed:
                segs.append((chain[-1], chain[0]))
            if any(on_seg(p, *s) for s in segs):
                return True
        return False

    for x, rid1 in enumerate(ids):
        for rid2 in ids[x + 1:]:
            for ring1 in (out[rid1].outer,) + out[rid1].holes:
                for ring2 in (out[rid2].outer,) + out[rid2].holes:
                    for e1 in edges(ring1):
                        for e2 in edges(ring2):
                            kind = seg_kind(*e1, *e2)
                            assert kind != "cross", (label, rid1, rid2)
                            if kind == "point":
                                contacts = {p for p in e1 + e2
                                            if on_seg(p, *e1) and on_seg(p, *e2)}
                                for p in contacts:
                                    assert (p in junctions
                                            or on_shared_chain(rid1, rid2, p)), (
                                        label, "new contact", rid1, rid2, p)
                            elif kind == "overlap":
                                check_shared_overlap(
                                    t, shared, rid1, rid2, e1, e2, label)
    for x, rid1 in enumerate(ids):
        for rid2 in ids[x + 1:]:
            assert not interiors_overlap(out[rid1], out[rid2]), (
                label, "interiors overlap", rid1, rid2)


def check_shared_overlap(t, shared, rid1, rid2, e1, e2, label):
    a, b = e1
    dx, dy = b[0] - a[0], b[1] - a[1]
    common = sorted(
        ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy, p)
        for p in e1 + e2 if on_seg(p, *e1) and on_seg(p, *e2))
    lo, hi = common[0][1], common[-1][1]
    mid = ((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2)
    for index, chain in shared:
        if not {rid1, rid2} <= {u[0] for u in t.arcs[index].users}:
            continue
        segs = list(zip(chain, chain[1:]))
        if t.arcs[index].closed:
            segs.append((chain[-1], chain[0]))
        if any(on_seg(mid, *s) for s in segs):
            return
    raise AssertionError((label, "unshared overlap", rid1, rid2, e1, e2))


# --- scenarios ----------------------------------------------------------------

def ring(coords):
    return tuple((F(x), F(y)) for x, y in coords)


def box(x0, y0, x1, y1):
    return ring([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def scenario_straight_shared_edge():
    a = Archive()
    a.put("west", box(0, 0, 2, 4))
    a.put("east", ring([(2, 0), (5, 0), (5, 4), (2, 4), (2, 2)]))
    t, out = check_pipeline(a, ["west", "east"], F(0), "straight")
    assert len(t.arcs) == 3
    shared = [arc for arc in t.arcs if len(arc.users) == 2][0]
    assert sorted(u[2] for u in shared.users) == [-1, 1]
    assert (F(2), F(2)) not in out["west"].outer
    assert len(out["west"].outer) == 4 and len(out["east"].outer) == 4
    assert out["west"].area == 8 and out["east"].area == 12
    print("ok: straight shared edge, different segmentation (3 arcs, "
          "midpoint removed at tolerance 0, both sides identical)")


def scenario_hole_island():
    b = Archive()
    b.put("A", ring([(0, 0), (6, 0), (6, 6), (0, 6)]),
          (ring([(2, 2), (3, 2), (4, 2), (4, 4), (2, 4)]),))
    b.put("B", ring([(2, 2), (4, 2), (4, 4), (2, 4)]))
    t, out = check_pipeline(b, ["A", "B"], F(0), "hole-island")
    assert len(t.arcs) == 2 and all(arc.closed for arc in t.arcs)
    shared = [arc for arc in t.arcs if len(arc.users) == 2][0]
    assert shared.users[0][2] == shared.users[1][2]
    assert len(out["A"].holes[0]) == 4 and len(out["B"].outer) == 4
    assert out["A"].area == 32 and out["B"].area == 4
    assert out["A"].locate((F(3), F(3))) == "exterior"
    assert out["B"].locate((F(3), F(3))) == "interior"
    print("ok: hole-island shared closed arc (2 closed arcs, same-direction "
          "users, (3,2) removed on both sides)")


def scenario_partition_sample():
    a = Archive.loads(Path("samples/partition.json").read_text(encoding="utf-8"))
    t, out = check_pipeline(a, ["east", "west"], F(0), "partition")
    assert len(t.arcs) == 4
    assert (F(2), F(1, 2)) not in out["east"].outer
    assert len(out["west"].outer) == 7 and len(out["east"].outer) == 7
    assert out["west"].area == 8 and out["east"].area == F(23, 2)
    t2, out2 = check_pipeline(a, ["east", "west"], F(1, 2), "partition-half")
    assert len(out2["west"].outer) == 4 and len(out2["east"].outer) == 4
    assert len(out2["west"].holes[0]) == 3
    print("ok: bent shared boundary in samples/partition.json (zero and "
          "1/2 tolerance, shared chain updated synchronously)")


def scenario_narrow_channel():
    a = Archive()
    a.put("ch", ring([(0, 0), (2, 0), (2, 1), (3, 1), (3, 0), (5, 0),
                      (5, 3), (3, 3), (3, 2), (2, 2), (2, 3), (0, 3)]))
    _, out = check_pipeline(a, ["ch"], F(1, 2), "channel")
    assert len(out["ch"].outer) == 12 and out["ch"].area == 13
    print("ok: narrow channel of width 1 keeps all 12 vertices at "
          "tolerance 1/2")


def scenario_near_hole():
    gap = F(1, 10 ** 12)
    a = Archive()
    a.put("nh", box(0, 0, 4, 4), (ring([(gap, gap), (1, gap), (1, 1), (gap, 1)]),))
    _, out = check_pipeline(a, ["nh"], F(1, 10 ** 13), "near-hole")
    assert len(out["nh"].outer) == 4 and len(out["nh"].holes[0]) == 4
    assert out["nh"].area == a.read("nh").area
    assert out["nh"].locate((gap / 2, F(1, 2))) == "interior"
    print("ok: near-hole with gap 1/10^12 untouched at tolerance 1/10^13")


def scenario_three_region_junction():
    a = Archive()
    a.put("A", box(0, 0, 1, 1))
    a.put("B", box(1, 0, 2, 1))
    a.put("C", box(0, 1, 1, 2))
    t, out = check_pipeline(a, ["A", "B", "C"], F(1), "junction")
    assert len(t.arcs) == 5
    for rid in ("A", "B", "C"):
        assert len(out[rid].outer) == 3
    print("ok: three-region junction at (1,1) keeps 4-way contact and "
          "cyclic order while corners are removed")


def scenario_point_contact():
    a = Archive()
    a.put("A", box(0, 0, 2, 2))
    a.put("B", ring([(2, 2), (4, 2), (3, 0)]))
    t, out = check_pipeline(a, ["A", "B"], F(1, 2), "contact")
    assert len(t.arcs) == 2
    assert all(not arc.closed and arc.points[0] == arc.points[-1]
               for arc in t.arcs)
    assert (F(2), F(2)) in out["A"].outer and (F(2), F(2)) in out["B"].outer
    print("ok: point contact produces two loop arcs (first==last, "
          "closed=False); contact preserved")


def scenario_closed_arc_reduction():
    a = Archive()
    a.put("box", box(0, 0, 4, 4))
    t, out = check_pipeline(a, ["box"], F(10), "closed")
    assert len(t.arcs) == 1 and t.arcs[0].closed
    assert len(out["box"].outer) == 3 and out["box"].area == 8
    print("ok: independent closed arc reduces to a triangle at large "
          "tolerance, never below three vertices")


def scenario_reversed_input():
    a = Archive()
    a.put("west", ring([(0, 4), (2, 4), (2, 0), (0, 0)]))
    a.put("east", ring([(2, 2), (2, 4), (5, 4), (5, 0), (2, 0)]))
    t, out = check_pipeline(a, ["west", "east"], F(0), "reversed")
    shared = [arc for arc in t.arcs if len(arc.users) == 2][0]
    assert sorted(u[2] for u in shared.users) == [-1, 1]
    assert twice_area(out["west"].outer) < 0
    assert out["west"].area == 8 and out["east"].area == 12
    print("ok: reversed input orientation flips user directions, output "
          "keeps the input direction")


def scenario_illegal_overlap():
    a = Archive()
    a.put("x", box(0, 0, 4, 4))
    a.put("y", box(2, 2, 6, 6))
    for selection in (None, ["x", "y"]):
        try:
            a.build_topology(selection)
        except ValueError:
            pass
        else:
            raise AssertionError("crossing overlap accepted")
    b = Archive()
    b.put("outer", box(0, 0, 4, 4))
    b.put("inside", box(1, 1, 2, 2))
    try:
        b.build_topology()
    except ValueError:
        pass
    else:
        raise AssertionError("contained overlap accepted")
    c = Archive()
    c.put("A", box(0, 0, 6, 6), (box(2, 2, 4, 4),))
    c.put("B", box(2, 2, 4, 4))
    check_pipeline(c, ["A", "B"], F(0), "hole-fill")
    big = Archive()
    big.put("r", ring([(i, 0) for i in range(61)] + [(60, 1), (60, 2)]
                      + [(i, 2) for i in range(58, -1, -1)]))
    try:
        big.build_topology()
    except ValueError:
        pass
    else:
        raise AssertionError("edge limit not enforced")
    print("ok: interior overlap (crossing, containment) and the 120-edge "
          "limit rejected; hole-filling island accepted")


def scenario_new_contact_rejection():
    a = Archive()
    a.put("A", ring([(0, 0), (1, 1), (2, 0), (2, 3), (0, 3)]))
    a.put("B", ring([(F(3, 4), -1), (F(5, 4), -1), (F(5, 4), 0), (F(3, 4), 0)]))
    t, out = check_pipeline(a, ["A", "B"], F(1), "new-contact")
    assert (F(1), F(1)) in out["A"].outer
    assert all(seg_kind(*e1, *e2) == "none"
               for e1 in edges(out["A"].outer) for e2 in edges(out["B"].outer))
    print("ok: tolerance-legal chord that would create a new shared edge "
          "is rejected; regions stay separated")


def scenario_cumulative_error():
    a = Archive()
    a.put("R1", ring([(0, 0), (1, 1), (2, 0), (3, -1), (3, 2), (0, 2)]))
    a.put("R2", ring([(0, 0), (3, -1), (1, -3)]))
    _, out = check_pipeline(a, ["R1", "R2"], F(1), "cumulative")
    assert (F(1), F(1)) not in out["R1"].outer
    assert (F(2), F(0)) in out["R1"].outer
    print("ok: cumulative tolerance keeps (2,0) after (1,1) is removed "
          "(squared distances 2/5 vs 8/5)")


def scenario_immutability_and_replay():
    a = Archive()
    a.put("west", box(0, 0, 2, 4))
    a.put("east", ring([(2, 0), (5, 0), (5, 4), (2, 4), (2, 2)]))
    before = a.dumps()
    t = a.build_topology(["west", "east"])
    snapshot = t.arcs
    first = simplify(t, F(1))
    second = simplify(t, F(1))
    assert first == second
    assert t.arcs == snapshot and a.dumps() == before
    empty = a.build_topology([])
    assert empty.arcs == () and simplify(empty, F(0)) == {}
    print("ok: archive and topology are never mutated; repeated calls "
          "restart from the build-time topology")


def main():
    suite = subprocess.run(
        [sys.executable, "-B", "-m", "unittest", "discover", "-s", "tests", "-v"],
        capture_output=True, text=True)
    sys.stderr.write(suite.stderr)
    if suite.returncode != 0:
        print(suite.stdout)
        raise SystemExit("baseline unittest suite failed")
    print(f"ok: baseline unittest suite "
          f"({suite.stderr.strip().splitlines()[-3].strip()})")
    scenario_straight_shared_edge()
    scenario_hole_island()
    scenario_partition_sample()
    scenario_narrow_channel()
    scenario_near_hole()
    scenario_three_region_junction()
    scenario_point_contact()
    scenario_closed_arc_reduction()
    scenario_reversed_input()
    scenario_illegal_overlap()
    scenario_new_contact_rejection()
    scenario_cumulative_error()
    scenario_immutability_and_replay()
    print("all reference checks passed")


if __name__ == "__main__":
    main()

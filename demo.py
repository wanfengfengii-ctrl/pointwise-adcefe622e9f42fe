"""Read a bounded archive, report exact values, and round-trip its raw rings."""
import argparse
from pathlib import Path
from fractions import Fraction
from polyarchive import Archive, _edges, _segment_kind
from toposimplify import simplify


def ring(coords):
    return tuple((Fraction(x), Fraction(y)) for x, y in coords)


def user_chain(points, j1, j2, direction):
    """Subchain of a ring between two junctions, in ring traversal order."""
    start, stop = (j1, j2) if direction == 1 else (j2, j1)
    n = len(points)
    i = points.index(start)
    chain = []
    while True:
        chain.append(points[i])
        if points[i] == stop:
            return chain
        i = (i + 1) % n


def show_shared_chain(topology, out, arc_index):
    arc = topology.arcs[arc_index]
    j1, j2 = arc.points[0], arc.points[-1]
    chains = []
    for rid, ring_index, direction in arc.users:
        region = out[rid]
        points = region.outer if ring_index == 0 else region.holes[ring_index - 1]
        chains.append((rid, user_chain(points, j1, j2, direction)))
    for rid, chain in chains:
        print(f"  {rid} shared chain {j1} -> {j2}: {chain}")
    first = chains[0][1]
    assert all(chain == first or chain == first[::-1] for _, chain in chains)
    print("  all users of the arc see one identical shared chain")


def topology_demo(archive):
    topology = archive.build_topology()
    print(f"shared topology: {len(topology.arcs)} arcs")
    for index, arc in enumerate(topology.arcs):
        if arc.closed:
            shape = f"closed points={len(arc.points)}"
        else:
            shape = f"{arc.points[0]}..{arc.points[-1]} points={len(arc.points)}"
        print(f"  arc{index}: {shape} users={list(arc.users)}")
    shared = [i for i, arc in enumerate(topology.arcs) if len(arc.users) > 1]
    for tolerance in (Fraction(0), Fraction(1, 2)):
        out = simplify(topology, tolerance)
        print(f"tolerance {tolerance}:")
        for rid in sorted(out):
            region = out[rid]
            print(f"  {rid}: outer={len(region.outer)} vertices "
                  f"holes={len(region.holes)} area={region.area}")
        for index in shared:
            show_shared_chain(topology, out, index)


def hole_island_demo():
    archive = Archive()
    archive.put("A", ring([(0, 0), (6, 0), (6, 6), (0, 6)]),
                (ring([(2, 2), (3, 2), (4, 2), (4, 4), (2, 4)]),))
    archive.put("B", ring([(2, 2), (4, 2), (4, 4), (2, 4)]))
    topology = archive.build_topology()
    print("hole-island topology:")
    for index, arc in enumerate(topology.arcs):
        print(f"  arc{index}: closed={arc.closed} users={list(arc.users)}")
    out = simplify(topology, Fraction(0))
    print(f"  A hole ring after tolerance 0: {out['A'].holes[0]}")
    print(f"  B outer ring after tolerance 0: {out['B'].outer}")
    print(f"  A area={out['A'].area} B area={out['B'].area}")
    probe = (Fraction(3), Fraction(3))
    print(f"  A.locate({probe})={out['A'].locate(probe)} "
          f"B.locate({probe})={out['B'].locate(probe)}")


def contact_guard_demo():
    archive = Archive()
    archive.put("A", ring([(0, 0), (1, 1), (2, 0), (2, 3), (0, 3)]))
    archive.put("B", ring([(Fraction(3, 4), -1), (Fraction(5, 4), -1),
                           (Fraction(5, 4), 0), (Fraction(3, 4), 0)]))
    out = simplify(archive.build_topology(), Fraction(1))
    print("new-contact guard at tolerance 1:")
    print(f"  A outer: {out['A'].outer}")
    print(f"  B outer: {out['B'].outer}")
    separated = all(_segment_kind(*e1, *e2) == "none"
                    for e1 in _edges(out["A"].outer)
                    for e2 in _edges(out["B"].outer))
    print(f"  regions still separated: {separated}")
    assert separated


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", nargs="?", default="samples/partition.json")
    args = parser.parse_args()
    archive = Archive.loads(Path(args.source).read_text(encoding="utf-8"))
    restored = Archive.loads(archive.dumps())
    assert restored.export() == archive.export()
    for region_id in archive.ids:
        region = restored.read(region_id)
        print(f"{region_id}: vertices={len(region.outer)} holes={len(region.holes)} area={region.area}")
    print("west (3/4,3/2):", restored.read("west").locate((Fraction(3, 4), Fraction(3, 2))))
    print("raw rational rings round-trip unchanged")
    topology_demo(restored)
    hole_island_demo()
    contact_guard_demo()


if __name__ == "__main__":
    main()

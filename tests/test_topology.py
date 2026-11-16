import unittest
from fractions import Fraction as F
from pathlib import Path

from polyarchive import Archive, Region
from toposimplify import simplify


def ring(*coords):
    return tuple((F(x), F(y)) for x, y in coords)


def box(x0=0, y0=0, x1=4, y1=4):
    return ring((x0, y0), (x1, y0), (x1, y1), (x0, y1))


def cyclic(ring_value):
    """Canonical cyclic form: rotation-invariant ring comparison."""
    best = min(range(len(ring_value)),
               key=lambda i: ring_value[i:] + ring_value[:i])
    return ring_value[best:] + ring_value[:best]


class BuildTopologyTests(unittest.TestCase):
    def test_straight_shared_edge_arcs(self):
        a = Archive()
        a.put("west", box(0, 0, 2, 4))
        a.put("east", ring((2, 0), (5, 0), (5, 4), (2, 4), (2, 2)))
        t = a.build_topology(["west", "east"])
        self.assertEqual(len(t.arcs), 3)
        shared = [arc for arc in t.arcs if len(arc.users) == 2]
        self.assertEqual(len(shared), 1)
        arc = shared[0]
        self.assertFalse(arc.closed)
        self.assertEqual(set(arc.points),
                         {(F(2), F(0)), (F(2), F(2)), (F(2), F(4))})
        self.assertEqual({arc.points[0], arc.points[-1]},
                         {(F(2), F(0)), (F(2), F(4))})
        self.assertEqual({(u[0], u[1]) for u in arc.users},
                         {("west", 0), ("east", 0)})
        self.assertEqual(sorted(u[2] for u in arc.users), [-1, 1])
        self.assertEqual(len(t.rings["west"]), 1)
        self.assertEqual(len(t.rings["east"]), 1)
        self.assertEqual(len(t.rings["west"][0]), 2)
        self.assertEqual(len(t.rings["east"][0]), 2)

    def test_hole_island_closed_arcs(self):
        b = Archive()
        b.put("A", ring((0, 0), (6, 0), (6, 6), (0, 6)),
              (ring((2, 2), (3, 2), (4, 2), (4, 4), (2, 4)),))
        b.put("B", ring((2, 2), (4, 2), (4, 4), (2, 4)))
        t = b.build_topology()
        self.assertEqual(len(t.arcs), 2)
        self.assertTrue(all(arc.closed for arc in t.arcs))
        shared = [arc for arc in t.arcs if len(arc.users) == 2]
        self.assertEqual(len(shared), 1)
        users = shared[0].users
        self.assertEqual({(u[0], u[1]) for u in users}, {("A", 1), ("B", 0)})
        self.assertEqual(users[0][2], users[1][2])
        self.assertEqual(len(t.rings["A"]), 2)
        self.assertEqual(len(t.rings["B"]), 1)

    def test_selection_rules(self):
        a = Archive()
        a.put("x", box())
        a.put("y", box(5, 5, 6, 6))
        self.assertEqual(len(a.build_topology().arcs), 2)
        self.assertEqual(len(a.build_topology(["x"]).arcs), 1)
        empty = a.build_topology([])
        self.assertEqual(empty.arcs, ())
        self.assertEqual(empty.rings, {})
        for bad in (["x", "x"], "x", [1], [""], [None]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                a.build_topology(bad)
        with self.assertRaises(KeyError):
            a.build_topology(["missing"])

    def test_overlap_rejected(self):
        a = Archive()
        a.put("x", box(0, 0, 4, 4))
        a.put("y", box(2, 2, 6, 6))
        with self.assertRaises(ValueError):
            a.build_topology()
        b = Archive()
        b.put("outer", box(0, 0, 4, 4))
        b.put("inside", box(1, 1, 2, 2))
        with self.assertRaises(ValueError):
            b.build_topology()
        c = Archive()
        c.put("same", box(0, 0, 4, 4))
        c.put("same2", box(0, 0, 4, 4))
        with self.assertRaises(ValueError):
            c.build_topology()
        # regions strictly inside another's hole are legal
        d = Archive()
        d.put("A", box(0, 0, 8, 8), (box(2, 2, 6, 6),))
        d.put("B", box(2, 2, 4, 4))
        d.put("C", box(4, 4, 5, 5))
        self.assertEqual(len(d.build_topology().arcs), 6)
        # the archive itself still stores overlapping regions independently
        self.assertEqual(a.read("x").area, 16)
        self.assertEqual(a.read("y").area, 16)

    def test_edge_limit(self):
        big = Archive()
        big.put("r", ring(*([(i, 0) for i in range(61)]
                            + [(60, 1), (60, 2)]
                            + [(i, 2) for i in range(58, -1, -1)])))
        self.assertEqual(len(big.read("r").outer), 122)
        with self.assertRaises(ValueError):
            big.build_topology()
        ok = Archive()
        ok.put("r", ring(*([(i, 0) for i in range(60)]
                           + [(59, 1), (59, 2)]
                           + [(i, 2) for i in range(57, -1, -1)])))
        self.assertEqual(len(ok.read("r").outer), 120)
        self.assertEqual(len(ok.build_topology().arcs), 1)

    def test_archive_unmodified_and_topology_detached(self):
        a = Archive()
        a.put("west", box(0, 0, 2, 4))
        a.put("east", ring((2, 0), (5, 0), (5, 4), (2, 4), (2, 2)))
        before = a.dumps()
        t = a.build_topology()
        a.put("west", box(0, 0, 1, 1))
        self.assertEqual(len(t.arcs), 3)
        out = simplify(t, F(0))
        self.assertEqual(out["west"].area, 8)
        self.assertNotEqual(a.dumps(), before)
        a2 = Archive()
        a2.put("west", box(0, 0, 2, 4))
        a2.put("east", ring((2, 0), (5, 0), (5, 4), (2, 4), (2, 2)))
        snapshot = a2.dumps()
        t2 = a2.build_topology()
        simplify(t2, F(1))
        self.assertEqual(a2.dumps(), snapshot)


class SimplifyTests(unittest.TestCase):
    def test_spec_example_straight_edge(self):
        a = Archive()
        a.put("west", box(0, 0, 2, 4))
        a.put("east", ring((2, 0), (5, 0), (5, 4), (2, 4), (2, 2)))
        before = a.dumps()
        t = a.build_topology(["west", "east"])
        out = simplify(t, F(0))
        self.assertEqual(a.dumps(), before)
        self.assertEqual(out["west"].area, 8)
        self.assertEqual(out["east"].area, 12)
        self.assertEqual(len(out["west"].outer), 4)
        self.assertEqual(len(out["east"].outer), 4)
        self.assertIsInstance(out["west"], Region)
        self.assertEqual(set(out), {"west", "east"})

    def test_spec_example_hole_island(self):
        b = Archive()
        b.put("A", ring((0, 0), (6, 0), (6, 6), (0, 6)),
              (ring((2, 2), (3, 2), (4, 2), (4, 4), (2, 4)),))
        b.put("B", ring((2, 2), (4, 2), (4, 4), (2, 4)))
        t = b.build_topology()
        out = simplify(t, F(0))
        self.assertEqual(out["A"].area, 32)
        self.assertEqual(out["B"].area, 4)
        self.assertEqual(out["A"].locate((F(3), F(3))), "exterior")
        self.assertEqual(out["B"].locate((F(3), F(3))), "interior")
        self.assertEqual(len(out["A"].holes[0]), 4)
        self.assertEqual(len(out["B"].outer), 4)
        # both sides still reference one identical closed chain
        self.assertEqual(cyclic(out["A"].holes[0]), cyclic(out["B"].outer))

    def test_partition_sample(self):
        source = Path("samples/partition.json").read_text(encoding="utf-8")
        a = Archive.loads(source)
        t = a.build_topology()
        self.assertEqual(len(t.arcs), 4)
        shared = [arc for arc in t.arcs if len(arc.users) == 2]
        self.assertEqual(len(shared), 1)
        self.assertIn((F(2), F(1, 2)), shared[0].points)
        out = simplify(t, F(0))
        self.assertNotIn((F(2), F(1, 2)), out["east"].outer)
        self.assertEqual(len(out["west"].outer), 7)
        self.assertEqual(len(out["east"].outer), 7)
        self.assertEqual(out["west"].area, 8)
        self.assertEqual(out["east"].area, F(23, 2))
        # positive tolerance straightens the shared bend on both sides
        out = simplify(t, F(1, 2))
        self.assertEqual(len(out["west"].outer), 4)
        self.assertEqual(len(out["east"].outer), 4)
        west = out["west"].outer
        i, j = west.index((F(2), F(0))), west.index((F(2), F(4)))
        self.assertEqual((j - i) % len(west), 1)
        east = out["east"].outer
        i, j = east.index((F(2), F(4))), east.index((F(2), F(0)))
        self.assertEqual((j - i) % len(east), 1)

    def test_bad_tolerance(self):
        a = Archive()
        a.put("x", box())
        t = a.build_topology()
        for bad in (0, 1, F(-1), F(-1, 2), 0.5, "0", None, True):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                simplify(t, bad)
        simplify(t, F(0))
        simplify(t, F(1, 2))

    def test_bad_topology(self):
        for bad in (None, {}, [], "x", Region(box()), 0):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                simplify(bad, F(0))

    def test_replay_no_accumulation(self):
        a = Archive()
        a.put("west", box(0, 0, 2, 4))
        a.put("east", ring((2, 0), (5, 0), (5, 4), (2, 4), (2, 2)))
        t = a.build_topology()
        snapshot = t.arcs
        first = simplify(t, F(1))
        self.assertEqual(t.arcs, snapshot)
        second = simplify(t, F(1))
        self.assertEqual(first, second)
        again = simplify(t, F(0))
        self.assertEqual(len(again["west"].outer), 4)

    def test_new_contact_rejected(self):
        a = Archive()
        a.put("A", ring((0, 0), (1, 1), (2, 0), (2, 3), (0, 3)))
        a.put("B", ring((F(3, 4), -1), (F(5, 4), -1), (F(5, 4), 0), (F(3, 4), 0)))
        out = simplify(a.build_topology(), F(1))
        self.assertIn((F(1), F(1)), out["A"].outer)
        self.assertEqual(len(out["A"].outer), 5)
        from polyarchive import _edges, _segment_kind
        self.assertTrue(all(
            _segment_kind(*e1, *e2) == "none"
            for e1 in _edges(out["A"].outer) for e2 in _edges(out["B"].outer)))

    def test_cumulative_error(self):
        a = Archive()
        a.put("R1", ring((0, 0), (1, 1), (2, 0), (3, -1), (3, 2), (0, 2)))
        a.put("R2", ring((0, 0), (3, -1), (1, -3)))
        out = simplify(a.build_topology(), F(1))
        self.assertNotIn((F(1), F(1)), out["R1"].outer)
        self.assertIn((F(2), F(0)), out["R1"].outer)
        self.assertEqual(len(out["R1"].outer), 5)

    def test_empty_selection(self):
        a = Archive()
        a.put("x", box())
        self.assertEqual(simplify(a.build_topology([]), F(0)), {})

    def test_reversed_input(self):
        a = Archive()
        a.put("west", ring((0, 4), (2, 4), (2, 0), (0, 0)))
        a.put("east", ring((2, 2), (2, 4), (5, 4), (5, 0), (2, 0)))
        t = a.build_topology()
        shared = [arc for arc in t.arcs if len(arc.users) == 2][0]
        self.assertEqual(sorted(u[2] for u in shared.users), [-1, 1])
        out = simplify(t, F(0))
        self.assertEqual(len(out["west"].outer), 4)
        self.assertEqual(len(out["east"].outer), 4)
        from polyarchive import _twice_area
        self.assertLess(_twice_area(out["west"].outer), 0)
        self.assertEqual(out["west"].area, 8)
        self.assertEqual(out["east"].area, 12)

    def test_loop_arc_point_contact(self):
        a = Archive()
        a.put("A", ring((0, 0), (2, 0), (2, 2), (0, 2)))
        a.put("B", ring((2, 2), (4, 2), (3, 0)))
        t = a.build_topology()
        self.assertEqual(len(t.arcs), 2)
        for arc in t.arcs:
            self.assertFalse(arc.closed)
            self.assertEqual(arc.points[0], arc.points[-1])
        out = simplify(t, F(1, 2))
        self.assertIn((F(2), F(2)), out["A"].outer)
        self.assertIn((F(2), F(2)), out["B"].outer)

    def test_closed_arc_reduction(self):
        a = Archive()
        a.put("box", box())
        t = a.build_topology()
        self.assertEqual(len(t.arcs), 1)
        self.assertTrue(t.arcs[0].closed)
        out = simplify(t, F(10))
        self.assertEqual(len(out["box"].outer), 3)
        self.assertEqual(out["box"].area, 8)
        unchanged = simplify(t, F(0))
        self.assertEqual(cyclic(unchanged["box"].outer), cyclic(box()))


if __name__ == "__main__":
    unittest.main()

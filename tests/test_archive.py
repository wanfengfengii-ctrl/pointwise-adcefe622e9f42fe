import dataclasses
from fractions import Fraction as F
import json
from pathlib import Path
import random
import unittest

from polyarchive import Archive, Region


def ring(*coords):
    return tuple((F(x), F(y)) for x, y in coords)


def box(x0=0, y0=0, x1=4, y1=4):
    return ring((x0, y0), (x1, y0), (x1, y1), (x0, y1))


class RegionTests(unittest.TestCase):
    def test_area_orientation_collinear_vertices(self):
        outer = ring((0, 0), (2, 0), (4, 0), (4, 4), (0, 4))
        holes = (box(1, 1, 2, 2), box(3, 1, F(7, 2), 2))
        a = Region(outer, holes)
        b = Region(outer[::-1], tuple(h[::-1] for h in holes[::-1]))
        self.assertEqual(a.area, F(29, 2))
        self.assertEqual(a.area, b.area)
        self.assertEqual(a.outer, outer)
        self.assertEqual(a.holes, holes)

    def test_rectangle_oracle(self):
        rng = random.Random(2903)
        for _ in range(20):
            x, y = rng.randrange(-5, 5), rng.randrange(-5, 5)
            w, h = rng.randrange(3, 8), rng.randrange(3, 8)
            poly = Region(box(x, y, x + w, y + h), (box(x+1, y+1, x+2, y+2),))
            self.assertEqual(poly.area, w * h - 1)
            for i in range(-1, 2*w + 2):
                for j in range(-1, 2*h + 2):
                    px, py = F(i, 2), F(j, 2)
                    boundary = ((px in (0, w) and 0 <= py <= h)
                                or (py in (0, h) and 0 <= px <= w)
                                or (px in (1, 2) and 1 <= py <= 2)
                                or (py in (1, 2) and 1 <= px <= 2))
                    inside = 0 < px < w and 0 < py < h and not (1 < px < 2 and 1 < py < 2)
                    expected = "boundary" if boundary else "interior" if inside else "exterior"
                    self.assertEqual(poly.locate((x+px, y+py)), expected)

    def test_concavity_narrow_channel_and_near_hole(self):
        channel = ring((0,0),(2,0),(2,1),(3,1),(3,0),(5,0),
                       (5,3),(3,3),(3,2),(2,2),(2,3),(0,3))
        poly = Region(channel)
        self.assertEqual(poly.area, 13)
        self.assertEqual(poly.locate((F(5,2), F(3,2))), "interior")
        self.assertEqual(poly.locate((F(5,2), F(1,2))), "exterior")
        gap = F(1, 10**30)
        near = Region(box(), (box(gap, gap, 1, 1),))
        self.assertEqual(near.locate((gap/2, F(1,2))), "interior")
        self.assertEqual(near.area, 16 - (1-gap)**2)

    def test_invalid_rings(self):
        invalid = [
            (), ring((0,0),(1,1)),
            ring((0,0),(1,0),(2,0)),
            ring((0,0),(2,0),(2,2),(0,0)),
            ring((0,0),(2,0),(1,0),(1,2),(0,2)),
            ring((0,0),(4,0),(0,3),(3,3)),
            ring((0,0),(4,0),(4,4),(2,0),(0,4)),
        ]
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                Region(value)
        for value in (None, "ring", [(0, F(0)), (F(1), F(0)), (F(0), F(1))],
                      [(F(0), F(0), F(0))] * 3):
            with self.subTest(value=value), self.assertRaises(ValueError):
                Region(value)

    def test_invalid_holes(self):
        invalid = [
            (box(5,5,6,6),), (box(0,1,1,2),),
            (box(-1,1,1,2),), (box(1,1,3,3), box(2,2, F(7,2), F(7,2))),
            (box(1,1,3,3), box(F(3,2), F(3,2),2,2)),
            (box(1,1,2,2), box(2,2,3,3)),
            (box(1,1,2,2), box(1,1,2,2)),
        ]
        for holes in invalid:
            with self.subTest(holes=holes), self.assertRaises(ValueError):
                Region(box(), holes)

    def test_input_aliases_and_immutable_values(self):
        outer = [list(p) for p in box()]
        hole = [list(p) for p in box(1,1,2,2)]
        holes = [hole]
        poly = Region(outer, holes)
        outer[0][0] = F(99)
        hole[0][0] = F(99)
        holes.clear()
        self.assertEqual(poly.outer, box())
        self.assertEqual(poly.holes, (box(1,1,2,2),))
        with self.assertRaises(dataclasses.FrozenInstanceError):
            poly.outer = box(0,0,2,2)
        with self.assertRaises(ValueError):
            poly.locate((0, 0))


class ArchiveTests(unittest.TestCase):
    def test_replacement_and_errors_are_atomic(self):
        archive = Archive()
        self.assertIsNone(archive.put("a", box()))
        saved = archive.read("a")
        before = archive.dumps()
        for region_id, outer, holes in [
            ("a", box(), (box(0,1,1,2),)), ("", box(), ()),
            (True, box(), ()), ("new", ring((0,0),(1,0),(2,0)), ()),
        ]:
            with self.subTest(region_id=region_id), self.assertRaises(ValueError):
                archive.put(region_id, outer, holes)
            self.assertEqual(archive.dumps(), before)
        archive.put("a", box(0,0,2,2))
        self.assertEqual(saved.area, 16)
        self.assertEqual(archive.read("a").area, 4)
        with self.assertRaises(KeyError):
            archive.read("missing")

    def test_independent_regions_need_no_common_segmentation(self):
        archive = Archive()
        archive.put("west", box(0,0,2,4))
        east = ring((2,0),(5,0),(5,4),(2,4),(2,2))
        archive.put("east", east)
        archive.put("point", box(5,4,6,5))
        archive.put("overlapping", box(1,1,3,3))
        self.assertEqual(archive.read("east").outer, east)
        self.assertEqual(archive.read("west").area, 8)
        self.assertEqual(archive.read("point").locate((F(5),F(4))), "boundary")

    def test_wire_sample_and_detached_export(self):
        source = Path("samples/partition.json").read_text(encoding="utf-8")
        archive = Archive.loads(source)
        self.assertEqual(archive.ids, ("east", "west"))
        self.assertEqual(archive.read("west").area, 8)
        self.assertEqual(archive.read("east").area, F(23,2))
        original = json.loads(source)
        by_id = {r["id"]: r for r in original["regions"]}
        self.assertEqual(archive.export()["regions"], [by_id[k] for k in archive.ids])
        self.assertEqual(Archive.loads(archive.dumps()).dumps(), archive.dumps())
        exported = archive.export()
        exported["regions"][0]["outer"][0][0] = "99"
        self.assertEqual(archive.read("east").outer[0], (F(2), F(0)))
        self.assertEqual(archive.export([]), {"version":1, "regions":[]})
        self.assertEqual(len(archive.export(["west"])["regions"]), 1)
        for selection in (["west","west"], "west", [1]):
            with self.assertRaises(ValueError):
                archive.export(selection)
        with self.assertRaises(KeyError):
            archive.export(["unknown"])

    def test_malformed_wire(self):
        valid = Archive()
        valid.put("r", box())
        raw = valid.export()
        texts = ['{', 'null', '{"version":1,"version":1,"regions":[]}',
                 '{"version":1,"regions":[],"extra":0}',
                 '{"version":true,"regions":[]}',
                 '{"version":2,"regions":[]}', '{"version":1,"regions":null}',
                 '{"version":1,"regions":[],"x":NaN}']
        for coord in (0, True, "0/1", "2/2", "1/0", "NaN", "1.0", "-0"):
            data = json.loads(json.dumps(raw))
            data["regions"][0]["outer"][0][0] = coord
            texts.append(json.dumps(data))
        for alter in ("duplicate", "extra", "missing", "holes", "point"):
            data = json.loads(json.dumps(raw))
            if alter == "duplicate":
                data["regions"] *= 2
            elif alter == "extra":
                data["regions"][0]["extra"] = 1
            elif alter == "missing":
                del data["regions"][0]["id"]
            elif alter == "holes":
                data["regions"][0]["holes"] = None
            else:
                data["regions"][0]["outer"][0] = ["0"]
            texts.append(json.dumps(data))
        for text in texts:
            with self.subTest(text=text), self.assertRaises(ValueError):
                Archive.loads(text)
        with self.assertRaises(ValueError):
            Archive.loads(b"{}")
        self.assertEqual(Archive.loads('{"version":1,"regions":[]}').ids, ())


if __name__ == "__main__":
    unittest.main()

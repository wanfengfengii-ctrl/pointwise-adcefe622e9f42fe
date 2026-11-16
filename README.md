# Rational region archive

Python 3.11 standard library. Run from the project directory:

    python3.11 -B -m unittest discover -s tests -v
    python3.11 -B demo.py samples/partition.json

The archive stores independent planar regions with exact rational coordinates.
No projection or coordinate rounding occurs. The JSON fixture contains two
adjacent regions with a bent common boundary, unequal segmentation and a hole.

## Python API

    from fractions import Fraction as F
    from polyarchive import Archive, Region

    a = Archive()
    a.put("district", [(F(0),F(0)), (F(4),F(0)),
                       (F(4),F(4)), (F(0),F(4))])
    region = a.read("district")
    assert region.area == F(16)
    assert region.locate((F(1),F(1))) == "interior"
    text = a.dumps()
    assert Archive.loads(text).read("district") == region

Region(outer, holes=()) accepts list/tuple rings and copies them into immutable
tuples. A point is exactly two Fraction coordinates; ints, floats and booleans
are rejected. A ring is implicitly closed, has at least three vertices and
nonzero area. Do not repeat the first point at the end. Rings must be simple:
nonadjacent edges cannot meet; adjacent edges meet only at their common endpoint.
Consecutive collinear vertices are allowed unless they backtrack. Orientation,
starting vertex, collinear points and hole order are preserved.

Each hole is a simple ring strictly inside the outer ring. Hole closures must be
pairwise disjoint; touching, overlapping and nested holes are rejected. There
are no islands inside a Region's holes; store an island as a separate region.
area subtracts the unsigned hole areas independently of orientation.
locate(point) returns "interior", "boundary" or "exterior"; hole edges are
boundary and hole interiors are exterior. Calculations use Fraction throughout.

Archive() starts empty. put(id, outer, holes=()) inserts or replaces one
region and returns None. IDs are nonempty strings. Validation completes before
replacement; failure raises ValueError and preserves all entries. read(id)
returns an immutable Region; a retained value stays unchanged after replacement.
Unknown valid IDs raise KeyError. ids returns IDs as a sorted tuple.

Different regions are stored independently: overlapping interiors, point
contacts, shared edges and different boundary segmentation are all accepted.
Insertion does not certify that a group of regions forms a partition.

export(region_ids=None) returns fresh JSON-compatible data for all IDs, or a
list/tuple selection. Duplicate selections raise ValueError; unknown IDs raise
KeyError. Records are sorted by ID regardless of selection order. An empty
selection exports an empty archive. dumps(region_ids=None) serializes that
data deterministically with a final newline. These operations never change the
archive. Archive.loads(text) constructs a new archive from version 1 text;
malformed input raises ValueError and cannot affect an existing archive.

## Version 1 interchange

The root has exactly "version" (integer 1) and "regions" (array). Each record has
exactly "id", "outer" and "holes". Rings are arrays of two-element coordinate
arrays; each coordinate is the canonical string produced by str(Fraction),
such as "0", "-3", or "5/2". Decimal spellings, non-reduced fractions, duplicate
object keys, duplicate IDs and unknown fields are rejected. Record order may
vary on input; export sorts IDs. Vertex and hole order do not vary on export.
The format is illustrated by samples/partition.json.

The Python API handles geometry and text; callers own any filesystem access.
The demo reads only the supplied project sample by default and prints an exact
area summary after a serialization round trip. Tests include exact
rectangle-and-hole point-location oracles, narrow channels, near-touching
rational geometry, invalid shapes, immutable reads and format rejection.
Validation deliberately checks all edge pairs.

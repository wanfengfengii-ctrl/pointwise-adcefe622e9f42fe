"""Read a bounded archive, report exact values, and round-trip its raw rings."""
import argparse
from pathlib import Path
from fractions import Fraction
from polyarchive import Archive


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


if __name__ == "__main__":
    main()

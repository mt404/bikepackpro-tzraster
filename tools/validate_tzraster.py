#!/usr/bin/env python3
"""
Validate a generated BikepackPro time zone raster against its source polygons.

The raster quantises the world onto a grid, so it is wrong in exactly one way: for
a cell that a zone boundary passes through, the whole cell takes one zone and the
sliver on the other side of the boundary is misreported. This script measures that,
and asserts the only property that actually matters -- *every* disagreement with the
source polygons lies within half a cell diagonal of a real boundary. A mismatch
further out than that is not quantisation, it is a bug in the generator.

Three sampling strategies, because they answer different questions:

  uniform   -- points spread evenly over the globe. Mostly ocean, so it mainly
               proves the raster is not broken somewhere nobody looks.
  land      -- points rejected until they miss the Etc/* maritime zones. Closer to
               where riders are, and the denominator most comparable to published
               figures from other libraries.
  route     -- dense samples along polylines that deliberately cross zone
               boundaries. This is the number to care about: people, roads and bike
               routes cluster along rivers and political lines, and political lines
               are frequently the same lines as zone boundaries, so an area-weighted
               error rate systematically understates what a rider experiences.

The reader here is written from the format documentation rather than shared with the
generator, so a disagreement between them catches a spec error rather than hiding it.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import struct
import sys
import time

import numpy as np
import shapely
from shapely.geometry import Point, shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from generate_tzraster import LEGACY_ZONE_NAMES  # noqa: E402

EARTH_R = 6_371_000.0


# --------------------------------------------------------------------------- #
# An independent reader for the .bin format
# --------------------------------------------------------------------------- #

class TzRaster:
    def __init__(self, path: str):
        with open(path, "rb") as handle:
            self.blob = handle.read()
        (magic, version, _, self.width, self.height, zone_count, _,
         palette_off, row_index_off, self.runs_off) = struct.unpack_from(">4sHHIIHHIII", self.blob, 0)
        if magic != b"BPTZ":
            raise ValueError(f"bad magic {magic!r}")
        if version != 1:
            raise ValueError(f"unsupported version {version}")

        self.palette = []
        pos = palette_off
        for _ in range(zone_count):
            (length,) = struct.unpack_from(">H", self.blob, pos)
            pos += 2
            self.palette.append(self.blob[pos:pos + length].decode("utf-8"))
            pos += length
        if pos != row_index_off:
            raise ValueError(f"palette ends at {pos}, row index starts at {row_index_off}")

        self.row_index = np.frombuffer(
            self.blob, dtype=">u4", count=self.height, offset=row_index_off
        ).astype(np.int64)

    def zone_at(self, lat: float, lon: float) -> str | None:
        col = int((lon + 180.0) / 360.0 * self.width)
        row = int((90.0 - lat) / 180.0 * self.height)
        col = min(max(col, 0), self.width - 1)
        row = min(max(row, 0), self.height - 1)

        pos = self.runs_off + int(self.row_index[row])
        end = (self.runs_off + int(self.row_index[row + 1])
               if row + 1 < self.height else len(self.blob))
        x = 0
        blob = self.blob
        while pos < end:
            length = 0
            shift = 0
            while True:
                byte = blob[pos]
                pos += 1
                length |= (byte & 0x7F) << shift
                shift += 7
                if not byte & 0x80:
                    break
            index = (blob[pos] << 8) | blob[pos + 1]
            pos += 2
            x += length
            if col < x:
                return self.palette[index] if index != 0xFFFF else None
        return None


# --------------------------------------------------------------------------- #
# Truth
# --------------------------------------------------------------------------- #

class Truth:
    """Ground truth straight from the source polygons, with the generator's rules."""

    def __init__(self, geojson_path: str):
        with open(geojson_path, encoding="utf-8") as handle:
            features = json.load(handle)["features"]
        self.tzids = [f["properties"]["tzid"] for f in features]
        # What the generator would have emitted for this zone.
        self.emitted = [LEGACY_ZONE_NAMES.get(t, t) for t in self.tzids]
        self.geoms = [shape(f["geometry"]) for f in features]
        self.areas = np.array([g.area for g in self.geoms])
        # Several of these polygons carry millions of vertices (Russia, Canada, the
        # Etc/* bands). Without prepared geometries every containment test walks all
        # of them and the validation run takes hours instead of minutes.
        shapely.prepare(self.geoms)
        self.tree = STRtree(self.geoms)
        # Extracting .boundary from a MultiPolygon with millions of vertices costs
        # real time, and the displacement pass asks for the same handful of zones
        # over and over, so build each at most once.
        self._boundaries: dict[int, object] = {}
        self._by_emitted: dict[str, list[int]] = {}
        for i, name in enumerate(self.emitted):
            self._by_emitted.setdefault(name, []).append(i)

    def _boundary(self, index: int):
        cached = self._boundaries.get(index)
        if cached is None:
            cached = self.geoms[index].boundary
            self._boundaries[index] = cached
        return cached

    def zones_for(self, lons: np.ndarray, lats: np.ndarray) -> list[str | None]:
        """Vectorised containment. Smallest-area polygon wins, matching the burn order."""
        points = shapely.points(lons, lats)
        # STRtree evaluates the predicate as input.PREDICATE(tree_geometry), not the
        # other way round, so the point has to be the subject: "covered_by", not
        # "covers". Getting this backwards returns zero pairs and silently reports a
        # perfect score, which is exactly how a validator lies to you.
        # "covered_by" rather than "within" so a point exactly on a boundary counts.
        pairs = self.tree.query(points, predicate="covered_by")
        best: dict[int, int] = {}
        for point_i, geom_i in zip(pairs[0].tolist(), pairs[1].tolist()):
            current = best.get(point_i)
            if current is None or self.areas[geom_i] < self.areas[current]:
                best[point_i] = geom_i
        return [self.emitted[best[i]] if i in best else None for i in range(len(points))]

    def distance_to_boundary_m(self, lat: float, lon: float, zone_emitted: str) -> float:
        """Metres from the point to the nearest edge of any polygon emitting this zone."""
        point = Point(lon, lat)
        best = math.inf
        for i in self._by_emitted.get(zone_emitted, ()):
            if self.geoms[i].distance(point) > 5.0:  # degrees; cheap reject
                continue
            near = nearest_points(self._boundary(i), point)[0]
            best = min(best, _metres(lat, lon, near.y, near.x))
        return best


def _metres(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Local equirectangular metres; exact enough at the sub-kilometre scale here."""
    mean_lat = math.radians((lat1 + lat2) / 2.0)
    dx = math.radians(lon2 - lon1) * math.cos(mean_lat) * EARTH_R
    dy = math.radians(lat2 - lat1) * EARTH_R
    return math.hypot(dx, dy)


# --------------------------------------------------------------------------- #
# Sample sets
# --------------------------------------------------------------------------- #

def uniform_samples(n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = random.Random(seed)
    lats = np.array([rng.uniform(-90, 90) for _ in range(n)])
    lons = np.array([rng.uniform(-180, 180) for _ in range(n)])
    return lats, lons


def land_samples(n: int, seed: int, truth: Truth) -> tuple[np.ndarray, np.ndarray]:
    """Rejection-sample until the point is not in an Etc/* maritime zone."""
    rng = random.Random(seed)
    lats: list[float] = []
    lons: list[float] = []
    rounds = 0
    while len(lats) < n:
        rounds += 1
        if rounds > 200:
            raise RuntimeError(
                f"land sampling stalled at {len(lats)}/{n}; the containment query is "
                f"probably returning nothing (check the STRtree predicate direction)")
        batch = 4 * (n - len(lats)) + 64
        cand_lat = np.array([rng.uniform(-60, 72) for _ in range(batch)])
        cand_lon = np.array([rng.uniform(-180, 180) for _ in range(batch)])
        zones = truth.zones_for(cand_lon, cand_lat)
        for lat, lon, zone in zip(cand_lat, cand_lon, zones):
            if zone and not zone.startswith("Etc/"):
                lats.append(float(lat))
                lons.append(float(lon))
                if len(lats) == n:
                    break
    return np.array(lats), np.array(lons)


# Representative bikepacking-shaped routes that cross a zone boundary. The repo has
# no GPX fixtures of any kind (verified: `find -name '*.gpx'` returns nothing), so
# these are constructed. Each is a straight leg between two real places chosen so the
# leg crosses at least one boundary that changes the clock.
ROUTES = [
    ("Idaho panhandle, Pacific to Mountain",
     (46.3800, -116.9800), (46.3800, -114.9000)),
    ("South Dakota, Mountain to Central",
     (44.0800, -102.5000), (44.0800, -99.5000)),
    ("Arizona to Utah, Phoenix to Denver rules",
     (36.5000, -111.5000), (38.5000, -111.5000)),
    ("Spain to Portugal, Madrid to Lisbon",
     (41.0000,   -6.4000), (41.0000,  -8.2000)),
    ("Pyrenees, France to Spain (control, no change)",
     (43.2000,    0.5000), (42.3000,   0.5000)),
    ("South Australia to New South Wales (half-hour step)",
     (-34.0000, 140.0000), (-34.0000, 142.5000)),
    ("Oregon coast, inland to territorial water",
     (45.0000, -123.9000), (45.0000, -124.6000)),
    ("Kazakhstan, Qostanay remap region",
     (52.0000,   62.0000), (50.0000,   66.0000)),
]


def route_samples(spacing_m: float) -> list[tuple[str, np.ndarray, np.ndarray]]:
    out = []
    for name, (lat1, lon1), (lat2, lon2) in ROUTES:
        length = _metres(lat1, lon1, lat2, lon2)
        steps = max(2, int(length / spacing_m))
        fractions = np.linspace(0.0, 1.0, steps)
        lats = lat1 + (lat2 - lat1) * fractions
        lons = lon1 + (lon2 - lon1) * fractions
        out.append((name, lats, lons))
    return out


# --------------------------------------------------------------------------- #
# Comparison
# --------------------------------------------------------------------------- #

def compare(label: str, raster: TzRaster, truth: Truth,
            lats: np.ndarray, lons: np.ndarray, cell_diag_m: float,
            report_worst: int = 5) -> dict:
    expected = truth.zones_for(lons, lats)
    mismatches = []
    counted = 0
    for lat, lon, want in zip(lats.tolist(), lons.tolist(), expected):
        if want is None:
            continue
        counted += 1
        got = raster.zone_at(lat, lon)
        if got != want:
            mismatches.append((lat, lon, want, got))

    result = {
        "label": label,
        "samples": counted,
        "mismatches": len(mismatches),
        "mismatch_rate": len(mismatches) / counted if counted else 0.0,
    }

    worst = 0.0
    beyond_bound = []
    detail = []
    for lat, lon, want, got in mismatches:
        distance = truth.distance_to_boundary_m(lat, lon, want)
        detail.append((distance, lat, lon, want, got))
        worst = max(worst, distance)
        if distance > cell_diag_m:
            beyond_bound.append((distance, lat, lon, want, got))

    result["worst_displacement_m"] = worst
    result["beyond_half_diagonal"] = len(beyond_bound)

    print(f"  {label:52s} n={counted:>8,}  mismatches={len(mismatches):>6,} "
          f"({100 * result['mismatch_rate']:6.3f}%)  worst={worst:7.1f} m")
    if detail:
        detail.sort(reverse=True)
        for distance, lat, lon, want, got in detail[:report_worst]:
            flag = "  <-- BEYOND BOUND" if distance > cell_diag_m else ""
            print(f"      {distance:8.1f} m from a {want} boundary at "
                  f"({lat:.5f}, {lon:.5f}); raster said {got}{flag}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raster", required=True)
    parser.add_argument("--geojson", required=True)
    parser.add_argument("--uniform", type=int, default=200_000)
    parser.add_argument("--land", type=int, default=100_000)
    parser.add_argument("--route-spacing-m", type=float, default=50.0)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--json-out")
    args = parser.parse_args()

    started = time.time()
    raster = TzRaster(args.raster)
    resolution = 360.0 / raster.width
    # Latitude degrees are ~111.32 km everywhere; longitude degrees shrink with
    # cos(lat), so the equator is the worst case for the diagonal.
    cell_lat_m = resolution * 111_320.0
    cell_diag_m = math.hypot(cell_lat_m, cell_lat_m) / 2.0
    print(f"raster {args.raster}: {raster.width} x {raster.height} "
          f"({resolution} deg), {len(raster.palette)} zones")
    print(f"half-diagonal bound at the equator: {cell_diag_m:.1f} m\n")

    truth = Truth(args.geojson)
    print(f"truth: {len(truth.geoms)} source polygons loaded "
          f"({time.time() - started:.0f}s)\n")

    results = []

    print("uniform global sampling")
    lats, lons = uniform_samples(args.uniform, args.seed)
    results.append(compare("uniform (all surfaces)", raster, truth, lats, lons, cell_diag_m))

    print("\nland sampling (Etc/* maritime zones rejected)")
    lats, lons = land_samples(args.land, args.seed + 1, truth)
    results.append(compare("land only", raster, truth, lats, lons, cell_diag_m))

    print(f"\nroute sampling at {args.route_spacing_m:.0f} m spacing")
    for name, lats, lons in route_samples(args.route_spacing_m):
        results.append(compare(name, raster, truth, lats, lons, cell_diag_m))

    total_beyond = sum(r["beyond_half_diagonal"] for r in results)
    print("\n" + "=" * 78)
    print(f"samples          {sum(r['samples'] for r in results):,}")
    print(f"mismatches       {sum(r['mismatches'] for r in results):,}")
    print(f"worst any set    {max(r['worst_displacement_m'] for r in results):.1f} m "
          f"(bound {cell_diag_m:.1f} m)")
    print(f"beyond the bound {total_beyond}")
    print("=" * 78)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as handle:
            json.dump({"cell_diagonal_bound_m": cell_diag_m, "results": results},
                      handle, indent=2)

    if total_beyond:
        sys.exit(f"FAILED: {total_beyond} mismatches lie further than half a cell "
                 f"diagonal from any boundary -- that is a generator bug, not quantisation")
    print("PASS: every mismatch lies within half a cell diagonal of a real boundary")


if __name__ == "__main__":
    main()

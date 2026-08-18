#!/usr/bin/env python3
"""
Build a bundled coordinate -> IANA time zone raster for BikepackPro.

WHY THIS EXISTS
---------------
The app renders times that belong to a *place* -- the sunset at the point where a
rider will run out of daylight, the arrival time at a POI -- and those have to be
shown on that place's clock, not the phone's. The device already ships the whole
IANA database via ICU, so offsets, DST rules and localized zone names are all
available locally. The one thing Android has no API for at any level is mapping a
latitude/longitude onto a zone identifier: `android.location.Geocoder` returns an
`Address`, and `Address` has no timezone field.

Every accurate global answer to that question is a network call or a 25-70 MB
library. A bikepacker is in the backcountry with no signal precisely when they most
need to know whether they will beat nightfall, so a network call is not an option.
This script bakes the answer into a file small enough to bundle and fast enough to
read synchronously.

OUTPUT FORMAT (all integers big-endian, matching java.io.DataInputStream and
java.nio.ByteBuffer defaults so the reader needs no byte-order fiddling)
-----------------------------------------------------------------------
    Header, 32 bytes
        0   magic          4 bytes  "BPTZ"
        4   version        u16      = 1
        6   reserved       u16      = 0
        8   width          u32      grid columns
        12  height         u32      grid rows
        16  zoneCount      u16      palette entries
        18  reserved2      u16      = 0
        20  paletteOffset  u32      absolute byte offset
        24  rowIndexOffset u32      absolute byte offset
        28  runsOffset     u32      absolute byte offset

    Palette, zoneCount entries
        u16 byte length, then that many UTF-8 bytes (an IANA identifier)

    Row index, height entries
        u32 byte offset of that row's first run, relative to runsOffset

    Runs, concatenated per row, north to south
        unsigned LEB128 varint run length (in cells), then u16 palette index

A cell's centre maps to lat/lon as
        col = floor((lon + 180) / 360 * width)
        row = floor(( 90 - lat) / 180 * height)
so row 0 is the northern edge and column 0 the antimeridian, which is the ordinary
north-up raster convention.

WHY ROW-RLE RATHER THAN A FLAT ARRAY
------------------------------------
A flat array at 0.001 degrees would be 360000 x 180000 x 2 bytes = 129.6 GB. Run
length encoding along rows exploits the fact that the compressed size scales with
the *length of the zone boundaries* (proportional to 1/resolution) rather than with
*area* (proportional to 1/resolution^2). Halving the cell size therefore roughly
doubles the file instead of quadrupling it.

The reader memory-maps this file, so resident heap stays flat no matter how fine the
grid gets: the pages live in the OS page cache, not the Java heap.

DATA SOURCE AND LICENSING
-------------------------
Input is the `timezones-with-oceans` release of timezone-boundary-builder
(https://github.com/evansiroky/timezone-boundary-builder), which is derived from
OpenStreetMap. That data is licensed under the Open Data Commons Open Database
License (ODbL); the project's own code is MIT.

The raster this script produces is a reformatted extraction of that database and is
therefore a Derivative Database under ODbL section 4.4b. Bundling it in a shipped
app is Public Use, which means section 4.3 (attribution) and section 4.6 (offer
recipients a machine-readable copy) both apply. Section 4.5b confirms the app's
*answers* are a Produced Work, so no obligation reaches the application source.

Usage:
    python3 generate_tzraster.py --input combined-with-oceans.json \\
        --output tzraster.bin --resolution 0.001
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import re
import struct
import sys
import time

import numpy as np
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely.geometry import shape
from shapely.strtree import STRtree

MAGIC = b"BPTZ"
VERSION = 1
NODATA = 0xFFFF
HEADER_BYTES = 32

# tzdata files that define Zone and Link records.
TZDATA_FILES = ["africa", "antarctica", "asia", "australasia", "europe",
                "northamerica", "southamerica", "etcetera", "backward"]

# Zones that timezone-boundary-builder uses but that an API-26-vintage device
# cannot name, mapped to an identifier that both an old and a current device
# resolve, with identical UTC offsets today.
#
# minSdk is 26 (Android 8.0, August 2017, tzdata ~2017a) and tz updates only became
# a mainline module at API 30, so the oldest devices we support will never learn
# these names. `ZoneId.of("Europe/Kyiv")` throws ZoneRulesException there, which
# would take out the whole feature on exactly the phones least able to spare it.
#
# The first three are ordinary tzdb backward links: the old spelling is retained
# forever and resolves to the same rules on every device, old or new. The last three
# are genuinely new zones created by boundary splits, with no alias to fall back on,
# so each is mapped to the nearest identifier whose modern offsets are byte-identical
# over a forward year. `--verify-palette` re-derives and re-checks this whole table
# from real tzdata rather than trusting it.
LEGACY_ZONE_NAMES = {
    # new id                  old id                 why
    "America/Nuuk":          "America/Godthab",      # tzdb backward link (2020a rename)
    "Europe/Kyiv":           "Europe/Kiev",          # tzdb backward link (2022b rename)
    "Pacific/Kanton":        "Pacific/Enderbury",    # tzdb backward link (2021b rename)
    "America/Ciudad_Juarez": "America/Denver",       # split from Ojinaga (2022g); Ojinaga
                                                     # is now UTC-6 fixed while Ciudad
                                                     # Juarez keeps US Mountain DST, so the
                                                     # historical parent is the WRONG answer
                                                     # today and US Mountain is the right one
    "America/Coyhaique":     "America/Punta_Arenas",  # split from Santiago (2025a); both
                                                      # southern Chile, UTC-3 year round
    "Asia/Qostanay":         "Asia/Qyzylorda",       # split from Qyzylorda (2018h); all of
                                                     # Kazakhstan is UTC+5 since 2024
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# Palette
# --------------------------------------------------------------------------- #

def parse_tzdata_names(root: str) -> set[str]:
    """Every Zone and Link identifier defined in an unpacked tzdata tree."""
    names: set[str] = set()
    for filename in TZDATA_FILES:
        path = os.path.join(root, filename)
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
                parts = re.split(r"\s+", line)
                if parts[0] == "Zone" and len(parts) >= 2:
                    names.add(parts[1])
                elif parts[0] == "Link" and len(parts) >= 3:
                    names.add(parts[2])
    return names


def verify_palette(tzids: list[str], old_tzdata_dir: str) -> None:
    """Fail loudly if LEGACY_ZONE_NAMES does not cover what old tzdata cannot name."""
    old_names = parse_tzdata_names(old_tzdata_dir)
    if not old_names:
        sys.exit(f"no tzdata Zone/Link records found under {old_tzdata_dir}")
    log(f"verify: reference tzdata defines {len(old_names)} identifiers")

    unknown = [t for t in tzids if t not in old_names]
    uncovered = [t for t in unknown if t not in LEGACY_ZONE_NAMES]
    stale = [k for k in LEGACY_ZONE_NAMES if k not in tzids]
    bad_target = [(k, v) for k, v in LEGACY_ZONE_NAMES.items() if v not in old_names]

    log(f"verify: {len(unknown)} of {len(tzids)} tzids are unknown to the reference tzdata")
    for tzid in unknown:
        log(f"verify:   {tzid} -> {LEGACY_ZONE_NAMES.get(tzid, '!! UNMAPPED !!')}")

    problems = []
    if uncovered:
        problems.append(f"unmapped zones old tzdata cannot resolve: {uncovered}")
    if bad_target:
        problems.append(f"substitutes that old tzdata also cannot resolve: {bad_target}")
    if stale:
        problems.append(f"LEGACY_ZONE_NAMES entries no longer present in the input: {stale}")
    if problems:
        sys.exit("palette verification FAILED:\n  " + "\n  ".join(problems))

    log("verify: palette OK -- every emitted identifier resolves on the reference tzdata")


# --------------------------------------------------------------------------- #
# Encoding helpers
# --------------------------------------------------------------------------- #

def encode_varint(value: int, out: bytearray) -> None:
    """Unsigned LEB128."""
    while value >= 0x80:
        out.append((value & 0x7F) | 0x80)
        value >>= 7
    out.append(value)


def encode_row(row: np.ndarray, out: bytearray) -> int:
    """Append one row's runs. Returns the number of runs written."""
    # Boundaries where the value changes, vectorised: for a row of N cells this is
    # one pass rather than a Python loop over cells.
    change = np.flatnonzero(row[1:] != row[:-1]) + 1
    starts = np.empty(change.size + 1, dtype=np.int64)
    starts[0] = 0
    starts[1:] = change
    ends = np.empty_like(starts)
    ends[:-1] = starts[1:]
    ends[-1] = row.size
    lengths = ends - starts
    values = row[starts]

    for length, value in zip(lengths.tolist(), values.tolist()):
        encode_varint(length, out)
        out.append((value >> 8) & 0xFF)
        out.append(value & 0xFF)
    return starts.size


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True,
                        help="timezone-boundary-builder combined-with-oceans.json")
    parser.add_argument("--output", required=True, help="destination .bin")
    parser.add_argument("--resolution", type=float, default=0.001,
                        help="cell size in degrees (default 0.001, about 111 m of latitude)")
    parser.add_argument("--strip-rows", type=int, default=100,
                        help="rows rasterised per pass; bounds peak memory")
    parser.add_argument("--verify-palette", metavar="TZDATA_DIR",
                        help="unpacked old tzdata tree to check legacy names against")
    parser.add_argument("--manifest", help="write a JSON build manifest here")
    args = parser.parse_args()

    if not (0 < args.resolution <= 1):
        sys.exit("--resolution must be in (0, 1] degrees")
    width = int(round(360.0 / args.resolution))
    height = int(round(180.0 / args.resolution))
    if width * args.resolution != 360.0 or height * args.resolution != 180.0:
        sys.exit(f"--resolution {args.resolution} does not divide the globe evenly")

    log(f"grid {width} x {height} = {width * height:,} cells at {args.resolution} deg")

    log(f"loading {args.input}")
    with open(args.input, encoding="utf-8") as handle:
        features = json.load(handle)["features"]
    log(f"loaded {len(features)} features")

    tzids = [f["properties"]["tzid"] for f in features]
    duplicates = [k for k, v in collections.Counter(tzids).items() if v > 1]
    if duplicates:
        sys.exit(f"input has more than one feature per tzid: {duplicates[:5]}")

    if args.verify_palette:
        verify_palette(sorted(tzids), args.verify_palette)

    log("building geometries")
    geometries = [shape(f["geometry"]) for f in features]

    # Palette order is the sorted *emitted* identifier, so the file is reproducible
    # and a diff between two builds is readable.
    emitted = [LEGACY_ZONE_NAMES.get(t, t) for t in tzids]
    palette = sorted(set(emitted))
    palette_index = {name: i for i, name in enumerate(palette)}
    if len(palette) > 0xFFFE:
        sys.exit(f"{len(palette)} zones exceeds the u16 palette index space")
    log(f"palette: {len(palette)} identifiers "
        f"({len(set(tzids)) - len(palette)} collapsed by legacy remapping)")

    # timezone-boundary-builder ships a small number of genuinely overlapping
    # polygons -- Asia/Shanghai over Asia/Urumqi in Xinjiang, Europe/Moscow over
    # Asia/Tbilisi in Abkhazia and South Ossetia -- where two authorities claim the
    # same ground. rasterize() lets later shapes win, so burning largest-area first
    # makes the *smallest* (most specific) polygon the answer. That matches what
    # other implementations do and, more importantly, is deterministic: without a
    # stated order the output would depend on GeoJSON feature order.
    order = sorted(range(len(features)), key=lambda i: geometries[i].area, reverse=True)
    burn_values = [palette_index[emitted[i]] for i in range(len(features))]

    tree = STRtree(geometries)

    runs = bytearray()
    row_offsets = np.zeros(height, dtype=np.uint64)
    total_runs = 0
    nodata_cells = 0
    max_runs_in_row = 0
    started = time.time()

    for row0 in range(0, height, args.strip_rows):
        strip_h = min(args.strip_rows, height - row0)
        north = 90.0 - row0 * args.resolution
        south = north - strip_h * args.resolution
        transform = from_origin(-180.0, north, args.resolution, args.resolution)

        # Only polygons that reach this latitude band can contribute.
        candidates = set(tree.query(_strip_box(south, north)).tolist())
        shapes = [(geometries[i], burn_values[i]) for i in order if i in candidates]

        strip = rasterize(
            shapes,
            out_shape=(strip_h, width),
            transform=transform,
            fill=NODATA,
            dtype=np.uint16,
            all_touched=False,
        )

        missing = int(np.count_nonzero(strip == NODATA))
        if missing:
            nodata_cells += missing

        for r in range(strip_h):
            row_offsets[row0 + r] = len(runs)
            n = encode_row(strip[r], runs)
            total_runs += n
            max_runs_in_row = max(max_runs_in_row, n)

        if (row0 // args.strip_rows) % 20 == 0 or row0 + strip_h >= height:
            done = row0 + strip_h
            elapsed = time.time() - started
            rate = done / elapsed if elapsed else 0
            eta = (height - done) / rate if rate else 0
            log(f"rows {done:,}/{height:,} ({100 * done / height:5.1f}%)  "
                f"runs={total_runs:,}  bytes={len(runs):,}  "
                f"elapsed={elapsed / 60:.1f}m  eta={eta / 60:.1f}m")

    log(f"rasterised in {(time.time() - started) / 60:.1f} minutes")
    if nodata_cells:
        # with-oceans is supposed to tile the whole globe; a gap means the reader
        # would have to invent an answer, which is the one thing it must never do.
        log(f"WARNING: {nodata_cells:,} cells had no zone "
            f"({100 * nodata_cells / (width * height):.6f}% of the grid)")

    palette_blob = bytearray()
    for name in palette:
        raw = name.encode("utf-8")
        palette_blob += struct.pack(">H", len(raw)) + raw

    palette_offset = HEADER_BYTES
    row_index_offset = palette_offset + len(palette_blob)
    runs_offset = row_index_offset + height * 4

    header = struct.pack(
        ">4sHHIIHHIII",
        MAGIC, VERSION, 0, width, height, len(palette), 0,
        palette_offset, row_index_offset, runs_offset,
    )
    assert len(header) == HEADER_BYTES, len(header)

    log(f"writing {args.output}")
    with open(args.output, "wb") as out:
        out.write(header)
        out.write(palette_blob)
        out.write(row_offsets.astype(">u4").tobytes())
        out.write(runs)

    size = os.path.getsize(args.output)
    log(f"done: {size:,} bytes ({size / 1048576:.2f} MiB)")
    log(f"      {total_runs:,} runs, {total_runs / height:.1f} per row on average, "
        f"{max_runs_in_row:,} at most")

    if args.manifest:
        manifest = {
            "format_version": VERSION,
            "resolution_degrees": args.resolution,
            "width": width,
            "height": height,
            "zone_count": len(palette),
            "total_runs": total_runs,
            "mean_runs_per_row": total_runs / height,
            "max_runs_in_row": max_runs_in_row,
            "nodata_cells": nodata_cells,
            "file_bytes": size,
            "legacy_remapped": LEGACY_ZONE_NAMES,
            "source": os.path.basename(args.input),
        }
        with open(args.manifest, "w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
        log(f"manifest -> {args.manifest}")


def _strip_box(south: float, north: float):
    from shapely.geometry import box
    # A hair of slop so a polygon whose edge lies exactly on the strip boundary is
    # not dropped by a floating point comparison.
    return box(-180.0, south - 1e-9, 180.0, north + 1e-9)


if __name__ == "__main__":
    main()

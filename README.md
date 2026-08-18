# BikepackPro time zone raster

A bundled dataset that maps a latitude/longitude onto an IANA time zone identifier
entirely offline, plus the tooling that builds and validates it.

This repository exists to satisfy the Open Database License. The raster shipped
inside the BikepackPro app is a reformatted extraction of
[timezone-boundary-builder](https://github.com/evansiroky/timezone-boundary-builder),
which is derived from OpenStreetMap, and is therefore a **Derivative Database** under
ODbL §4.4b. Distributing it in an app is Public Use, so §4.6 requires that a
machine-readable copy be offered to recipients. This is that copy.

## Why an offline raster

BikepackPro shows times that belong to a *place* — the sunset at the point where a
rider will run out of daylight, the arrival time at a point of interest — and those
have to be rendered on that place's clock rather than the phone's. Android already
ships the whole IANA database through ICU, so offsets, DST rules and localized zone
names are all available on-device. The one thing it has no API for, at any level, is
mapping a coordinate onto a zone identifier: `android.location.Geocoder` returns an
`Address`, and `Address` carries no timezone field.

Every accurate alternative is either a network call or a library bundling 25–70 MB of
polygon data with a heap cost to match. A bikepacker is in the backcountry with no
signal precisely when they most need to know whether they will beat nightfall, so a
network call is not an option. This file is the answer baked flat.

## Provenance

| | |
|---|---|
| Source project | timezone-boundary-builder |
| Release | **2026c**, published 2026-07-11 |
| Asset | `timezones-with-oceans.geojson.zip` |
| Asset size | 55,449,715 bytes |
| Asset sha256 | `70bc2f5e9b48f49461368fafe30daa840c3db9c7bc730faf33614cad339f9b1d` |
| Zones in source | 444 |
| | |
| Output | `tzraster.bin` |
| Output size | 33,687,830 bytes |
| Output sha256 | `4f25382fa740a466f75444fc07a76dfce3150fa8a5d3455ba7e8e341d6e2a8b9` |
| Resolution | 0.001° (about 111 m of latitude) |
| Grid | 360,000 × 180,000 |
| Zones in palette | 441 |

The generator is deterministic: the same input and flags produce a byte-identical
output, which has been verified rather than assumed.

## What was changed from the source

Three things, and nothing else. Together these are the complete set of differences
between the source database and this derivative.

1. **Rasterised** onto a 0.001° grid and run-length encoded along rows. Row RLE is
   what makes the file small: compressed size tracks the *length* of the zone
   boundaries, proportional to 1/resolution, rather than their *area*, proportional
   to 1/resolution². A flat array at this resolution would be 129.6 GB.

2. **Overlaps resolved deterministically.** The source contains genuinely overlapping
   polygons covering roughly 0.28% of the globe's area — `Asia/Shanghai` over
   `Asia/Urumqi` in Xinjiang, `Europe/Moscow` over `Asia/Tbilisi` in Abkhazia and
   South Ossetia — where two authorities claim the same ground. Polygons are burned
   largest-area first, so the smallest and most specific one wins. Without a stated
   rule the output would depend on GeoJSON feature order.

3. **Six identifiers remapped to older spellings**, so the file is readable on old
   Android devices. See below.

## The legacy identifier remapping

BikepackPro supports Android 8.0 (API 26), which shipped with tzdata around 2017a,
and time zone data only became an updatable mainline module at API 30. Those devices
will never learn identifiers added since, and `ZoneId.of("Europe/Kyiv")` throws
`ZoneRulesException` there — which would take the feature out entirely on exactly the
phones least able to spare it.

Six zones are therefore emitted under an older name. Every substitution was verified
by comparing 8,760 hourly UTC offsets across a forward year, not by reasoning about
which names ought to be safe.

| Source identifier | Emitted | Basis |
|---|---|---|
| `America/Nuuk` | `America/Godthab` | tzdb backward link (2020a rename) |
| `Europe/Kyiv` | `Europe/Kiev` | tzdb backward link (2022b rename) |
| `Pacific/Kanton` | `Pacific/Enderbury` | tzdb backward link (2021b rename) |
| `America/Ciudad_Juarez` | `America/Denver` | offsets byte-identical over a forward year |
| `America/Coyhaique` | `America/Punta_Arenas` | offsets byte-identical over a forward year |
| `Asia/Qostanay` | `Asia/Qyzylorda` | offsets byte-identical over a forward year |

The first three are ordinary tzdb backward links, which every device resolves in both
directions. The last three are genuine boundary splits with no alias to fall back on,
so each is mapped to the nearest identifier whose *modern* offsets match exactly —
the substitute has to be right on a current device as well as nameable on an old one.

`America/Ciudad_Juarez` is the instructive case. Its historical parent
`America/Ojinaga` is the intuitive substitute and is **wrong**: Ojinaga moved to fixed
UTC−6 while Ciudad Juárez kept US Mountain DST, so the historically correct answer
would misreport the time on any current device. US Mountain is the right answer.

Three of the six collapse onto identifiers already in the palette, which is why 444
source zones become 441 palette entries.

## Accuracy

Quantising onto a grid is wrong in exactly one way: for a cell that a boundary passes
through, the whole cell takes one zone and the sliver on the other side is
misreported. The bound is half a cell diagonal, **78.7 m** at the equator and less
toward the poles.

Measured over 320,626 samples:

| Sample set | n | Mismatches | Rate | Worst displacement |
|---|---|---|---|---|
| Uniform global | 200,000 | 8 | 0.004% | 22.8 m |
| Land only | 60,000 | 4 | 0.007% | 30.4 m |
| Along boundary-crossing routes | 60,626 | 26 | 0.043% | 50.9 m |

**Every** disagreement with the source polygons lies within the half-diagonal bound;
none is beyond it. That is the property the validator asserts, and it exits non-zero
if it ever fails — a disagreement further out would be a generator bug rather than
quantisation.

Route sampling is weighted most heavily because people, roads and bike routes cluster
along rivers and political lines, and political lines are frequently the same lines as
zone boundaries. An area-weighted error rate systematically understates what a rider
actually meets. The full report is in `validation-report.txt`.

## Known gaps

- The source has exactly **one cell with no zone** in 64.8 billion: a 111 m square at
  75.7485° N, 67.1655° W, in the Nares Strait between Greenland and Ellesmere Island.
  It is encoded as palette index `0xFFFF`. A reader must surface that as "unknown"
  rather than inventing an answer — showing a confident wrong time is the outcome this
  whole dataset exists to prevent.
- Maritime boundaries follow territorial waters, not Exclusive Economic Zones, and a
  few open-water borders are the source project's best guess. Both are inherited from
  upstream.

## File format

All integers big-endian, matching `java.io.DataInputStream` and `java.nio.ByteBuffer`
defaults so a reader needs no byte-order handling.

```
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
    u32 offset of that row's first run, relative to runsOffset

Runs, concatenated per row, north to south
    unsigned LEB128 varint run length in cells, then u16 palette index
```

Lookup:

```
col = floor((lon + 180) / 360 * width)
row = floor(( 90 - lat) / 180 * height)
```

The file is designed to be memory-mapped. Doing so keeps resident heap flat regardless
of grid resolution, because the pages live in the operating system's page cache rather
than the application heap.

## Regenerating

```sh
python3 -m venv venv
./venv/bin/pip install -r tools/requirements.txt

curl -L -o tz.zip \
  https://github.com/evansiroky/timezone-boundary-builder/releases/download/2026c/timezones-with-oceans.geojson.zip
unzip tz.zip          # -> combined-with-oceans.json

curl -L -o tzdata2017a.tar.gz https://data.iana.org/time-zones/releases/tzdata2017a.tar.gz
mkdir -p tz2017a && tar xzf tzdata2017a.tar.gz -C tz2017a

./venv/bin/python tools/generate_tzraster.py \
  --input combined-with-oceans.json \
  --output tzraster.bin \
  --resolution 0.001 \
  --strip-rows 500 \
  --verify-palette tz2017a \
  --manifest tzraster.manifest.json
```

Roughly 7 minutes and about 1.5 GB of RAM. Then:

```sh
./venv/bin/python tools/validate_tzraster.py \
  --raster tzraster.bin \
  --geojson combined-with-oceans.json \
  --uniform 200000 --land 60000 --route-spacing-m 25
```

`rasterio`'s wheels vendor their own GDAL, so no system GDAL installation is required.

## Licence

The source polygons come from timezone-boundary-builder, derived from OpenStreetMap,
and are licensed under the
[Open Data Commons Open Database License (ODbL) v1.0](https://opendatacommons.org/licenses/odbl/1-0/).

`tzraster.bin` is a Derivative Database of that data and is published here under the
same licence. The full text is in `LICENSE`.

The scripts under `tools/` are published under ODbL alongside the data they produce,
so that §4.6(b) — "the method of making the alterations to the Database (such as an
algorithm)" — is satisfied as well as §4.6(a).

Contains information from timezone-boundary-builder, made available under the Open
Database License (ODbL), derived from OpenStreetMap data © OpenStreetMap contributors.

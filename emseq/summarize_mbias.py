
#!/usr/bin/env python3
import re
import os
import sys

def parse(path):
    with open(path) as f:
        s = f.read()

    # Read numeric y-axis tick labels and their pixel coordinates.
    ticks = re.findall(
        r'<text x="55" y="([\d.]+)"[^>]*>([\d.]+)</text>', s
    )
    ticks = [(float(y), float(v)) for y, v in ticks]

    if len(ticks) < 2:
        raise ValueError(f"Cannot read y-axis ticks in {path}")

    # Derive the linear pixel-to-value conversion from the tick labels.
    (y1, v1), (y2, v2) = ticks[0], ticks[-1]
    slope = (v2 - v1) / (y2 - y1)

    # Use the SVG plot's left and bottom coordinates.
    series = []
    for d, attr in re.findall(r'<path d="([^"]*)"([^>]*)>', s):
        if 'stroke=' not in attr:
            continue

        xy = re.findall(r'([ML])\s+([\d.]+)\s+([\d.]+)', d)
        if not xy:
            continue

        # Convert pixel coordinates to read position and plotted percentage.
        xvals = [float(x) for _, x, _ in xy]
        xmin, xmax = min(xvals), max(xvals)
        positions = {}

        for _, x, y in xy:
            pos = round((float(x) - 80) / (500 / (len(xy) + 4)))
            value = (float(y) - y1) * slope + v1
            positions[pos] = value * 100

        series.append(positions)

    sug = re.findall(r'--O[TB] [\d,]+', s)
    return series, (sug[0] if sug else "NA")


print("\t".join([
    "library", "strand", "read", "plateau",
    "first6", "last6", "suggested"
]))

for path in sys.argv[1:]:
    base = os.path.basename(path)
    lib, strand = base[:-4].rsplit("_", 1)

    series, sug = parse(path)

    for name, pts in zip(["R1", "R2"], series):
        if not pts:
            continue

        n = max(pts)
        mid = [pts[i] for i in range(20, n - 19) if i in pts]
        first = " ".join(
            "%.0f" % pts[i] for i in range(1, 7) if i in pts
        )
        last = " ".join(
            "%.0f" % pts[i] for i in range(n - 5, n + 1) if i in pts
        )
        plateau = sum(mid) / len(mid) if mid else float("nan")

        print("\t".join([
            lib, strand, name, "%.2f" % plateau,
            first, last, sug
        ]))

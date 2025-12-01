from Bio import SeqIO
import csv
import math
import sys

INPUT_FASTA = sys.argv[1]
OUTPUT_CSV = sys.argv[2]

PROBE_LEN = 120
UNIVERSAL_PRIMER = "CGAAGAGCCCTATAGTGAGTCGTATTAGAA"

def tile_sequence(seq, window=120):
    """Tile a sequence with 1x tiling density without overlapping windows."""
    tiles = []
    total = len(seq)
    num_tiles = math.ceil(total / window)

    for i in range(num_tiles):
        start = i * window
        end = min(start + window, total)
        tiles.append(str(seq[start:end]))
    
    return tiles


records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

output_rows = []

for record in records:
    exon_name = record.id
    seq = record.seq

    tiles = tile_sequence(seq, window=PROBE_LEN)

    for i, tile in enumerate(tiles, start=1):

        probe_name = f"{exon_name}_probe{i}"

        probe_120 = tile.upper()

        twist_oligo = probe_120 + UNIVERSAL_PRIMER

        output_rows.append([
            probe_name,
            probe_120,
            UNIVERSAL_PRIMER,
            twist_oligo
        ])


# Write output CSV

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.writer(f, delimiter=",")
    writer.writerow([
        "Probe name",
        "Probe sequence (120 nt)",
        "Nt.BspQI sequence (30 nt)",
        "Twist oligo sequence (150 nt)"
    ])
    writer.writerow(output_rows)

print(f"Done! Wrote {len(output_rows)} probes to {OUTPUT_CSV}")


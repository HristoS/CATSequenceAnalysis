"""Origin rotation at scale, second key: the skew-minimum origin of every rotated chromosome (origin_dnaa_check.cs's CSV) against
DoriC 12.0's oriC and, where one exists, an experimentally verified oriC.

Both keys come from Zenodo record 10888581 (van Meel 2024, CC-BY-4.0; DoriC 12.0 from Dong, Luo and Gao, NAR 2023):
  DoriC_complete_circular.csv  md5 6fa4e6ec45927a6b613c9e6052bb5e37  (refseq = assembly, nc = chromosome accession; has oriC DNA,
                               only coordinates are read; delete the file afterwards)
  experimental_dataset.csv     md5 9692dcff6a0abace9745e65b7422f747
Records match a DoriC chromosome by accession without NZ_ or version (an NZ_ RefSeq copy has the GenBank coordinates). An oriC whose
end is below its start runs past position 0. DoriC comes from Ori-Finder, which itself weighs a Z-curve GC skew, DnaA boxes and dnaA,
so agreement with it is not fully independent; the experimental set is. Distances in windows of 32 kbp, the level rotated chromosomes
are scored at; the stop line (plan, 2026-11-12) is 80% within +-5.

Usage: python -I origin_doric_check.py S32_OriginDnaA.csv DoriC_complete_circular.csv experimental_dataset.csv
"""
import collections
import csv
import statistics
import sys

window = 32_000
csv.field_size_limit(10**9)


def key(accession):
    return accession.split('.')[0].removeprefix('NZ_')


def middle(start, end, length):
    return (start + end) // 2 if end >= start else ((start + end + length) // 2) % length


def distance(a, b, length):
    d = abs(a - b) % length
    return min(d, length - d)


def share(values, windows):
    return sum(v <= windows * window for v in values) / len(values) if values else float('nan')


def median_kbp(values):
    return statistics.median(values) / 1000 if values else float('nan')


def line(label, rows):
    skew = [r['skew_to_doric'] for r in rows]
    dnaa = [r['dnaa_to_doric'] for r in rows if r['dnaa_to_doric'] is not None]
    print(f"  {label:<24} n={len(rows):5}  skew-DoriC {share(skew, 1):6.1%} {share(skew, 5):6.1%} median {median_kbp(skew):7.1f} kbp"
          f" | dnaA-DoriC {share(dnaa, 1):6.1%} {share(dnaa, 5):6.1%} median {median_kbp(dnaa):7.1f} kbp")


def main(origin_path, doric_path, experimental_path):
    records = list(csv.DictReader(open(origin_path, newline=''), skipinitialspace=True))
    doric = collections.defaultdict(list)
    with open(doric_path, newline='') as f:
        for r in csv.DictReader(f):
            doric[key(r['nc'])].append((int(r['oric_start']), int(r['oric_end'])))
    verified = collections.defaultdict(list)
    with open(experimental_path, newline='') as f:
        for r in csv.DictReader(f):
            verified[key(r['accession'])].append((int(r['begin']), int(r['end']), int(r['seq_len']), r['organism']))

    matched = []
    for r in records:
        accession, length, origin = key(r['Accession']), int(r['LengthBp']), int(r['OriginBp'])
        dnaa = int(r['DnaAStartBp']) if r['DnaAStartBp'] else None
        oris = doric.get(accession)
        if not oris or any(max(o) > length for o in oris):
            continue
        middles = [middle(s, e, length) for s, e in oris]
        matched.append(dict(phylum=r['Phylum'], genus=r['Genus'], skew=float(r['SkewBpPerKbp']),
                            skew_to_doric=min(distance(origin, m, length) for m in middles),
                            dnaa_to_doric=None if dnaa is None else min(distance(dnaa, m, length) for m in middles)))

    print(f"rotated records {len(records):,}; matched to a DoriC chromosome {len(matched):,}")
    print("within +-1 / +-5 windows of 32 kbp; median distance")
    line("all", matched)
    by_phylum = collections.defaultdict(list)
    for m in matched:
        by_phylum[m['phylum']].append(m)
    for phylum, rows in sorted(by_phylum.items(), key=lambda kv: -len(kv[1]))[:12]:
        line(phylum, rows)
    for low, high in [(3, 5), (5, 10), (10, 20), (20, float('inf'))]:
        line(f"skew [{low}, {high})", [m for m in matched if low <= m['skew'] < high])

    print("experimentally verified oriCs (same chromosome accession and length):")
    for r in records:
        accession, length = key(r['Accession']), int(r['LengthBp'])
        oris = verified.get(accession)
        if not oris or any(o[2] != length for o in oris):
            continue
        middles = [middle(s, e, length) for s, e, _, _ in oris]
        skew = min(distance(int(r['OriginBp']), m, length) for m in middles)
        dnaa = min(distance(int(r['DnaAStartBp']), m, length) for m in middles) if r['DnaAStartBp'] else None
        print(f"  {oris[0][3][:45]:<45} {r['Phylum']:<17} skew {skew / 1000:8.1f} kbp   dnaA "
              f"{'-' if dnaa is None else f'{dnaa / 1000:.1f}'} kbp   skew strength {float(r['SkewBpPerKbp']):.1f}")


if __name__ == '__main__':
    main(*sys.argv[1:4])

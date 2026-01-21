import argparse
import re
from collections import defaultdict
from statistics import median

## This script corrects misassemblies around telomeres based on coverage
##      Gets rid of low coverage tails outside of terminal telomeres
##      and splits contigs at internal telomeres on the side with low coverage

TEL_RE = re.compile(r'(TTAGGG){6,}|(CCCTAA){6,}', re.IGNORECASE)

def read_fasta(fa):
    seqs = {}
    name = None
    buf = []
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)
    return seqs

def read_coverage(tsv):
    cov = defaultdict(dict)
    with open(tsv) as f:
        for line in f:
            c, pos, d = line.strip().split()
            cov[c][int(pos)] = int(d)
    return cov

def mean_cov(cov_dict, start, end):
    vals = [cov_dict.get(i, 0) for i in range(start, end)]
    return sum(vals) / max(1, (end - start))

## find telomeric repeats in a sequence and merge nearby hits
def find_telomeres(seq, merge_distance=50):
    hits = []
    for m in TEL_RE.finditer(seq):
        hits.append((m.start(), m.end()))
    if not hits:
        return []
    ## sort hits by start coordinate
    hits.sort(key=lambda x: x[0])
    ## merge nearby hits
    merged_hits = []
    current_start, current_end = hits[0]
    for start, end in hits[1:]:
        if start - current_end <= merge_distance:
            current_end = max(current_end, end)
        else:
            merged_hits.append((current_start, current_end))
            current_start, current_end = start, end
    merged_hits.append((current_start, current_end))
    return merged_hits


ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--coverage", required=True)
ap.add_argument("--window", type=int, default=5000)
ap.add_argument("--break_bed", default="breaks.bed")
ap.add_argument("--keep_bed", default="keep_regions.bed")
args = ap.parse_args()

seqs = read_fasta(args.fasta)
cov = read_coverage(args.coverage)

with open(args.keep_bed, "w") as keep_out:
    for contig, seq in seqs.items():
        clen = len(seq)
        print(f"processing: {contig}")

        ## skip contigs with no coverage info
        if contig not in cov or not cov[contig]:
            keep_out.write(f"{contig}\t0\t{clen}\n")
            print(f"   No coverage info for {contig}. Skipping")
            continue

        cov_vals = list(cov[contig].values())
        med_cov = median(cov_vals)
        low = 0.1 * med_cov
        high = 0.5 * med_cov

        telos = find_telomeres(seq)
        if telos:
            print("   Telomeres:")
            for telo in telos:
                print(f"      {contig}\t{telo[0]}\t{telo[1]}")

        ## keep entire contig by default
        keep_start, keep_end = 0, clen
        trimmed = False
        cut = None

        ## process individual telomeric repeats
        for start, end in telos:
            left_cov = mean_cov(cov[contig], max(1, start - args.window), start)
            right_cov = mean_cov(cov[contig], end, min(clen, end + args.window))

            ## internal telomere -> break
            if end <= clen * 0.95 and start >= clen * 0.05:
                ## keep telomere on the more supported side
                if left_cov >= right_cov:
                    cut = end      ## telomere stays on left contig
                    print(f"   Internal telomere break: {contig}, {cut}, telomere stays on left contig")
                else:
                    cut = start    ## telomere stays on right contig
                    print(f"   Internal telomere break: {contig}, {cut}, telomere stays on right contig")
            ## terminal telomere with low coverage tail -> trim
            else:
                if left_cov > high and right_cov < low and end > clen * 0.95 and end != clen and clen - end > 6:
                    keep_start, keep_end = keep_start, end
                    trimmed = True
                    num_trimmed = clen - end
                    print(f"   Trim outside terminal telomere: {contig}, {keep_end}, right side (trimmed {num_trimmed} bases)")

                if right_cov > high and left_cov < low and start < clen * 0.05 and start != 0 and start > 6:
                    keep_start, keep_end = start, keep_end
                    trimmed = True
                    num_trimmed = start
                    print(f"   Trim outside terminal telomere: {contig}, {keep_start}, left side (trimmed {num_trimmed} bases)")

        if cut:
            keep_out.write(f"{contig}\t{keep_start}\t{cut}\n")
            keep_out.write(f"{contig}\t{cut+1}\t{keep_end}\n")
        else:
            keep_out.write(f"{contig}\t{keep_start}\t{keep_end}\n")


import sys


def parse_paf(path, minimap2=False):
    result = {}

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue

            query = parts[0]
            q_len = int(parts[1])
            q_start = int(parts[2])
            q_end = int(parts[3])
            strand = parts[4]
            target = parts[5]
            t_start = int(parts[7])
            t_end = int(parts[8])
            block_len = int(parts[10])
            mapq = int(parts[11])

            tags = parts[12:]
            is_primary = (not minimap2) or ("tp:A:P" in tags)

            if minimap2 and not is_primary:
                continue

            rec = {
                "query": query,
                "q_len": q_len,
                "q_start": q_start,
                "q_end": q_end,
                "strand": strand,
                "target": target,
                "t_start": t_start,
                "t_end": t_end,
                "block_len": block_len,
                "mapq": mapq,
            }

            if query not in result:
                result[query] = rec
            else:
                old = result[query]
                if (rec["mapq"], rec["block_len"]) > (old["mapq"], old["block_len"]):
                    result[query] = rec

    return result


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


if len(sys.argv) != 3:
    print("Usage: python compare_paf.py <our.paf> <minimap2.paf>")
    sys.exit(1)

our = parse_paf(sys.argv[1], minimap2=False)
mm2 = parse_paf(sys.argv[2], minimap2=True)

common = sorted(set(our) & set(mm2))

strand_match = 0
target_match = 0
overlap_50 = 0
overlap_80 = 0

both_strand_and_overlap80 = 0
strand_only = 0
overlap80_only = 0
neither = 0


for q in common:
    a = our[q]
    b = mm2[q]

    ov = overlap(a["t_start"], a["t_end"], b["t_start"], b["t_end"])
    shorter = min(a["t_end"] - a["t_start"], b["t_end"] - b["t_start"])
    frac = ov / shorter if shorter > 0 else 0.0

    same_strand = a["strand"] == b["strand"]
    good_overlap80 = frac >= 0.8

    if same_strand and good_overlap80:
        both_strand_and_overlap80 += 1
    elif same_strand and not good_overlap80:
        strand_only += 1
    elif not same_strand and good_overlap80:
        overlap80_only += 1
    else:
        neither += 1

    if a["strand"] == b["strand"]:
        strand_match += 1

    if a["target"] == b["target"]:
        target_match += 1

    if frac >= 0.5:
        overlap_50 += 1

    if frac >= 0.8:
        overlap_80 += 1

print(f"Our mapped: {len(our)}")
print(f"Minimap2 mapped primary: {len(mm2)}")
print(f"Mapped by both: {len(common)}")

if common:
    print(f"Target match: {target_match}/{len(common)} = {target_match / len(common) * 100:.2f}%")
    print(f"Strand match: {strand_match}/{len(common)} = {strand_match / len(common) * 100:.2f}%")
    print(f"Target overlap >= 50%: {overlap_50}/{len(common)} = {overlap_50 / len(common) * 100:.2f}%")
    print(f"Target overlap >= 80%: {overlap_80}/{len(common)} = {overlap_80 / len(common) * 100:.2f}%")
    print()
    print("Strand/overlap contingency table using overlap >= 80%:")
    print(f"Both strand match and overlap >= 80%: {both_strand_and_overlap80}")
    print(f"Strand match only: {strand_only}")
    print(f"Overlap >= 80% only: {overlap80_only}")
    print(f"Neither: {neither}")
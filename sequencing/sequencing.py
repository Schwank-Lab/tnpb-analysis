"""FASTQ anchoring and QC. Extracted from the Cell_Screen sequencing package; only the functions reached from this repository's notebooks are kept."""

from pathlib import Path
from Bio import SeqIO
import pandas as pd
import gzip
import subprocess
import numpy as np


def open_fastq(fastq_path, mode='r'):
    if isinstance(fastq_path, Path): 
        fastq_path = str(fastq_path)
    if fastq_path.endswith('.gz'):
        return gzip.open(fastq_path, mode+'t')
    else:
        return open(fastq_path, mode)


def run_fastqc(fastq: Path, fastqc_dir: Path):
    """Run FastQC on the specified file."""
    cmd = f"fastqc -o {str(fastqc_dir)} {str(fastq)}"
    print(cmd)
    # Replace ShellOperation with subprocess.run
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)
    print(f"FastQC analysis completed for {fastq}")


def decompose_read_into_parts(read_id, r, found_anchor_pos, anchor_idx, part_lens, part_names):
    """Splits read into parts accoring to the layout."""

    layout_anchor_start = sum(part_lens[:anchor_idx])
    part_seqs = []
    read_len = len(r)
    for i in range(len(part_lens)):
        layout_part_start = sum(part_lens[:i])
        offset = layout_part_start - layout_anchor_start
        pos_in_read = found_anchor_pos + offset
        part_len = part_lens[i]
        if pos_in_read < 0:
            print(f'Warning: read {read_id} part {part_names[i]} start is outside of the read, make sure the part_offset is correct')
            return None
        elif pos_in_read + part_len > read_len:
            print(f'Warning: read {read_id} part {part_names[i]} end is outside of the read, make sure the part_len is correct')
            return None
        part_seq = r[pos_in_read: pos_in_read + part_len]
        part_seqs.append(part_seq)
    return part_seqs


def anchor_reads(sample, layout, anchor_seq, fastq: Path, out_csv: Path, reverse_complement=False, part_anchor='anchor'):
    """
    Anchors reads by finding the anchor sequence and extracting sequence parts.

    params: 
        parts: tuples (part_len, part_name)
        anchor_seq: anchor sequence to find in the read.
        reverse_complement: if True, reverse complement the read before anchoring.

        part_anchor: name of the part to use an an anchor.

    """
    num_anchor_error = 0 
    num_reads = 0 
    part_names = [part_name for _, part_name in layout]
    part_lens = [part_len for part_len, _ in layout]
    try:
        anchor_idx = part_names.index(part_anchor)
    except ValueError:
        raise ValueError(f"Anchor part '{part_anchor}' not found in layout.")
    
    with open(out_csv, 'w') as out_f, open_fastq(fastq, 'r') as in_f:
         out_f.write(','.join(['read_id'] + part_names) + '\n')
         count=0
         for read in SeqIO.parse(in_f, 'fastq'): 
            num_reads += 1
            r = str(read.seq) if not reverse_complement else str(read.seq.reverse_complement())
            found_anchor_pos = r.find(anchor_seq)

            if found_anchor_pos == -1 or r.count(anchor_seq) > 1:
                num_anchor_error += 1
                continue
            part_seqs = decompose_read_into_parts(read.id, r, found_anchor_pos, anchor_idx, part_lens, part_names)
            if part_seqs is None: 
                num_anchor_error += 1
                continue
            out_f.write(','.join([read.id] + part_seqs) + '\n')
            # count += 1
            # if count >100000:  # for quick testing, remove later
            #     break
           

    print(f"Processed {sample} (reverse={reverse_complement}): reads {num_reads}, anchor errors {num_anchor_error}")

    frac = num_anchor_error / num_reads if num_reads > 0 else np.nan
    stats = pd.Series({
        'sample': sample, 
        'num_reads': num_reads,
        'frac_anchor_error': frac,
    })
    return stats

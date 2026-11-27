from Bio import SeqIO

"""
Reference sequence
"""
ref_genome_file = "ref/TnpB-omega_amplicon.fasta"
with open(ref_genome_file, "r") as ref_file:
    ref_record = next(SeqIO.parse(ref_file, "fasta"))

ref_seq = ref_record.seq
ref_aa_seq = ref_record.seq.translate()
ref_aa_seq_str = str(ref_aa_seq)

min_reads_th = 100 # all samples that have few reads will be discusrded during analysis.

"""
Output files where processed read-level data is stored. 
"""
from pathlib import Path 

fastq_dir = Path('data/fastqs_renamed')
def fastq_file(sample):
    return fastq_dir / (sample + ".fastq.gz")

sam_dir = Path('out/sam')
sam_dir.mkdir(exist_ok=True)

def sam_file(sample):
    return sam_dir / (sample + ".sam")

aligned_dir = Path('out/aligned')
aligned_dir.mkdir(exist_ok=True)

def aligned_file(sample):
    return aligned_dir / (sample + '.aligned.fasta')

base_dist_dir = Path('out/base_dist')
base_dist_dir.mkdir(exist_ok=True)

def base_dist_file(sample): 
    return base_dist_dir / (sample + '.base_dist.csv')

filtered_dir = Path('out/filtered')
filtered_dir.mkdir(exist_ok=True)

def filtered_subs_file(sample): 
    return filtered_dir / (sample + '.filter_subs.fasta')

gt_dir = Path('out/gt_dist')
gt_dir.mkdir(exist_ok=True)

def gt_dist_file(sample):
    return gt_dir / (sample + ".gt_dist.csv")

aa_dist_dir = Path('out/aa_dist')
aa_dist_dir.mkdir(exist_ok=True)

def aa_dist_file(sample):
    return aa_dist_dir / (sample + '.aa_dist.csv')

""" 
Paths where to store important results
"""

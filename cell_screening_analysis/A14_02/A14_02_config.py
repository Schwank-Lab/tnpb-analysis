
import os
import pandas as pd 
from pathlib import Path 
from Bio.Seq import Seq
from Bio.Data import IUPACData





# Read layout for this library.
# Each read part is taken at a fixed offset from a constant anchor sequence;
# R2 is reverse complemented first, so its parts are listed 5'->3' on that strand.
LIB_LAYOUT = {
    'r1_anchor': 'GGTTCAA',   # 3' end of the reRNA scaffold
    'r2_anchor': 'CTGAGTT',   # plasmid backbone, given on the forward strand
    'r1_layout': (
        (7, 'r1_anchor'),
        (20, 'spacer'),
    ),
    'r2_layout': (
        (20, 'target_site_rc'),
        (5, 'TAM_rc'),
        (7, 'barcode'),
        (7, 'r2_anchor'),
    ),
}
def expand_pam_pattern(pattern):
    iupac_dict = IUPACData.ambiguous_dna_values
    sequences = ['']

    for char in pattern:
        if char in iupac_dict:
            new_sequences = []
            for seq in sequences:
                for base in iupac_dict[char]:
                    new_sequences.append(seq + base)
            sequences = new_sequences
        else:
            sequences = [seq + char for seq in sequences]

    return sequences




LOW_READ_TH = 100  # minimum reads per target for a point to be used in any analysis
LOW_SITE_TH = 3 # mininum number of edit site in a site library to be considered a valid point in any analysis

PAMS = expand_pam_pattern('NNNNN')

# Data and output locations on the lab storage.
DATADIR = Path('/mnt/G_Pharma/Users/tfang/Cell_Screen/PW_A14-02_Screen')
WORKDIR = Path('/mnt/mnemo1/Cell_Screen/A14_02')  # local scratch; writing to the NAS is slow
if not (DATADIR / 'meta').is_dir():
    raise SystemExit(
        f"Sequencing data not found at {DATADIR}. "
        "Edit DATADIR at the top of this file to point at the folder "
        "holding meta/ and the fastq files."
    )

FASTQ_DIR = DATADIR / 'fastq'



def get_samples(): 
    samples=pd.read_excel(DATADIR / 'meta/A14-02_screen_novaseq_samplesheet_PW.xlsx',
                          index_col=False,engine='openpyxl')
    return samples

def get_library(target_rc_poitions=(136,161),barcode_positions=(161,168)):
    screen_libraries_path = Path(f'{DATADIR}/meta/A14-01_librarybuildoverview.xlsx')
    library = pd.read_excel(screen_libraries_path,index_col=None, dtype={'length': int})
    # print(library.head())
    library["barcode_verified"] = library["oligofinal"].apply(lambda x: x[barcode_positions[0]:barcode_positions[1]])
    library["target_site_rc_PlusTAM"] = library["oligofinal"].apply(lambda x: x[target_rc_poitions[0]:target_rc_poitions[1]])
    library["target_site_PlusTAM"] = library["target_site_rc_PlusTAM"].apply(lambda x: str(Seq(x).reverse_complement()))
    library["target_site"] = library["target_site_PlusTAM"].apply(lambda x: x[-20:])
    library["TAM"] = library["target_site_PlusTAM"].apply(lambda x: x[0:5])
    assert library["target_site"].equals(library["spacer"]), "The columns 'target_site' and 'spacer' are different."

    return library

def fastq_r1_file(fastq_file): 
    return Path(FASTQ_DIR / fastq_file)

def fastq_r2_file(fastq_r1_file):
    if "R1" in fastq_r1_file:
        return fastq_r1_file.parent / fastq_r1_file.name.replace('R1', 'R2')
    else:
        return Path(FASTQ_DIR / fastq_r1_file)
    




result_dir=WORKDIR  / 'out/results'
result_dir.mkdir(exist_ok=True, parents=True)


parts_dir = WORKDIR  / 'out/01_processed'
parts_dir.mkdir(exist_ok=True, parents=True)
def parts(sample):
    return parts_dir / (sample + '.parts.csv')
def parts_r1(sample): 
    return parts_dir / (sample + '.parts_R1.csv')
def parts_r2(sample): 
    return parts_dir / (sample + '.parts_R2.csv')

filtered_dir = WORKDIR / 'out/02_filtered'
filtered_dir.mkdir(exist_ok=True, parents=True)
def filtered(sample): 
    return filtered_dir / (sample + '.filtered.csv')

aligned_dir = WORKDIR / 'out/03_aligned'
aligned_dir.mkdir(exist_ok=True, parents=True)
def aligned(sample): 
    return aligned_dir / (sample + '.aligned.csv')

outcome_dir = WORKDIR / 'out/05_outcomes'
outcome_dir.mkdir(exist_ok=True, parents=True)
def outcomes(sample): 
    return outcome_dir / f'{sample}.outcomes.csv'

edit_rate_dir = WORKDIR / 'out/06_edit_rate'
edit_rate_dir.mkdir(exist_ok=True, parents=True)
def edit_rate(sample):
    return edit_rate_dir / f'{sample}.edit_rate.csv'

pam_profiles = WORKDIR / 'out/07_pam_profiles'
pam_profiles.mkdir(exist_ok=True, parents=True)
def pam_profile(sample, lib):
    return pam_profiles / f'{sample}_{lib}.pam_profile.png'


editing_windows = WORKDIR / 'out/08_editing_windows'
editing_windows.mkdir(exist_ok=True, parents=True)
def editing_window(editor, sublib):
    return editing_windows / f'{editor}_{sublib}.editing_window.png'


logoplots = WORKDIR / 'out/logoplots'
logoplots.mkdir(exist_ok=True, parents=True)
def logoplot(sample, lib):
    return logoplots / f'{sample}_{lib}.logoplot.png'


FASTQC_DIR = WORKDIR / 'out/fastqc'
FASTQC_DIR.mkdir(exist_ok=True, parents=True)
MULTIQC_DIR = WORKDIR / 'out/multiqc'
MULTIQC_DIR.mkdir(exist_ok=True, parents=True) 
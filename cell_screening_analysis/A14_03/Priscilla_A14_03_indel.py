import A14_03_config as config

import sys
import pathlib
import os
import pandas as pd
import logging  # Added for logging configuration
import seaborn as sns
import matplotlib.pyplot as plt
from prefect import flow, task, unmapped

# Configure logging to write to a file
from prefect.logging import get_logger
logger = get_logger()
file_handler = logging.FileHandler('pipeline.log')  # Log file name
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# Make the analysis packages importable without installing them.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from pipelines import common_tasks
from screening import screen_utils


@task
def process_reads(config,sample_name_column='sample_name'):
    samples = config.get_indel_samples()


    res = common_tasks.process_sample.map(samples[sample_name_column], samples.file_name_R1, unmapped(config), samples.file_name_R2)
    res = [res.result() for res in res]


@task
def filter_reads_per_sample(
    sample: str,
    lib_section: str,
    lib: pd.DataFrame,
    config,
    usecols: list = ['spacer', 'Random_filler_10bp_rc', 'barcode', 'TAM', 'target_site']
) -> pd.Series:
    """Filter out reads that don't match library members.

    Args:
        sample: Sample name
        lib_section: Library section identifier
        lib: Library dataframe with reference sequences
        config: Configuration object with paths and parameters
        usecols: List of columns to read from CSV. If None, uses default columns.

    Returns:
        Series containing filtering statistics for the sample
    """
    print(f"Filtering reads for sample: {sample}")


    # Load sample data with only necessary columns to save memory
    reads = pd.read_csv(
        config.parts(sample),
        index_col=None,
        usecols=usecols
    )
    num_reads = len(reads)

    # Calculate initial library membership fractions
    frac_spacer_from_lib = reads["spacer"].isin(lib['spacer'].unique()).sum() / num_reads
    frac_barcode_from_lib = reads["barcode"].isin(lib['barcode'].unique()).sum() / num_reads

    # Merge with library to retain only reads with known barcodes
    lib_cols = ['spacer', 'barcode', 'TAM', 'target_site']
    reads = pd.merge(
        reads,
        lib[lib_cols],
        left_on=['barcode'],
        right_on=['barcode'],
        how='inner',
        suffixes=('', '_lib')
    )

    num_reads_barcode_match = len(reads)
    frac_barcode_match = num_reads_barcode_match / num_reads

    # Calculate match statistics for reads with matching barcodes
    frac_spacer_match_when_barcode_match = (
        (reads["spacer"] == reads["spacer_lib"]).sum() / num_reads_barcode_match
    )

    frac_spacer_recombination_when_barcode_match = (
        (reads["spacer"] != reads["spacer_lib"]) &
        reads["spacer"].isin(lib['spacer'].unique())
    ).sum() / num_reads_barcode_match

    frac_TAM_match_when_barcode_match = (
        (reads["TAM"] == reads["TAM_lib"]).sum() / num_reads_barcode_match
    )

    # Calculate fraction with valid target positions (>= 18 matching positions)
    # Allow A in library to match C or G in reads
    def count_matching_positions(row):
        return sum(
            lib_base == read_base or (lib_base == 'A' and read_base in ['C', 'G'])
            for lib_base, read_base in zip(row['target_site_lib'], row['target_site'])
        )

    frac_valid_target_pos_when_barcode_match = (
        reads.apply(count_matching_positions, axis=1) >= 18
    ).sum() / num_reads_barcode_match

    # Filter to retain only reads with matching spacer and TAM
    reads = reads[
        (reads["spacer"] == reads["spacer_lib"]) &
        (reads["TAM"] == reads["TAM_lib"])
    ]

    # Save filtered reads
    reads.to_csv(config.filtered(sample))
    del reads

    return pd.Series({
        'sample_name': sample,
        'num_reads': num_reads,
        'frac_spacer_from_lib': frac_spacer_from_lib,
        'frac_barcode_from_lib': frac_barcode_from_lib,
        'num_reads_barcode_match': num_reads_barcode_match,
        'frac_barcode_match': frac_barcode_match,
        'frac_spacer_match_when_barcode_match': frac_spacer_match_when_barcode_match,
        'frac_spacer_recombination_when_barcode_match': frac_spacer_recombination_when_barcode_match,
        'frac_TAM_match_when_barcode_match': frac_TAM_match_when_barcode_match,
        'frac_valid_target_pos_when_barcode_match': frac_valid_target_pos_when_barcode_match
    })



@task(log_prints=True)
def filter_reads(config, sample_name_column='sample_name', is_indel=True, usecols=None):
    samples = config.get_indel_samples()
    lib = config.get_library()
    res = filter_reads_per_sample.map(samples[sample_name_column], samples["lib_section"], unmapped(lib), unmapped(config), unmapped(usecols))
    res = [res.result() for res in res]
    res = pd.DataFrame(res)
    if is_indel:
        res.to_excel(config.result_dir / 'indel_library_stats.xlsx', index=False)
    else:
        res.to_excel(config.result_dir / 'library_stats.xlsx', index=False)



#     d.loc[:, 'target_site_aligned'] = d.apply(lambda row: screen_utils.align_target_site(row.target_site, row.target_site_lib, insertion_side='five_prime'), axis=1)



@task(log_prints=True)
def count_outcomes_per_sample(sample,config,usecols=['Random_filler_10bp_rc','barcode','target_site', 'TAM', 'spacer','target_site_lib',]):
    """Count editing outcomes for a given sample."""
    print(f'Counting outcomes for {sample}')
    d_filterd=pd.read_csv(config.filtered(sample), index_col=0)
    # Aggregate reads sharing the same library member (barcode), edited target site and TAM.
    # 'spacer' and 'target_site_lib' are constant within a member and are carried through for the output.
    outcomes = d_filterd.groupby(usecols).size().reset_index(name='num_reads')
    outcomes.to_csv(config.outcomes(sample))

@task(log_prints=True)
def count_outcomes(config,sample_name_column='sample_name'):
    samples = config.get_indel_samples()
    res = count_outcomes_per_sample.map(samples[sample_name_column], unmapped(config))
    _ = [res.result() for res in res]



def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))


def build_tnpb_window(target_site, random_filler_10bp, tnpb_half_window_size, cut_site_position=15):
    """Build TNPB window centered around the cut site with configurable half-window size."""
    target_right_len = max(len(target_site) - cut_site_position, 0)
    left_start = max(cut_site_position - tnpb_half_window_size, 0)
    filler_bases_needed = max(tnpb_half_window_size - target_right_len, 0)
    return target_site[left_start:] + random_filler_10bp[:filler_bases_needed]

@task(log_prints=True,tags=["calculate_indel_edit_rate"])
def calculate_indel_edit_rate(sample, config, tnpb_half_window_size=10):
    print(f'Analyzing indel rates for {sample}')

    lib=config.get_library()
    lib["Random_filler_10bp_rc"]=lib["oligofinal"].str[126:126+10]
    lib["Random_filler_10bp"]=lib["Random_filler_10bp_rc"].apply(lambda x: reverse_complement(x))
    # https://www.nature.com/articles/s41586-021-04058-1: run-off sequencing of the cleavage products confirmed a staggered cleavage pattern at 15–21 bp from the TAM that resulted in 5′ overhangs (Fig. 2g).
    # left and right window around tnpb cut sites (default 10 bp each side), tnpb target site is chosen as 15 bp from TAM
    lib["tnpb_window"] = lib.apply(
        lambda row: build_tnpb_window(
            row["target_site"],
            row["Random_filler_10bp"],
            tnpb_half_window_size,
            cut_site_position=15,
        ),
        axis=1,
    )

    o = pd.read_csv(config.outcomes(sample), index_col=0)
    o["Random_filler_10bp"]=o["Random_filler_10bp_rc"].apply(lambda x: reverse_complement(x))
    # https://www.nature.com/articles/s41586-021-04058-1: run-off sequencing of the cleavage products confirmed a staggered cleavage pattern at 15–21 bp from the TAM that resulted in 5′ overhangs (Fig. 2g).
    # left and right window around tnpb cut sites (default 10 bp each side), tnpb target site is chosen as 15 bp from TAM
    o["tnpb_window"] = o.apply(
        lambda row: build_tnpb_window(
            row["target_site"],
            row["Random_filler_10bp"],
            tnpb_half_window_size,
            cut_site_position=15,
        ),
        axis=1,
    )
    o=o.merge(lib[["barcode","Random_filler_10bp","tnpb_window"]], left_on="barcode", right_on="barcode", how="left", suffixes=('', '_lib'))
    e1 = screen_utils.calculate_indel_rate_per_read_group(o, cols_groupby=['barcode','spacer',"TAM"])
    e1.to_csv(config.edit_rate(f'{sample}-reads-indelrates-unfiltered')) # here for quality check for each replicates


def calculate_indel_rates(config,sample_name_column='sample_name', tnpb_half_window_size=10): # tnpb_half_window_size=10
    indel_samples = config.get_indel_samples()
    res = calculate_indel_edit_rate.map(
        indel_samples[sample_name_column],
        unmapped(config),
        unmapped(tnpb_half_window_size),
    )
    res = [r.result() for r in res]


@task(log_prints=True)
def calculate_library_distribution(sample, config):
    """Calculate library member distribution for a sample."""
    print(f'Calculating library distribution for {sample}')
    lib = config.get_library()
    d = pd.read_csv(config.filtered(sample), index_col=0)

    # Count reads per barcode
    member_counts = d.groupby(['barcode']).size().reset_index(name='reads')

    # Merge with library info to get sublibrary information
    # Adjust column names based on your actual library structure
    lib_cols = ['barcode', 'spacer','TAM']

    member_counts = pd.merge(
        member_counts,
        lib[lib_cols],
        on='barcode',
        how='left'
    )

    member_counts['sample_name'] = sample
    return member_counts

@task(log_prints=True)
def collect_library_distributions(config, sample_name_column='sample_name'):
    """Collect library distributions for all samples."""
    print('Collecting library distributions for all samples')
    samples = config.get_samples()

    res = calculate_library_distribution.map(samples[sample_name_column], unmapped(config))
    res = [r.result() for r in res]

    library_stats = pd.concat(res, ignore_index=True)
    library_stats.to_excel(config.result_dir / 'library_member_stats.xlsx', index=False)
    return library_stats

@task(log_prints=True)
def plot_library_stats(config):
    """Plot library distribution histograms."""
    print('Plotting library distribution statistics')
    library_member_stats = pd.read_excel(config.result_dir / 'library_member_stats.xlsx')

    # Plot histograms
    g = sns.FacetGrid(library_member_stats, col='sample_name', height=4, aspect=1.5)

    g.map_dataframe(sns.histplot, x='reads', binwidth=10)
    g.set_titles('{col_name}')
    g.set_axis_labels('Number of reads', 'Count')

    plt.tight_layout()
    plt.savefig(config.result_dir / 'library_member_reads_dist.png', dpi=300)
    plt.close()

    print(f'Library distribution plot saved to {config.result_dir / "library_member_reads_dist.png"}')







@flow
def run_pipeline():
    print('Starting the pipeline...')
    # process_reads(config,sample_name_column='sample_name')
    # filter_reads(config,sample_name_column='sample_name')




    # count_outcomes(config,sample_name_column='sample_name')
    calculate_indel_rates(config,sample_name_column='sample_name')


    print('Pipeline completed.')

if __name__ == "__main__":
    run_pipeline()



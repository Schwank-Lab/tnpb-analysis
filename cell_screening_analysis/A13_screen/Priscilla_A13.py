import A13_config as config

import sys
import pathlib
import os
import pandas as pd
from prefect import flow, task, unmapped
from prefect_dask import DaskTaskRunner  # task runner used for these two memory-bound flows
import logging
from prefect.logging import get_logger

# Make the analysis packages importable without installing them.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from pipelines import common_tasks
from screening import screen_utils

# Configure logging to write to a file
logger = get_logger()
file_handler = logging.FileHandler('pipeline.log')  # Log file name
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

@task(log_prints=True)
def process_reads(config,sample_name_column='sample_name'):
    samples = config.get_samples()
    res = common_tasks.process_sample.map(samples[sample_name_column], samples.file_name_R1, unmapped(config), samples.file_name_R2)
    res = [res.result() for res in res]
    res = pd.DataFrame(res)
    res.to_excel(config.result_dir / 'read_stats.xlsx')

@task(log_prints=True)
def filter_reads_per_sample(sample: str, lib: pd.DataFrame, config) -> pd.Series:
    """ Filter out reads that don't match library members."""
    print(f"Filter out reads that don't match library members for {sample}")
    d = pd.read_csv(config.parts(sample), index_col=None, usecols=['spacer', 'barcode', 'TAM', 'target_site']) #use necessary columns to save memroy later
    num_reads = len(d.index)

    # Optimize data types for memory usage
    d['barcode'] = d['barcode'].astype('category')
    d['spacer'] = d['spacer'].astype('category')
    lib['barcode'] = lib['barcode'].astype('category')
    lib['spacer'] = lib['spacer'].astype('category')

    #frac_spacer_from_lib=sum(d["spacer"].isin(lib['spacer'].unique()))/num_reads
    frac_spacer_from_lib = sum(d["spacer"].str[:10].isin(lib['spacer'].str[:10].unique())) / num_reads
    frac_barcode_from_lib=sum(d["barcode"].isin(lib['barcode'].unique()))/num_reads


    #frac_barcode_match = d.barcode.isin(lib['barcode'].unique()).sum() / num_reads

    # retain reads with known sub-library and member barcodes
    #d = pd.merge(d, lib[['spacer', 'barcode','TAM_rc','TAM','target_site_rc','target_site']],
    d = pd.merge(d, lib[['spacer', 'barcode','TAM','target_site']],
                    left_on=['barcode', ],
                    right_on=['barcode'],
                    how='inner', suffixes=('', '_lib'))
    num_reads_barcode_match = len(d.index)
    frac_barcode_match = num_reads_barcode_match / num_reads

    # Spacers are compared on their first 10 bases, which is less sensitive to
    # sequencing error than the full spacer. Fractions are computed without
    # materialising intermediate columns.
    frac_spacer_match_when_barcode_match = (d["spacer"].str[:10] == d["spacer_lib"].str[:10]).sum() / num_reads_barcode_match
    frac_spacer_recombination_when_barcode_match = (
        (d["spacer"].str[:10] != d["spacer_lib"].str[:10]) &
        d["spacer"].str[:10].isin(lib['spacer'].str[:10].unique())
    ).sum() / num_reads_barcode_match
    frac_TAM_match_when_barcode_match = (d["TAM"] == d["TAM_lib"]).sum() / num_reads_barcode_match
    frac_valid_target_pos_when_barcode_match = (
        d.apply(lambda row: sum(tsl == ts or (tsl == 'A' and ts in ['C', 'G'])
                for tsl, ts in zip(row['target_site_lib'], row['target_site'])), axis=1) >= 18
    ).sum() / num_reads_barcode_match


    # retain reads with matching spacer
    #d = d[d.spacer_match & d.TAM_match]
    d = d[(d["spacer"].str[:10] == d["spacer_lib"].str[:10]) & (d["TAM"] == d["TAM_lib"])]
    # here as in A13 experiment, some spacers are 16nt and some are 20nt, so after check if spacer and spacer_lib matches in first 10 bases, replace it
    d['spacer'] = d["spacer_lib"]

    d.to_csv(config.filtered(sample))
    del d

    return pd.Series({
        'sample_name': sample, 'num_reads': num_reads,
        'frac_spacer_from_lib':frac_spacer_from_lib,'frac_barcode_from_lib':frac_barcode_from_lib,
        'num_reads_barcode_match': num_reads_barcode_match,
        'frac_barcode_match': frac_barcode_match,
        'frac_spacer_match_when_barcode_match': frac_spacer_match_when_barcode_match, 'frac_spacer_recombination_when_barcode_match': frac_spacer_recombination_when_barcode_match,
        'frac_TAM_match_when_barcode_match': frac_TAM_match_when_barcode_match,
        'frac_valid_target_pos_when_barcode_match':frac_valid_target_pos_when_barcode_match
    })

@task(log_prints=True)
def filter_reads(config,sample_name_column='sample_name'):
    samples = config.get_samples()
    lib = config.get_library()
    res = filter_reads_per_sample.map(samples[sample_name_column], unmapped(lib), unmapped(config))
    res = [res.result() for res in res]
    res = pd.DataFrame(res)
    res.to_excel(config.result_dir / 'library_stats.xlsx')



@task(log_prints=True)
def align_target_sites(sample: str, config):
    print(f'Aligning target sites for {sample}')
    d = pd.read_csv(config.filtered(sample), index_col=0)
    # PAM is on the 3' end of the target site, so replace bases with <isnertion> tag on the 5' end to make sure to not disrupt the PAM sequence.
    d.loc[:, 'target_site_aligned'] = d.apply(lambda row: screen_utils.align_target_site(row.target_site, row.target_site_lib, insertion_side='five_prime'), axis=1)
    d.to_csv(config.aligned(sample))

@task(log_prints=True)
def process_target_sites(config,sample_name_column='sample_name'):
    samples = config.get_samples()
    res = align_target_sites.map(samples[sample_name_column], unmapped(config))
    _ = [res.result() for res in res]


@task(log_prints=True)
def count_outcomes_per_sample(sample,config):
    """Count editing outcomes for a given sample."""
    print(f'Counting outcomes for {sample}')
    d_filterd=pd.read_csv(config.filtered(sample), index_col=0)
    # Aggregate reads sharing the same library member (barcode), edited target site and TAM.
    # 'spacer' and 'target_site_lib' are constant within a member and are carried through for the output.
    outcomes = d_filterd.groupby(['barcode','target_site', 'TAM', 'spacer','spacer_lib','target_site_lib', ]).size().reset_index(name='num_reads')
    outcomes.to_csv(config.outcomes(sample))

@task(log_prints=True)
def count_outcomes(config,sample_name_column='sample_name'):
    samples = config.get_samples()
    res = count_outcomes_per_sample.map(samples[sample_name_column], unmapped(config))
    _ = [res.result() for res in res]



@task(log_prints=True,tags=["calculate_abe_edit_rate"])
def calculate_abe_edit_rate(sample, config):
    print(f'Analyzing bystander rates for {sample}')
    o = pd.read_csv(config.outcomes(sample), index_col=0)
    e1 = screen_utils.calculate_edit_rate_per_read_group(o.copy(),  edited_spacer="target_site",col_ref_spacer='target_site_lib', cols_groupby=['barcode',"TAM"])
    e1.to_csv(config.edit_rate(f'{sample}-reads-editedrates-unfiltered')) # here for quality check for each replicates

    o = screen_utils.extract_abe_edit_positions(o, col_ref_spacer='target_site_lib',col_edit_position='edit_position') # spacer are/should same as target_site
    o = o[o['edit_position'].notna()]# whend do ABE; not all libraries could have A in thier target site
    o['is_target_base_edited'] = o.apply(lambda row: screen_utils.is_ontarget_a_to_g_conversion(row["target_site"], row["target_site_lib"], row["edit_position"]), axis=1)
    e = screen_utils.calculate_edit_rate_per_group(o, col_is_outcome_edited='is_target_base_edited', cols_groupby=['barcode',"TAM", 'edit_position',])# here calculate for each library at different edit_positions.
    e.to_csv(config.edit_rate(f'{sample}-unfiltered'))
    #e.to_csv(config.edit_rate(f'{sample}-unfiltered', sublib))

def get_abe_samples(config,sample_name_column='sample_name',editor_column='editor'):
    samples = config.get_samples()
    return samples.loc[samples[editor_column].str.contains('ABE'),:]
    #return samples.loc[samples[editor_column].str.contains('ABE'),sample_name_column]
    # controls are shared between ABE and nucleases.



# async with get_client() as client:


def calculate_abe_edit_rates(config,sample_name_column='sample_name',editor_column='editor'):
    abe_samples = get_abe_samples(config,sample_name_column,editor_column)
    res = calculate_abe_edit_rate.map(abe_samples[sample_name_column], unmapped(config))
    res = [r.result() for r in res]

@task(log_prints=True)
def filter_edit_rate(sample, config):
    print(f'Filtering edit rates for {sample}')
    d = pd.read_csv(config.edit_rate(f'{sample}-unfiltered'), index_col=0)
    # filter out low reads.
    d.loc[d.num_reads < config.LOW_READ_TH, 'edit_rate'] = pd.NA
    frac_low_reads = d['edit_rate'].isna().sum() / len(d.index)
    d.to_csv(config.edit_rate(sample))
    return pd.Series({'sample_name': sample,  'frac_low_reads_per_edit_postion': frac_low_reads})

@task(log_prints=True)
def filter_edit_rates(config,sample_name_column='sample_name'):
    samples = config.get_samples()
    res = filter_edit_rate.map(samples[sample_name_column], unmapped(config))
    res = [res.result() for res in res]
    res = pd.DataFrame(res)
    res.to_excel(config.result_dir / 'filter_edit_rate_stats.xlsx')

@task(log_prints=True)
def collect_abe_edit_rates(config,sample_name_column='sample_name',editor_column='editor',control_column='control_yes_no'):
    """ Collect edit rates for all samples"""
    print('Collecting ABE edit rates for all samples')
    samples = get_abe_samples(config,sample_name_column,editor_column)
    all_sample_rates = pd.concat(
        [
            pd.read_csv(config.edit_rate(f'{sample[sample_name_column]}'), index_col=0).assign(
                editor=sample[editor_column],
                sample_name=sample[sample_name_column],
                is_ctrl=sample[control_column],
                replicate=sample["replicate"],
            )
            for _,sample in samples.iterrows()
        ]
    ).reset_index(drop=True)
    all_sample_rates.loc[:, 'is_ctrl'] = all_sample_rates.is_ctrl == 'yes'
    all_sample_rates.to_csv(config.edit_rate('all-abe-samples-not-normalized'))

@task(log_prints=True)
def normalize_abe_editing_rate(config):
    d = pd.read_csv(config.edit_rate('all-abe-samples-not-normalized'), index_col=0)
    d,d_final = screen_utils.normalize_editing_rate(d,
        # Normalising within barcode/TAM/edit_position also absorbs target-site recombination.
        col_normalize_over=['barcode','TAM','edit_position'],
        col_rate='edit_rate',
        col_control='is_ctrl',
        # A member failing at any edit position is dropped for all of them.
        ignore_control=['barcode', 'TAM'],
        th_ctrl_missing=0.05,
        ctrl_edit_rate_thres=0.05
    )
    d.drop(columns=['is_ctrl'])
    d.to_csv(config.edit_rate('all-abe-samples-normalized'))
    d_final.to_csv(config.edit_rate('all-abe-samples-normalized-PlusCtrl'))


@task(log_prints=True)
def aggregate_abe_replicates(config):
    print("aggregating ABE replicates")
    d = pd.read_csv(config.edit_rate('all-abe-samples-normalized'), index_col=0)
    d=d.loc[~d["edit_rate"].isna(),:]
    # Group by the specified columns and filter groups with exactly all replicates
    d = d.groupby(['editor', 'barcode', 'TAM', 'edit_position']).filter(lambda x: len(x) == 3)
    d = screen_utils.aggregate_replicates(d,
        cols_group=['editor','barcode','TAM','edit_position'], # here editor and barcoe(library) combination is a biological replicate, PAM/TAM is not consider here ?
        col_rate='edit_rate'                             # PAM/TAM is not consider here ?
    ).drop(columns=['sample_name', 'replicate'])
    d.to_csv(config.edit_rate('all-abe-editors'))

    d_plusctrl = pd.read_csv(config.edit_rate('all-abe-samples-normalized-PlusCtrl'), index_col=0)
    d_plusctrl=d_plusctrl.loc[~d_plusctrl["edit_rate"].isna(),:]
    # Group by the specified columns and filter groups with exactly all replicates
    d_plusctrl = d_plusctrl.groupby(['editor', 'barcode', 'TAM', 'edit_position']).filter(lambda x: len(x) == 3)
    d_plusctrl = screen_utils.aggregate_replicates(d_plusctrl,
        cols_group=['editor','barcode','TAM','edit_position'], # here editor and barcoe(library) combination is a biological replicate, PAM/TAM is not consider here ?
        col_rate='edit_rate'                             # PAM/TAM is not consider here ?
    ).drop(columns=['sample_name', 'replicate'])



    d_plusctrl.to_csv(config.edit_rate('all-abe-editors-PlusCtrl'))






@flow(log_prints=True,task_runner=DaskTaskRunner(cluster_kwargs={"n_workers": 6,"memory_limit": "15GB"}))
def run_pipeline():
    print('Starting the pipeline...')
    # process_reads(config,sample_name_column='long_name')
    # filter_reads(config,sample_name_column='long_name')
    # count_outcomes(config,sample_name_column='long_name')
    calculate_abe_edit_rates(config,sample_name_column='long_name',editor_column='editor')
    filter_edit_rates(config,sample_name_column='long_name')
    collect_abe_edit_rates(config,sample_name_column='long_name',editor_column='editor')
    normalize_abe_editing_rate(config)
    aggregate_abe_replicates(config)
    print('Pipeline completed.')

if __name__ == "__main__":
    run_pipeline()

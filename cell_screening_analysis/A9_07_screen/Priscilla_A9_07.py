import A9_07_config as config

import sys
import pathlib
import os
import pandas as pd
import logging  # Added for logging configuration
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
    samples = config.get_samples()

    skipped_samples = samples[samples[sample_name_column].apply(lambda sample: (config.result_dir / f'{sample}_read_stats.csv').exists())]
    print(f"Skipped samples: {skipped_samples[sample_name_column].tolist()}")
    samples = samples[~samples[sample_name_column].apply(lambda sample: (config.result_dir / f'{sample}_read_stats.csv').exists())]

    res = common_tasks.process_sample.map(samples[sample_name_column], samples.file_name_R1, unmapped(config), samples.file_name_R2)
    res = [res.result() for res in res]

@task
def filter_reads_per_sample(sample: str,lib_section: str, lib: pd.DataFrame, config) -> pd.Series:


    """ Filter out reads that don't match library members."""
    d = pd.read_csv(config.parts(sample), index_col=None) #use necessary columns to save memroy later ,nrows=1000
    num_reads = len(d.index)


    d["spacer_from_lib"]=d["spacer"].str[:10].isin(lib.loc[lib["lib_section"]==lib_section,'spacer'].str[:10].unique()) # in A09, there is only one spacer per sample
    frac_spacer_from_lib = d["spacer_from_lib"].sum() / num_reads
    d=d[d["spacer_from_lib"]]


    d = pd.merge(d, lib[['spacer', 'target_site']],
                    left_on=['spacer', ],    # here as barcode is missing so use spacer as barcode instead
                    right_on=['spacer'],
                    how='inner', suffixes=('', '_lib'))
    num_reads_barcode_match = len(d.index)
    frac_barcode_match = num_reads_barcode_match / num_reads


    frac_valid_target_pos_when_barcode_match = (
        d.apply(lambda row: sum(tsl == ts or (tsl == 'A' and ts in ['C', 'G'])
                for tsl, ts in zip(row['target_site_lib'], row['target_site'])), axis=1) >= 18
    ).sum() / num_reads_barcode_match

    d.to_csv(config.filtered(sample))
    del d

    return pd.Series({
        'sample_name': sample, 'num_reads': num_reads,
        'frac_spacer_from_lib':frac_spacer_from_lib,
        'num_reads_barcode_match': num_reads_barcode_match,
        'frac_barcode_match': frac_barcode_match,
        'frac_valid_target_pos_when_barcode_match':frac_valid_target_pos_when_barcode_match
    })


@task(log_prints=True)
def filter_reads(config,sample_name_column='sample_name'):
    samples = config.get_samples()
    lib = config.get_library()
    res = filter_reads_per_sample.map(samples[sample_name_column],samples["lib_section"], unmapped(lib), unmapped(config))
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
    outcomes = d_filterd.groupby(['target_site', 'TAM', 'spacer','target_site_lib',]).size().reset_index(name='num_reads')
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
    e1 = screen_utils.calculate_edit_rate_per_read_group(o.copy(),  edited_spacer="target_site",col_ref_spacer='target_site_lib', cols_groupby=['spacer',"TAM"],editor="tnpb")
    e1.to_csv(config.edit_rate(f'{sample}-reads-editedrates-unfiltered')) # here for quality check for each replicates

    o = screen_utils.extract_abe_edit_positions(o, col_ref_spacer='target_site_lib',col_edit_position='edit_position') # spacer are/should same as target_site
    o = o[o['edit_position'].notna()]# when do ABE; not all libraries could have A in thier target site
    o['is_target_base_edited'] = o.apply(lambda row: screen_utils.is_ontarget_a_to_g_conversion(row["target_site"], row["target_site_lib"], row["edit_position"]), axis=1)
    e = screen_utils.calculate_edit_rate_per_group(o, col_is_outcome_edited='is_target_base_edited', cols_groupby=['spacer',"TAM", 'edit_position',])# here calculate for each library at different edit_positions.
    e.to_csv(config.edit_rate(f'{sample}-unfiltered'))
    #e.to_csv(config.edit_rate(f'{sample}-unfiltered', sublib))

def get_abe_samples(config,sample_name_column='sample_name',editor_column='editor'):
    samples = config.get_samples()
    return samples.loc[samples[editor_column].str.contains('ABE'),:]
    #return samples.loc[samples[editor_column].str.contains('ABE'),sample_name_column]
    # controls are shared between ABE and nucleases.



# async with get_client() as client:


def calculate_abe_edit_rates(config,sample_name_column='sample_name',editor_column='editor'):
    #abe_samples = get_abe_samples(config,sample_name_column,editor_column)
    abe_samples = config.get_samples()
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
    #samples = get_abe_samples(config,sample_name_column,editor_column)
    samples = config.get_samples()
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









@flow
def run_pipeline():
    print('Starting the pipeline...')
    process_reads(config,sample_name_column='sample_name')
    filter_reads(config,sample_name_column='sample_name')
    count_outcomes(config,sample_name_column='sample_name')
    calculate_abe_edit_rates(config,sample_name_column='sample_name',editor_column='editor')
    filter_edit_rates(config,sample_name_column='sample_name')
    collect_abe_edit_rates(config,sample_name_column='sample_name',editor_column='editor',control_column='is_ctrl')

    print('Pipeline completed.')

if __name__ == "__main__":
    run_pipeline()



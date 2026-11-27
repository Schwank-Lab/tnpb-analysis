"""Per-sample read processing task. Extracted from the Cell_Screen pipelines package; only the functions used by the notebooks in this repository are kept."""

import pandas as pd
from Bio.Seq import Seq
from prefect.tasks import task
from sequencing import sequencing


@task(log_prints=True)
def process_sample(sample, r1_fastq, config, r2_fastq=None,TAM_tag="TAM"):
    print(f'Processing {sample}')
    """Anchor reads and extract read parts."""
    r1_fastq = config.fastq_r1_file(r1_fastq)
    assert r1_fastq.exists(), f'R1 fastq file {r1_fastq} does not exist'
    print(r1_fastq)
    if not r2_fastq: 
        r2_fastq = config.fastq_r2_file(r1_fastq)
    else: 
        r2_fastq = config.fastq_r2_file(r2_fastq)
    assert r2_fastq.exists(), f'R2 fastq file {r2_fastq} does not exist'
    
    sequencing.run_fastqc(r1_fastq, config.FASTQC_DIR)
    sequencing.run_fastqc(r2_fastq, config.FASTQC_DIR)
    
    r1_stats = sequencing.anchor_reads(sample, config.LIB_LAYOUT['r1_layout'], config.LIB_LAYOUT['r1_anchor'], 
        r1_fastq, config.parts_r1(sample), part_anchor='r1_anchor'
    )
    r2_stats = sequencing.anchor_reads(sample, config.LIB_LAYOUT['r2_layout'], config.LIB_LAYOUT['r2_anchor'],  
        r2_fastq, config.parts_r2(sample), 
        reverse_complement=True, part_anchor='r2_anchor')
    
    
    d = pd.merge(pd.read_csv(config.parts_r1(sample), index_col=0), pd.read_csv(config.parts_r2(sample), index_col=0), left_index=True, right_index=True)
    
    # the whole target site is added to the amplicon as a reverse complement, so reverse it back.
    d.loc[:, 'target_site'] = d.target_site_rc.map(lambda x: str(Seq(x).reverse_complement()))
    
    if TAM_tag:
        d.loc[:, 'TAM'] = d.TAM_rc.map(lambda x: str(Seq(x).reverse_complement()))
    

    
    d.to_csv(config.parts(sample))
    
    stats = pd.Series({
        'sample_name': sample, 
        'num_reads': r1_stats['num_reads'], 
        'r1_frac_anchor_error': r1_stats['frac_anchor_error'], 
        'r2_frac_anchor_error': r2_stats['frac_anchor_error'],
    })
    stats.to_csv(config.result_dir / f'{sample}_read_stats.csv')

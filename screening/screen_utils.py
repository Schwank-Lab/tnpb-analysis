"""Screening analysis helpers. Extracted from the Cell_Screen screening package; only the functions used by the notebooks in this repository are kept."""

import pandas as pd
import parasail

try:
    from Bio import pairwise2
except ImportError:  # pairwise2 is deprecated; only the alignment indel caller needs it
    pairwise2 = None


def align_target_site(target_site, ref_site, insertion_side='five_prime', mismatch_penalty=1, gap_penalty=1):
    """Align a target site to a reference site using parasail for speed.
    
    Returns a string representing the alignment, marking insertions with 'I'.
    """
    assert len(ref_site) > 0, 'Reference spacer is empty'
    matrix = parasail.matrix_create("ACGT", 1, -mismatch_penalty)
    # Global Needleman-Wunsch alignment with traceback; note gap penalties are passed as positive values.
    result = parasail.nw_trace_striped_16(target_site, ref_site, gap_penalty, gap_penalty, matrix)
    aligned_target = result.traceback.query
    aligned_ref = result.traceback.ref

    aligned_site = []
    insertion_flag = [False] * len(ref_site)
    aligned_index = 0
    for r, t in zip(aligned_ref, aligned_target):
        if r != '-':
            aligned_site.append(t)
            aligned_index += 1
        else:
            if insertion_side == 'five_prime':
                if aligned_index > 0:
                    insertion_flag[aligned_index - 1] = True 
            else:  # insertion_side == 'three_prime'
                if aligned_index < len(ref_site):
                    insertion_flag[aligned_index] = True
    final_alignment = ''.join('I' if flag and base != '-' else base for base, flag in zip(aligned_site, insertion_flag))
    return final_alignment


def is_indel_by_alignment(
    edited_seq,
    ref_seq,
    min_homology=0.7,
    ignore_substitutions=True,
    mismatch_penalty=1,
    gap_penalty=1, # when change to 2, indel becomes too low 
):
    """Determine whether an edited sequence should be counted as an indel event.

    This function is used for stricter indel calling than simple string comparison.
    It aligns ``edited_seq`` to ``ref_seq`` with Biopython pairwise alignment,
    filters low-quality alignments by a minimum homology threshold, and then
    decides whether the read is "indeled".

    Decision logic:
    1) Align ``ref_seq`` and ``edited_seq`` globally.
    2) Compute homology = matches / aligned reference length (excluding gaps).
    3) If homology < ``min_homology``, return ``False``.
    4) If ``ignore_substitutions=True``, return ``True`` only when the alignment
       contains at least one gap event (insertion or deletion).
    5) Otherwise, return ``edited_seq != ref_seq`` (counts substitutions too).

    Args:
        edited_seq (str): Observed sequence from read/outcome.
        ref_seq (str): Reference sequence to compare against.
        min_homology (float): Minimum alignment homology required to evaluate
            the read as a valid comparison (default 0.7).
        ignore_substitutions (bool): If True, mismatch-only changes are ignored
            to reduce false positives; only insertion/deletion events are counted.
        mismatch_penalty (int): Mismatch penalty used in pairwise alignment.
        gap_penalty (int): Gap open/extend penalty used in pairwise alignment.

    pairwise2.align.globalms settings used here:
        - seqA = ref_seq
        - seqB = edited_seq
        - match score = 1
        - mismatch score = -mismatch_penalty (default: -1)
        - gap open score = -gap_penalty (default: -1)
        - gap extend score = -gap_penalty (default: -1)
        - one_alignment_only = True

    Biopython default optional behavior (if not overridden):
        - one_alignment_only = False
        - score_only = False
        - penalize_extend_when_opening = False
        - penalize_end_gaps = True for global alignments

    Returns:
        bool: True if the read is classified as an indel under the selected
        strictness settings; otherwise False.
    """
    edited_seq = str(edited_seq)
    ref_seq = str(ref_seq)


    # pairwise2 returns aligned sequences with '-' for gaps.
    alignments = pairwise2.align.globalms(
        ref_seq,
        edited_seq,
        1,
        -mismatch_penalty,
        -gap_penalty,
        -gap_penalty,
        one_alignment_only=True,
    )
    if not alignments:
        return False
    aligned_ref, aligned_target = alignments[0].seqA, alignments[0].seqB

    matches = sum(1 for r, t in zip(aligned_ref, aligned_target) if r == t and r != '-')
    alignment_length = len([r for r in aligned_ref if r != '-'])
    homology = matches / alignment_length if alignment_length > 0 else 0
    if homology < min_homology:
        return False

    if ignore_substitutions:
        # Only count insertion/deletion events as edited.
        return any((r == '-' and t != '-') or (t == '-' and r != '-') for r, t in zip(aligned_ref, aligned_target))

    return edited_seq != ref_seq


def is_indel_by_hamming_threshold(
    edited_seq,
    ref_seq,
    similarity_threshold=0.9,
):
    """Classify as indel using Hamming similarity for equal-length sequences.

    If similarity >= similarity_threshold, return False (not indel);
    otherwise return True (indel).
    """
    edited_seq = str(edited_seq)
    ref_seq = str(ref_seq)

    if len(edited_seq) != len(ref_seq):
        raise ValueError("edited_seq and ref_seq must have the same length for Hamming distance")

    if len(ref_seq) == 0:
        return False

    mismatches = sum(1 for a, b in zip(edited_seq, ref_seq) if a != b)
    similarity = 1 - (mismatches / len(ref_seq))
    return similarity < similarity_threshold


def calculate_indel_rate_per_read_group(
    df,
    cols_groupby=['barcode', "TAM"],
    col_edit_rate='frac_reads_edited',
    similarity_threshold=0.85,
    indel_calling_method='hamming',
):
    """Calculate the fraction of reads carrying an indel, per read group.

    A read is scored as indel-containing when its window around the TnpB
    cleavage site ('tnpb_window') shares less than `similarity_threshold`
    identity with the designed window ('tnpb_window_lib'). At the default
    0.85 over a 20-bp window this means four or more mismatching positions.

    The window is extracted at a fixed offset from the constant backbone, so
    an insertion or deletion at the cleavage site displaces sequence within it.
    """
    if indel_calling_method == 'hamming':
        df['read_indeled'] = df.apply(
            lambda row: is_indel_by_hamming_threshold(
                row["tnpb_window"],
                row["tnpb_window_lib"],
                similarity_threshold=similarity_threshold,
            ),
            axis=1,
        )
    elif indel_calling_method == 'alignment':
        # Align the window instead of comparing position by position, and count only
        # gap events. `similarity_threshold` is then the minimum alignment homology.
        df['read_indeled'] = df.apply(
            lambda row: is_indel_by_alignment(
                row["tnpb_window"],
                row["tnpb_window_lib"],
                min_homology=similarity_threshold,
                ignore_substitutions=True,
            ),
            axis=1,
        )
    else:
        raise ValueError("indel_calling_method must be 'hamming' or 'alignment'")

    def calculate_edit_rate(group):
        total_reads = group['num_reads'].sum()
        edited_reads = group.loc[group['read_indeled'], 'num_reads'].sum()
        return edited_reads / total_reads
    
    # Calculate edit rate for each group
    df_grouped = df.groupby(cols_groupby).apply(
        lambda group: pd.Series({
            col_edit_rate: calculate_edit_rate(group),
            'num_reads': group['num_reads'].sum(),
            **{col: group[col].iloc[0] for col in df.columns if col not in cols_groupby + ['read_indeled', 'num_reads']}
        }),
        include_groups=False
    ).reset_index()
    
    return df_grouped


def calculate_edit_rate_per_read_group(df, edited_spacer="target_site",col_ref_spacer='target_site_lib',cols_groupby=['barcode',"TAM"],col_edit_rate = 'frac_reads_edited',editor="abe"):
    """ check for each target, how many of reads being edited. one read count once only  
    """ 


    df['edit_positions'] = df[col_ref_spacer].apply(lambda s: [i+1 for i, b in enumerate(s) if b == "A"])
    if editor=="abe":
        df['read_edited'] = df.apply(lambda row: any(row[edited_spacer][i-1] == "G" for i in row['edit_positions']), axis=1)
    elif editor=="tnpb":
        df['read_edited'] = df.apply(lambda row: any(2 <= i <= 12 and row[edited_spacer][i-1] == "G" for i in row['edit_positions']), axis=1) # for tnpb, the editing windows is 3-10
    
    
    def calculate_edit_rate(group):
        total_reads = group['num_reads'].sum()
        edited_reads = group.loc[group['read_edited'], 'num_reads'].sum()
        return edited_reads / total_reads
    
    # Calculate edit rate for each group                                                            
    df_grouped = df.groupby(cols_groupby).apply(
        lambda group: pd.Series({
            col_edit_rate: calculate_edit_rate(group),
            'num_reads': group['num_reads'].sum(),
            **{col: group[col].iloc[0] for col in df.columns if col not in cols_groupby + ['read_edited', 'num_reads']}
        }),
        include_groups=False
    ).reset_index()
    
    return df_grouped


def extract_abe_edit_positions(df, col_ref_spacer='ref_spacer', col_edit_position='edit_position'):
    """ Extract all A->G editing positions from df. 
    """ 
    df[col_edit_position] = df[col_ref_spacer].apply(lambda s: [i + 1 for i, b in enumerate(s) if b == "A"])
    return df.explode(col_edit_position)


def is_ontarget_a_to_g_conversion(edited_spacer, ref_spacer, position):
    expected_base = ref_spacer[position-1]
    
    if not expected_base == 'A':
        raise ValueError(f'Expected base at position {position} is not A. Got {expected_base}')
    target_base = edited_spacer[position-1]
    return target_base == 'G'


def calculate_edit_rate_per_group(outcomes: pd.DataFrame, col_is_outcome_edited: str, cols_groupby: list[str], 
                                    col_num_reads: str = 'num_reads', col_edit_rate: str = 'edit_rate') -> pd.DataFrame:
    """ Aggregate outcomes into edit rates.
    
    Args:
        col_is_outcome_edited (str): column name for outcome flag. If true, outcomes is counted as edited.
        cols_groupby (List[str]): columns by which to group outcomes together.
        col_num_reads (str): column name for number of reads
        col_edit_rate (str): column name for editing rate
    Returns:
        pd.DataFrame: df with columns
            - col_num_reads: number of total reads per group
            - col_edit_rate: editing rate per group 
            - rest of the columns in ooutocmes
    """
    
    # Outcomes are read groups, so the edited count is (number of reads) x (outcome is edited).
    outcomes['num_outcome_is_edited'] = outcomes.apply(
        lambda row: row[col_num_reads] if row[col_is_outcome_edited] else 0, axis=1
    )
    
    
    #for each cols_groupby, get the the most common target_site that are different as the target_site_lib
    idx = outcomes[outcomes['target_site'] != outcomes['target_site_lib']].groupby(cols_groupby)['num_reads'].idxmax()
    target_site_map = outcomes.loc[idx].set_index(cols_groupby)['target_site'].to_dict()
    target_site_num_map = outcomes.loc[idx].set_index(cols_groupby)['num_reads'].to_dict()
    outcomes['most_common_diff_target_site'] = outcomes[cols_groupby].apply(lambda row: target_site_map.get(tuple(row), None), axis=1)
    outcomes['most_common_diff_target_site_num'] = outcomes[cols_groupby].apply(lambda row: target_site_num_map.get(tuple(row), None), axis=1)
    
    # Within a group, 'target_site' retained below is one of several edited versions;
    # editing is scored per edit_position, so the representative sequence does not affect the rate.
    d_grouped = outcomes.groupby(cols_groupby).agg({
        col_num_reads: 'sum',
        'num_outcome_is_edited': 'sum',
        **{col: 'first' for col in outcomes.columns if col not in cols_groupby + [col_num_reads, 'num_outcome_is_edited']}
    }).reset_index()
    
    
    d_grouped[col_edit_rate] = d_grouped.num_outcome_is_edited / d_grouped[col_num_reads]
    d_grouped = d_grouped.drop(columns=['num_outcome_is_edited', col_is_outcome_edited])
    return d_grouped


def _drop_by_control(d_ctrl, d_no_ctrl, mask, ignore_control, label):
    """Drop every record whose ignore_control key appears in the masked control rows."""
    bad = d_ctrl[mask][ignore_control]
    n_before = d_ctrl.shape[0]
    key = lambda df: df[ignore_control].apply(tuple, axis=1)
    bad_keys = bad.apply(tuple, axis=1)
    d_ctrl = d_ctrl[~key(d_ctrl).isin(bad_keys)]
    d_no_ctrl = d_no_ctrl[~key(d_no_ctrl).isin(bad_keys)]
    print(f'{label}: {bad.shape[0] / n_before if n_before else 0}. Drop these records')
    return d_ctrl, d_no_ctrl


def _normalize(data, col_normalize_over, col_rate, col_control, ignore_control,
               th_ctrl_missing, ctrl_edit_rate_thres, on_missing_control):
    """Subtract the background editing measured in editor-free control samples.

    Library members whose mean control rate exceeds `ctrl_edit_rate_thres`, and
    members whose control value is missing because of too few reads, are dropped
    rather than corrected. For the rest the mean control rate is subtracted and
    the difference divided by one minus that rate, so the corrected value gives
    editing as a fraction of what was still attainable above background.
    Negative values are set to zero.

    on_missing_control: 'raise' aborts when more than `th_ctrl_missing` of sites
    have no control at all; 'warn' reports and continues.
    """
    d_ctrl_ori = data[data[col_control]].copy()
    d_ctrl = data[data[col_control]][col_normalize_over + [col_rate]].rename(
        columns={col_rate: 'edit_rate_ctrl'})
    d_no_ctrl = data[~data[col_control].astype(bool)].copy()

    d_ctrl, d_no_ctrl = _drop_by_control(
        d_ctrl, d_no_ctrl, d_ctrl['edit_rate_ctrl'].isna(), ignore_control,
        'frac_low_num_read_in_control')
    d_ctrl, d_no_ctrl = _drop_by_control(
        d_ctrl, d_no_ctrl, d_ctrl['edit_rate_ctrl'] > ctrl_edit_rate_thres, ignore_control,
        'frac_unintendededit_in_control')

    # there may be several control points per group
    d_ctrl = d_ctrl.groupby(col_normalize_over).edit_rate_ctrl.mean().reset_index()

    d_norm = pd.merge(d_no_ctrl, d_ctrl, on=col_normalize_over, how='left')
    ctrl_missing = d_norm[d_norm['edit_rate_ctrl'].isna()].groupby(col_normalize_over).first()
    n_missing = len(ctrl_missing)
    n_total = len(d_ctrl) + n_missing
    if n_total and n_missing / n_total > th_ctrl_missing:
        msg = f'{n_missing} out of {n_total} sites are missing control editing rate.'
        if on_missing_control == 'raise':
            raise ValueError(msg)
        print(msg)
    else:
        print(f'{n_missing} out of {n_total} sites are missing control editing rate. drop these sites')
        d_norm = d_norm.dropna(subset=['edit_rate_ctrl'])

    d_norm[col_rate] = (d_norm[col_rate] - d_norm['edit_rate_ctrl']) / (1 - d_norm['edit_rate_ctrl'])
    d_norm[col_rate] = d_norm[col_rate].clip(lower=0)
    d_norm = d_norm.drop(columns=['edit_rate_ctrl'])
    return d_norm, pd.concat([d_norm, d_ctrl_ori], ignore_index=True)


def normalize_editing_rate(data, col_normalize_over=['site_id', 'sublibrary'], col_rate='edit_rate',
                           col_control='is_ctrl', ignore_control=['barcode', 'TAM'],
                           th_ctrl_missing=0.05, ctrl_edit_rate_thres=0.05,
                           on_missing_control='raise') -> pd.DataFrame:
    """Normalize per-position editing rates against editor-free controls."""
    return _normalize(data, col_normalize_over, col_rate, col_control, ignore_control,
                      th_ctrl_missing, ctrl_edit_rate_thres, on_missing_control)


def normalize_editing_rate_reads(data, col_normalize_over=['site_id', 'sublibrary'],
                                 col_rate='edit_rate', col_control='is_ctrl',
                                 th_ctrl_missing=0.05, ctrl_edit_rate_thres=0.1,
                                 ignore_control=['barcode', 'TAM'],
                                 on_missing_control='warn') -> pd.DataFrame:
    """Normalize per-read editing or indel rates against editor-free controls."""
    return _normalize(data, col_normalize_over, col_rate, col_control, ignore_control,
                      th_ctrl_missing, ctrl_edit_rate_thres, on_missing_control)


def aggregate_replicates(data: pd.DataFrame, cols_group = ['sublibrary','site_id', 'PAM', 'bioreplicate'], col_rate='edit_rate') -> pd.DataFrame:
    """ Aggregate replicates by taking the mean editing rate. 

    Args:
        data (pd.DataFrame): data to aggregate
        cols_group (List[str]): columns to group by
        col_rate (str): column name for editing rate
    Returns:
        pd.DataFrame: aggregated data
    """    
    # columns outside the grouping keys take the first value within each group
    agg_funcs = {col: 'first' for col in data.columns if col not in cols_group + [col_rate]}
    agg_funcs[col_rate] = 'mean'
    return data.groupby(cols_group).agg(agg_funcs).reset_index()

"""Genotype and mutation calling. Extracted from the Schwank-Lab nerai-genotyping package; only the functions used by the notebooks in this repository are kept."""

import time
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import traceback
import pandas as pd


def compute_aa_dist(fasta_file, output_file, ref_len): 
    """
    Computes distribution of amino acids at each position.

    fasta_file: path to the fasta file. 
                Each sequence is assumed to be aligned to the refernce genome starting from the first position
    output_file: path where to store the distribution.
    ref_len:    length of the reference genome.
    """
    start_time = time.time()
    protein_sequences = []
    num_err = 0
    amino_acid_counter = Counter()

    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) != ref_len:
            num_err += 1
            continue
        protein_seq = str(record.seq.translate())
        protein_sequences.append(protein_seq)

        for i, aa in enumerate(protein_seq, start=1):
            amino_acid_counter[(i, aa)] += 1

   

    amino_acids = list("ACDEFGHIKLMNPQRSTVWY*")
    df = pd.DataFrame(0, index=range(1, ref_len // 3 + 1), columns=amino_acids)

    for (i, aa), count in amino_acid_counter.items():
        df.at[i, aa] = count

    df.to_csv(output_file)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Processed {fasta_file}: time {elapsed_time:.2f}s, wrong length in {num_err}/{len(protein_sequences)}")


def find_mutations(seq, ref_sequence):
    """Find all point-mutations between two sequences. 

    Returns a tuple of point mutations. 
    """
    mutations = []
    ref = ref_sequence
    # ignore parts of the sequence that are longer/shorter than reference.
    seq_len = min(len(seq), len(ref))
    for i in range(seq_len):
        if seq[i] != ref[i]:
            mutation = f"{ref[i]}{i+1}{seq[i]}"
            mutations.append(mutation)
    return tuple(mutations)


def get_genotype_distribution(fasta_file, output_file, ref_sequence, aa_level=True):
    """
    Get prevalence of each genotype among sequences.

    Genotype is a set of mutations.

    params: 
        - ref_sequence: reference sequence, should be either AAs or bases, to match aa_level parameter
        - aa_lavel: if True, mutations are called on amino acid level, otherwise on base level
    """

    def translate_sequence(nucleotide_seq):
        return Seq(nucleotide_seq).translate()
    
    genotype_counts = {}
    total = 0

    try: 
        for record in SeqIO.parse(fasta_file, "fasta"):
            target_seq =  translate_sequence(record.seq) if aa_level else record.seq
            target_seq = str(target_seq)
            genotype_counts[target_seq] = genotype_counts.get(target_seq, 0) + 1
            total += 1
    
        # Calculate fractions
        genotype_distribution = {genotype: count / total for genotype, count in genotype_counts.items()}
    
        # Convert to DataFrame with swapped columns and write to file
        df = pd.DataFrame(list(genotype_distribution.items()), columns=['Genotype', 'Fraction'])
        df = df[['Fraction', 'Genotype']]  # Swap columns
        df = df.sort_values(by='Fraction', ascending=False)
        df['Mutations'] = df.Genotype.map(lambda seq: find_mutations(seq, ref_sequence))  #get mutaions posiions like "'A50T' to the ref, here index start from 1
        df.to_csv(output_file, index=False)
        
        print(f"Processed {fasta_file}: {len(genotype_distribution)} genotypes in {total} reads")
    except Exception as ex:
        print(f'Error processing {fasta_file}')
        traceback.print_exc()


def load_aa_dist_per_sample(aa_dist_file): 
    """
    return amino acid distributions at all postions for each sample. looks like 
            A         D         E     ...(other amino acid)
    pos1    0.2     0.0     ... 0.7
    pos2    0.4     0.7     ... 0.6
    ...
    posN    0.0     1.0     ... 1.0
    
    """    

    muts = pd.read_csv(aa_dist_file, index_col=0)
    mut_fracs = muts.divide(muts.sum(axis=1), axis=0)
    return mut_fracs


def find_top_mutations_per_sample(mut_fracs, ref_seq): 
    """
    return top point mutation for all postions for each sample.  looks like :
            TopMutation	Fraction
    pos1       C         0.0
    pos2       A        0.001
    ...
    posN       R        0.006
    """      
    
    ref_seq_series = pd.Series(list(ref_seq), index=mut_fracs.index)
    
    top_mutations = []
    for position, ref_aa in ref_seq_series.items():
        position_data = mut_fracs.loc[position].drop(ref_aa)
        top_mutation = pd.NA if position_data.isna().all() else position_data.idxmax() 
        fraction = position_data.max()
        top_mutations.append((position, top_mutation, fraction))
        
    top_mutations_df = pd.DataFrame(top_mutations, columns=['Position', 'TopMutation', 'Fraction']).set_index('Position')
    return top_mutations_df


def find_top_mutations(samples2aadist, ref_seq, min_frac_th=0.2): 
    """ Identifies top mutations in a set of samples and their frequencines in each sample. 
    
    Function works as following: 
    1. For each sample, each position it idenfies the most prevalent non-WT mutation
    2. For each position it identifies the most prevalent mutation across all samples
    3. It filters positions for with the most prevalent mutation is above the threshold
    4. For each sample and each position it computes prevalence of top-mutation at that position in the concrete sample

    The resulting dataframe looks like this: 
            sample1 sample2 ... sampleN
    pos1    0.2     0.0     ... 0.7
    pos2    0.4     0.7     ... 0.6
    ...
    posN    0.0     1.0     ... 1.0

    
    parameters:
        samples2aadist: dictionary with keys as samples and values as paths to mutation distribution csv files. 
        ref_seq: reference sequence to know which AA value of each position is WT
        min_frac_th: minimum fraction of the mutation in at least one of the samples to be included in final outout.
    return:
        top_mutations_df - identified top mutations 
        mutation_prevalence_df - prevalence of top mutations in each sample.
    """

    muts_per_sample = {
        key:  load_aa_dist_per_sample(aa_dist) for key, aa_dist in samples2aadist.items()
    }
    top_muts_per_sample = {
        key: find_top_mutations_per_sample(df, ref_seq) for key, df in muts_per_sample.items()
    }
    
    top_muts = []
    #positions = list(next(iter(top_muts_per_sample.values())).index)
    top_muts_per_sample_first=list(top_muts_per_sample.values())[0]
    positions=list(top_muts_per_sample_first.index)
    for position in positions:
    # now for each postion, identifies top mutation across a set of samples
    # and get their frequencines in each sample
        max_fraction = 0
        max_mutation = None
        for df in top_muts_per_sample.values(): # this this specific postion,  identifies top mutation across a set of samples
            # Assuming 'TopMutation' and 'Fraction' are column names
            if position in df.index:
                fraction = df.loc[position, 'Fraction']
                mutation = df.loc[position, 'TopMutation']
                if fraction >= max_fraction:
                    max_fraction = fraction
                    max_mutation = mutation
                    
        # Append the max fraction and mutation for the position
        top_muts.append((position, max_mutation, max_fraction))

    # Convert the list to a DataFrame
    top_muts_df = pd.DataFrame(top_muts, columns=['Position', 'TopMutation', 'Fraction']).set_index('Position')
    top_muts_filtered_df =  top_muts_df[top_muts_df['Fraction'] > min_frac_th] 
    
    #now check the frequeny of this top mutions in all sample 
    top_mut_frac_per_sample = {}
    for pos, (mut, _) in top_muts_filtered_df.iterrows():
        top_mut_frac_per_sample[pos] = [df.loc[pos, mut] for df in muts_per_sample.values()]
        
        
    return top_muts_filtered_df, pd.DataFrame(top_mut_frac_per_sample, index=muts_per_sample.keys()).T


def filter_top_muts_per_sample(aa_dist, ref_aa_seq, min_frac_th=0.2): 

    muts_per_sample = load_aa_dist_per_sample(aa_dist) 

    top_muts_per_sample = find_top_mutations_per_sample(muts_per_sample, ref_aa_seq)
    
    top_muts_filtered_per_sample=top_muts_per_sample[top_muts_per_sample['Fraction'] > min_frac_th]
        
    return top_muts_filtered_per_sample


def find_mutations_per_sample(mut_fracs, ref_seq): 
    """
    return top point mutations for all positions for each sample.  looks like :
            Mutation	Fraction
    pos1       C         0.0
               G         0.002
    pos2       A        0.001
    ...
    posN       R        0.006
    """      
    
    ref_seq_series = pd.Series(list(ref_seq), index=mut_fracs.index)
    
    top_mutations = []
    for position, ref_aa in ref_seq_series.items():
        position_data = mut_fracs.loc[position].drop(ref_aa)
        non_zero_mutations = position_data[position_data > 0]
        for mutation, fraction in non_zero_mutations.items():
            top_mutations.append((position, mutation, fraction))
        
    top_mutations_df = pd.DataFrame(top_mutations, columns=['Position', 'Mutation', 'Fraction']).set_index('Position')
    return top_mutations_df


def filter_muts_per_sample(aa_dist, ref_aa_seq, min_frac_th=0.2): 

    muts_per_sample = load_aa_dist_per_sample(aa_dist) 

    top_muts_per_sample = find_mutations_per_sample(muts_per_sample, ref_aa_seq)
    
    top_muts_filtered_per_sample=top_muts_per_sample[top_muts_per_sample['Fraction'] > min_frac_th]
        
    return top_muts_filtered_per_sample

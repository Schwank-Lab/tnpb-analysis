"""PacBio demultiplexing. Extracted from the Schwank-Lab nerai-genotyping package; only the functions used by the notebooks in this repository are kept."""

import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import gzip


BARCODE_LEN = 6


def demultiplex_samples(fastq_path,multiplexed_id_mapping, meta_info, DEST_FASTQ_DIR, input_is_gz=False,barcode_len=BARCODE_LEN,):
    """
        demultiplex reads from fastq file based on barcodes. 

        Assumes the following structure of the read: 
        fwd_tp_barcode --- read --- reverse_complement(rev_tp_barcode)

        params:
            - barcode_to_sample_dict: maps a combination of (fwd_tp_barcode, rev_tp_barcode) to the sample
        returns: 
            dataframe with read counts for each demultiplexed sample. 
            Contains 'unmatched' sample for the number of reads that could not be attributed to any sample
        
    """
    assert os.path.exists(fastq_path), f'{fastq_path} does not exist'
    barcode_counts = {}
    
    ori_fastq_sample = Path(fastq_path).name.split(".fastq")[0]
    
    if ori_fastq_sample in multiplexed_id_mapping: #mapping multiplex_fastq filename to multiple_sample name (which contain the well information )
        fastq_sample = multiplexed_id_mapping[ori_fastq_sample]
        meta_info=meta_info.loc[meta_info["multiplexed_sample"]==fastq_sample,:]
        
        # Use single_sample_mode if sample_name equals multiplexed_sample for all rows
        if (meta_info["sample_name"] == meta_info["multiplexed_sample"]).all():
            single_sample =fastq_sample
            sample_to_fastq_dict = {single_sample: f'{DEST_FASTQ_DIR}/{single_sample}.fastq'}
            Path(sample_to_fastq_dict[single_sample]).touch(exist_ok=False)
            read_counts = {single_sample: 0}
            single_sample_mode = True
        else:
            single_sample_mode = False
            # Initialize the barcode2tp_dict dictionary
            barcode_to_sample_dict = {}
            # Group by barcode_forward and barcode_reverse, then check if there's exactly one unique timepoint/plate for each group
            for k, v in meta_info.groupby(['barcode_forward', 'barcode_reverse'])['timepoint']:
                # Assert that there is exactly one unique timepoint for each barcode pair
                assert len(v.unique()) == 1, f"Inconsistent timepoint mapping for barcode pair: {k}"
                # If the assertion passes, add to the dictionary
                barcode_to_sample_dict[k] = v.iloc[0]
            
            sample_to_fastq_dict={tp: f'{DEST_FASTQ_DIR}/{fastq_sample}_{tp}.fastq' for tp in barcode_to_sample_dict.values()} # now this step use all barcode but not sample use all barcode so some emtpy files will be created
            #here each sample refers to result of each timepoint from the mixed input sample  
            # so in the end, each  mixed input sample fastq file  are split into 24 time points of fastq ile 
            for file in sample_to_fastq_dict.values():
                Path(file).touch(exist_ok=False) # this step in needed, otherwise, as the same result will be apppended to same file again and again 
            
            # Initialize a dictionary to count reads per timepoint
            read_counts = {sample: 0 for sample in sample_to_fastq_dict.keys()} 
            
        count_unmatched = 0
        #with open(fastq_path, 'r') as f:
        with open_fastq(fastq_path, 'r') as f:  # Open the file using open_fastq
            for record in SeqIO.parse(f, 'fastq'):
                read = str(record.seq)
                # Extract forward and reverse barcodes
                fwd_barcode = read[:barcode_len]  # the forward barcode must be at the very start of the read
                rev_barcode = str(Seq(read[-barcode_len:]).reverse_complement())
                read_stripped = read[barcode_len:-barcode_len]
                barcode_counts[(fwd_barcode, rev_barcode)] = barcode_counts.get((fwd_barcode, rev_barcode), 0) + 1
                
                if single_sample_mode: # 
                    new_record = SeqRecord(Seq(read_stripped), id=record.id, description=record.description)
                    new_record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][barcode_len:-barcode_len]
                    with open(sample_to_fastq_dict[single_sample], 'a') as sample_fastq:
                        SeqIO.write(new_record, sample_fastq, 'fastq')
                    read_counts[single_sample] += 1
                else:
                    # Check if the barcode pair is in the dictionary
                    if (fwd_barcode, rev_barcode) in barcode_to_sample_dict:
                        sample = barcode_to_sample_dict[(fwd_barcode, rev_barcode)]
                        new_record = SeqRecord(Seq(read_stripped), id=record.id, description=record.description)
                        new_record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][barcode_len:-barcode_len]
                        with open(sample_to_fastq_dict[sample], 'a') as sample_fastq:
                            # NOTE: writes `record`, so forward-strand reads keep their barcodes,
                            # while reverse-strand reads below are written stripped. BWA soft-clips
                            # the retained barcode, so downstream reconstruction is unaffected.
                            SeqIO.write(record, sample_fastq, 'fastq')
                        # Increment the read count for the matched timepoint
                        read_counts[sample] += 1
                        
                        
                    elif (rev_barcode, fwd_barcode) in barcode_to_sample_dict:
                        sample = barcode_to_sample_dict[(rev_barcode, fwd_barcode)]
                        # read is on the reverse strand, compute reverse complement
                        read_stripped_rev_comp = str(Seq(read_stripped).reverse_complement())
                        new_record = SeqRecord(Seq(read_stripped_rev_comp), id=record.id, description=record.description)
                        new_record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][barcode_len:-barcode_len][::-1]
                        with open(sample_to_fastq_dict[sample], 'a') as tp_fastq:
                            SeqIO.write(new_record, tp_fastq, 'fastq')
                        read_counts[sample] += 1
                    else:
                        count_unmatched += 1
                    
                    
        read_counts = pd.DataFrame({sample: [read_count, sample_to_fastq_dict[sample]] for sample, read_count in read_counts.items()}, 
                                    index=['Read Count', 'fastq']).T
        fastq = os.path.basename(fastq_path)
        count_matched = read_counts['Read Count'].sum()
        stats = pd.DataFrame({'fastq': [fastq_sample], 'matched': [count_matched], 'unmatched': [count_unmatched], 
                            'matched_frac':[count_matched / (count_matched + count_unmatched)]})

        return read_counts, barcode_counts, stats


def open_fastq(fastq_path, mode='r'):
    if isinstance(fastq_path, Path): 
        fastq_path = str(fastq_path)
    if fastq_path.endswith('.gz'):
        return gzip.open(fastq_path, mode+'t')
    else:
        return open(fastq_path, mode)

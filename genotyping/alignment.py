"""Read alignment helpers. Extracted from the Schwank-Lab nerai-genotyping package; only the functions used by the notebooks in this repository are kept."""

import subprocess
import pysam
import pandas as pd


def run_bwa(
    fasta_file,
    ref_genome_file: str,
    output_file,
    num_threads=4,
    seed_length=19,
    match_score=1,                 # BWA default: 1   → Higher = more reward for matches
    mismatch_penalty=4,            # BWA default: 4   → Higher = mismatches discouraged
    gap_open_penalty=(20, 20),       # BWA default: 6,6 → Higher = fewer indels introduced
    gap_extension_penalty=(1, 1),  # BWA default: 1,1 → Higher = shorter indels preferred
    clipping_penalty=(5, 5)        # BWA default: 5,5 → Higher = less soft/hard clipping
):
    """
    Aligns reads to a reference genome using BWA-MEM.

    PARAMETERS & EFFECTS:
    ----------------------
    seed_length (-k)
        • Default: 19
        • Increase → Faster but less sensitive (misses small matches)
        • Decrease → More sensitive but slower

    match_score (-A)
        • Default: 1
        • Increase → Favors high-identity alignments, improves mapping in clean data

    mismatch_penalty (-B)
        • Default: 4
        • Increase → Penalizes mismatches more, intolerant to SNPs or noisy reads
        • Decrease → More tolerant to sequencing errors or variants

    gap_open_penalty (-O del,ins)
        • Default: (6,6)
        • Increase → Indels become unlikely; alignment avoids starting gaps
        • Decrease → BWA becomes more willing to call insertions/deletions

    gap_extension_penalty (-E del,ins)
        • Default: (1,1)
        • Increase → Indels become shorter (extensions costly)
        • Decrease → Allows long indels (useful for long-read data)

    clipping_penalty (-L 5',3')
        • Default: (5,5)
        • Increase → Avoids clipping; keeps full read
        • Decrease → Allows trimming of noisy ends (useful for adapters or low-quality bases)

    RETURN VALUE:
        Total alignment score = sum of AS tags from SAM file
        (This reflects how good the overall alignment was.)
    """

    # Build the BWA command
    bwa_command = [
        "bwa", "mem",
        "-t", str(num_threads),
        "-k", str(seed_length),
        "-v", str(0),
        "-A", str(match_score),
        "-B", str(mismatch_penalty),
        "-O", f"{gap_open_penalty[0]},{gap_open_penalty[1]}",
        "-E", f"{gap_extension_penalty[0]},{gap_extension_penalty[1]}",
        "-L", f"{clipping_penalty[0]},{clipping_penalty[1]}",
        ref_genome_file,
        str(fasta_file)
    ]

    print("Running BWA command:\n" + " ".join(bwa_command) + "\n")
    result = subprocess.run(
        " ".join(bwa_command) + f" > {output_file}",
        shell=True, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"bwa mem failed with exit code {result.returncode} for {fasta_file}. "
            "Without this check the alignment file is left empty and the failure only "
            "surfaces later as a confusing parse error. Is bwa on PATH? "
            "See the Environment section of README.md.\n"
            + (result.stderr or "").strip()[-600:]
        )
    return result


CIGAR_MATCH = 0


CIGAR_INS = 1


CIGAR_DEL = 2


CIGAR_SOFT_CLIP = 4


def reconstruct_sequences_checkMaping(sam_file, ref_seq, output_file, sample_name=None,return_simple=True,
                                      min_coverage_frac=0.0):
    """
    Reconstruct sequences from alignment file, ignoring indels.

    Deletions are filled with bases from reference genome. Insertions are discarded.

    min_coverage_frac: discard reads whose reconstruction spans less than this
        fraction of the reference. The published results were produced without
        this filter (0.0): short alignments were kept and the uncovered part
        filled with reference bases. Set to 0.8 to discard them instead.
    """
    
    try:
        sample_name = sample_name if sample_name is not None else sam_file
        # reconstructed_sequences = []

        with pysam.AlignmentFile(sam_file, "r", check_sq=False) as samfile:
            num_indels = 0
            total_reads = 0
            mapped_reads = 0 # reads taht are matched to the reference genome 
            full_reads = 0 # reads that are aligned from the beginning of the reference genome.
            num_reads_with_indels = 0
            num_reads_with_inserts=0
            num_reads_with_deletes=0
            all_del_pos_list_list=list()
            all_insert_pos_list_list=list()
            all_reconstructed_seq=dict()
            for read in samfile:
                del_pos_list=[]#[read.query_name]
                insert_pos_list=list()
                has_insert=False
                has_delete=False
                total_reads += 1
                
                if read.is_unmapped:
                    continue
                mapped_reads += 1
                
                reconstructed_seq = []
                
                if read.reference_start > 3:  # alignment must start at the reference start, allowing a few missing bases
                    continue 
                #full_reads +=1 
                
                ref_pos = read.reference_start #0
                if ref_pos>0: # insert missing bases at the begining.
                    reconstructed_seq=[ref_seq[i] for i in range(ref_pos)]
                
                query_pos = 0 # start position of each fragment in the query read 
                previous_s=None
                for op, length in read.cigar:# read is splited into multiple fragment, and each fragment is label by its matching status to the refence 
                    # demultiplexing orients reads to the forward strand, so only that case is handled
                    if op == CIGAR_MATCH:
                        reconstructed_seq.append(str(read.query_sequence[query_pos:query_pos + length]))
                        query_pos += length
                        ref_pos += length
                        previous_s=ref_pos-1
                    elif op == CIGAR_DEL: # Deletions are filled with bases from reference genome.
                        reconstructed_seq.append(str(ref_seq[ref_pos:ref_pos + length]))
                        del_pos_list.extend(list(range(ref_pos,ref_pos+length)))
                        ref_pos += length
                        num_indels += 1
                        has_delete=True
                        previous_s=ref_pos-1
                        
                    elif op == CIGAR_INS: # Insertions are discarded.
                        query_pos += length
                        num_indels += 1
                        #has_indel = True
                        has_insert=True
                        insert_pos_list.extend([previous_s]*length)
                    elif op == CIGAR_SOFT_CLIP: # soft_clips are also discarded, they are also treated like insertation . adapters ?!
                        # In our PACE setup, the essential gene III on the M13 selection phage (SP) is replaced by a catalytically dead CjCas9 (dCjCas9), which is N-terminally fused to the ω subunit of bacterial RNA polymerase. 
                        # https://www.nature.com/articles/s41589-023-01427-x
                        query_pos += length

                all_del_pos_list_list.append(del_pos_list)
                all_insert_pos_list_list.append(insert_pos_list)
                
                # insert missing bases at the end.
                reconstructed_seq = ''.join(reconstructed_seq)
                
                if len(reconstructed_seq) < min_coverage_frac * len(ref_seq):  # too little of the reference covered
                    continue
                #full_reads +=1  # update it to only consider a full read is read itself have large mapping with ref seq 
                
                if len(reconstructed_seq) < len(ref_seq):
                    # print(read)
                    # break
                    assert ref_pos == len(reconstructed_seq), f'{ref_pos} != {len(reconstructed_seq)}'
                    reconstructed_seq += ref_seq[ref_pos:len(ref_seq)]
                
                if read.query_name not in all_reconstructed_seq:
                    all_reconstructed_seq[read.query_name]=str(reconstructed_seq)
                else:
                    current_uppercase_count = sum(1 for char in all_reconstructed_seq[read.query_name] if char.isupper())
                    new_uppercase_count = sum(1 for char in str(reconstructed_seq) if char.isupper())
                    if new_uppercase_count>current_uppercase_count: # in the case of multiple mapping for same query read, keep the one need least replacement from ref_seq. 
                        all_reconstructed_seq[read.query_name]=str(reconstructed_seq)
            
                # reconstructed_sequences.append({
                #     'read_name': read.query_name,
                #     'reconstructed_seq': reconstructed_seq
                # })
                # assert len(reconstructed_seq) == len(ref_seq), f'{len(reconstructed_seq)} != {len(ref_seq)}\n{reconstructed_seq}\n{ref_seq}'
                num_reads_with_inserts += 1 if has_insert else 0
                num_reads_with_deletes += 1 if has_delete else 0
                num_reads_with_indels +=1 if (has_insert) or (has_delete) else 0

        
        full_reads=len(all_reconstructed_seq)
        with open(output_file, "w") as outfile:
            for k,v in all_reconstructed_seq.items():
                outfile.write(f">{k}\n")
                outfile.write(f"{v}\n")
        
        stats = pd.Series({
            'Sample': sample_name,  # capitalised: 'sample' collides with DataFrame.sample
            'total_reads': total_reads,
            'mapped_reads': mapped_reads,
            'full_reads': full_reads,
            'num_reads_with_inserts': num_reads_with_inserts,
            'num_reads_with_deletes': num_reads_with_deletes,
            'num_reads_with_indels': num_reads_with_indels,
            'num_indels': num_indels,
        })
        
        if return_simple:
            return stats
        else:
            return all_reconstructed_seq,all_insert_pos_list_list,all_del_pos_list_list,stats
    except Exception as e:
        print(f"Error processing {sam_file} with reference {ref_seq}. Output file: {output_file}, Sample name: {sample_name}")
        raise


def reconstruct_sequences(sam_file, ref_seq, output_file, sample_name=None, min_coverage_frac=0.0):
    """
    Reconstruct sequences from alignment file, ignoring indels.

    Deletions are filled with bases from reference genome. Insertions are discarded.

    min_coverage_frac: discard reads whose reconstruction spans less than this
        fraction of the reference. The published results were produced without
        this filter (0.0): short alignments were kept and the uncovered part
        filled with reference bases. Set to 0.8 to discard them instead.
    """
    sample_name = sample_name if sample_name is not None else sam_file
    reconstructed_sequences = []

    with pysam.AlignmentFile(sam_file, "r", check_sq=False) as samfile:
        num_indels = 0
        total_reads = 0
        mapped_reads = 0
        full_reads = 0 # reads that are aligned from the beginning of the reference genome.
        num_reads_with_indels = 0
        for read in samfile:
            has_indel = False
            total_reads += 1
            if read.is_unmapped:
                continue
            mapped_reads += 1
            reconstructed_seq = []
            
            if read.reference_start > 3:  # alignment must start at the reference start, allowing a few missing bases
                continue 
            #full_reads +=1 
            
            ref_pos = read.reference_start #0
            if ref_pos>0: # insert missing bases at the begining.
                reconstructed_seq=[ref_seq[i] for i in range(ref_pos)]
                
            # if read.reference_start > 0:
            #     continue 
            # full_reads +=1 
            
            # ref_pos = 0
            query_pos = 0
            for op, length in read.cigar:
                if op == CIGAR_MATCH:
                    reconstructed_seq.append(str(read.query_sequence[query_pos:query_pos + length]))
                    query_pos += length
                    ref_pos += length
                elif op == CIGAR_INS:
                    query_pos += length
                    num_indels += 1
                    has_indel = True
                elif op == CIGAR_DEL:
                    reconstructed_seq.append(str(ref_seq[ref_pos:ref_pos + length]))
                    ref_pos += length
                    num_indels += 1
                    has_indel = True
                elif op == CIGAR_SOFT_CLIP:
                    query_pos += length
               
            # insert missing bases at the end.
            reconstructed_seq = ''.join(reconstructed_seq)
            
            if len(reconstructed_seq) < min_coverage_frac * len(ref_seq):  # too little of the reference covered
                continue
            full_reads +=1  # update it to only consider a full read is read itself have large mapping with ref seq 
            
            if len(reconstructed_seq) < len(ref_seq):
                assert ref_pos == len(reconstructed_seq), f'{ref_pos} != {len(reconstructed_seq)}'
                reconstructed_seq += ref_seq[ref_pos:len(ref_seq)]
            
            reconstructed_sequences.append({
                'read_name': read.query_name,
                'reconstructed_seq': reconstructed_seq
            })
            # assert len(reconstructed_seq) == len(ref_seq), f'{len(reconstructed_seq)} != {len(ref_seq)}\n{reconstructed_seq}\n{ref_seq}'
            num_reads_with_indels += 1 if has_indel else 0

    with open(output_file, "w") as outfile:
        for entry in reconstructed_sequences:
            outfile.write(f">{entry['read_name']}\n")
            outfile.write(f"{entry['reconstructed_seq']}\n")
    
    stats = pd.Series({
        'sample': sample_name,
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'full_reads': full_reads,
        'num_reads_with_indels': num_reads_with_indels,
        'num_indels': num_indels,
    })
    return stats

# This script is an alternate way of running the pipeline
# This pipeline takes in a single metagenomic reads fastq file as input and runs the entire pipeline on that sample
# Also handles rerunning. To rerun samples it will pick up where it left off, see bottom of code for this


# --------------------------------------------------- #
                # Start of pipeline
# --------------------------------------------------- #
source /mfs/bens/miniconda3/etc/profile.d/conda.sh
sample=$1

# Check if there is a final pipeline.done file , if so, then exit 0. This sample is already complete
if [[ -f "/mfs/bens/rdt_plus_data/xx_pipeline_tracking/${sample}.pipeline.done" ]]; then
    echo "The pipeline is already complete for sample: $sample. Skipping $sample."
    exit 0
fi
echo -e "\033[0;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Starting assembly pipeline for sample: ${sample}\033[0m"
uname -a


# --------------------------------------------------- #
# Chopper: read cleaning and filtering by Q-score 20
# --------------------------------------------------- #
conda activate chopper_env
reads_file=/mfs/bens/rdt_plus_data/00_reads/${sample}.fastq.gz
output_path=/mfs/bens/rdt_plus_data/01_chopper

# Skip step if already done (by checking if the done file exists)
if [[ -f "$output_path/${sample}.chopper.done" && -s "$output_path/${sample}_cleaned.fastq.gz" ]]; then
    echo "Chopper - Skipping $sample — already completed."
else
    rm -f "$output_path/${sample}_cleaned.fastq.gz"
    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running Chopper for sample: ${sample}\033[0m"
    chopper --trim-approach trim-by-quality --cutoff 20 --minlength 500 -i $reads_file | gzip > $output_path/${sample}_cleaned.fastq.gz && echo "done" > "$output_path/${sample}.chopper.done"  
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished Chopper for sample: ${sample}\033[0m"
fi
conda deactivate
echo


# --------------------------------------------------- #
# Sra-human-scrubber: remove human reads from read file
# --------------------------------------------------- #
reads_file=/mfs/bens/rdt_plus_data/01_chopper/${sample}_cleaned.fastq.gz
output_path=/mfs/bens/rdt_plus_data/02_sra-human-scrubber

if [[ -f "$output_path/${sample}.sra-human-scrubber.done" && -s "$output_path/${sample}_decontaminated.fastq.gz" ]]; then
    echo "Sra-human-scrubber - Skipping $sample — already completed."
else
    # Check that the done file for the previous step exists first
    if [[ ! -f "/mfs/bens/rdt_plus_data/01_chopper/${sample}.chopper.done" ]]; then
        echo -e "\033[1;31mERROR: Sra-human-scrubber cannot run because Chopper failed for $sample\033[0m"
        echo "Expected file not found: /mfs/bens/rdt_plus_data/01_chopper/${sample}.chopper.done"
        exit 1
    fi
    rm -f "$output_path/${sample}_decontaminated.fastq.gz"

    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running sra-human-scrubber for sample: ${sample}\033[0m"
    gunzip -c $reads_file | /mfs/bens/software/sra-human-scrubber/scripts/scrub.sh -x -p 4 | gzip > ${output_path}/${sample}_decontaminated.fastq.gz && echo "done" > "$output_path/${sample}.sra-human-scrubber.done"
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished sra-human-scrubber for sample: ${sample}\033[0m"
fi
echo


# --------------------------------------------------- #
            # Flye: Metagenome assembly
# --------------------------------------------------- #
conda activate flye_env
reads_file=/mfs/bens/rdt_plus_data/02_sra-human-scrubber/${sample}_decontaminated.fastq.gz
output_path=/mfs/bens/rdt_plus_data/03_flye/$sample

if [[ -f "$output_path/${sample}.flye.done" && -s "$output_path/assembly.fasta" ]]; then
    echo "Flye - Skipping $sample — already completed."
else
    if [[ ! -f "/mfs/bens/rdt_plus_data/02_sra-human-scrubber/${sample}.sra-human-scrubber.done" ]]; then
        echo -e "\033[1;31mERROR: Flye cannot run because sra-human-scrubber failed for $sample\033[0m"
        echo "Expected file not found: /mfs/bens/rdt_plus_data/02_sra-human-scrubber/${sample}.sra-human-scrubber.done"
        exit 1
    fi
    rm -rf "$output_path"

    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running Flye v2.9.6 for sample: ${sample}\033[0m"
    mkdir -p $output_path
    flye --nano-hq $reads_file --out-dir $output_path --meta -t 4 && echo 'done' > $output_path/${sample}.flye.done
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished Flye for sample: ${sample}\033[0m"
fi
conda deactivate
echo


# --------------------------------------------------- #
        # Filter out contigs <1.5kb with SeqKit
# --------------------------------------------------- #
conda activate seqkit_env
metagenome_file=/mfs/bens/rdt_plus_data/03_flye/$sample/assembly.fasta
output_path=/mfs/bens/rdt_plus_data/03_flye/$sample # We will store our filtered contigs in the same directory, but under a new name

if [[ -f "$output_path/${sample}.seqkit.done" && -s "$output_path/${sample}_filtered_1.5kb.fna" ]]; then
    echo "Seqkit Filter - Skipping $sample — already completed."
else
    if [[ ! -f "/mfs/bens/rdt_plus_data/03_flye/$sample/${sample}.flye.done" ]]; then
        echo -e "\033[1;31mERROR: Seqkit Filter cannot run because Flye failed for $sample\033[0m"
        echo "Expected file not found: /mfs/bens/rdt_plus_data/03_flye/$sample/${sample}.flye.done"
        exit 1
    fi
    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running Seqkit Filter for sample: ${sample}\033[0m"
    seqkit seq -m 1500 $metagenome_file > ${output_path}/${sample}_filtered_1.5kb.fna && echo 'done' > $output_path/${sample}.seqkit.done
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished Seqkit Filter for sample: ${sample}\033[0m"
fi
conda deactivate
echo


# --------------------------------------------------- #
        # GeNomad: Identify viral contigs
# --------------------------------------------------- #
conda activate genomad_env
metagenome_file=/mfs/bens/rdt_plus_data/03_flye/$sample/${sample}_filtered_1.5kb.fna
output_path=/mfs/bens/rdt_plus_data/04_genomad/$sample
genomad_database=/mfs/databases/genomad_db_v1.9

if [[ -f "$output_path/${sample}.genomad.done" && -s "$output_path/${sample}_filtered_1.5kb_find_proviruses/${sample}_filtered_1.5kb_provirus.fna" ]]; then
    echo "GeNomad - Skipping $sample — already completed."
else
    if [[ ! -f "/mfs/bens/rdt_plus_data/03_flye/$sample/${sample}.seqkit.done" ]]; then
        echo -e "\033[1;31mERROR: GeNomad cannot run because Seqkit Filter failed for $sample\033[0m"
        echo "Expected file not found: /mfs/bens/rdt_plus_data/03_flye/$sample/${sample}.seqkit.done"
        exit 1
    fi
    rm -rf "$output_path"

    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running GeNomad for sample: ${sample}\033[0m"
    genomad end-to-end --cleanup --splits 8 "$metagenome_file" "$output_path" "$genomad_database" && echo 'done' > $output_path/${sample}.genomad.done
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished GeNomad for sample: ${sample}\033[0m"
fi
conda deactivate
echo


# --------------------------------------------------- #
# CheckV: Assess quality and completeness of viral contigs
# --------------------------------------------------- #
conda activate checkv_env
input_phages=/mfs/bens/rdt_plus_data/04_genomad/$sample/${sample}_filtered_1.5kb_find_proviruses/${sample}_filtered_1.5kb_provirus.fna # Only prophages
output_path=/mfs/bens/rdt_plus_data/05_checkv/$sample
export CHECKVDB=/mfs/databases/checkv-db-v1.5   # CheckV requires database directory to be exported

# TO DO 
    # Revisit my && statements, they dont work as intended. CheckV produes and internal error but still return exit code 0. 
    # Check to see if the quality_summary.tsv file is made AND the done file exists. If not, then delete the directory and restart
if [[ -f "$output_path/${sample}.checkv.done" && -s "$output_path/quality_summary.tsv" ]]; then
    echo "CheckV - Skipping $sample — already completed."
else
    if [[ ! -f "/mfs/bens/rdt_plus_data/04_genomad/$sample/${sample}.genomad.done" ]]; then
        echo -e "\033[1;31mERROR: CheckV cannot run because GeNomad failed for $sample\033[0m"

        echo "Expected file not found: /mfs/bens/rdt_plus_data/04_genomad/$sample/${sample}.genomad.done"
        exit 1
    fi
    # Remove partial results directory if it didn't complete fully
    rm -rf "$output_path"

    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running CheckV for sample: ${sample}\033[0m"
    checkv end_to_end $input_phages $output_path -t 4

    # Count data rows (excluding header)
    n_prophages=$(awk 'NR>1 {c++} END {print c+0}' "$output_path/quality_summary.tsv")
    if [[ "$n_prophages" -eq 0 ]]; then
        echo "CheckV completed: no prophages detected for $sample"
        # Mark CheckV as done
        echo 'done' > "$output_path/${sample}.checkv.done"
        # Optional marker for bookkeeping
        touch "$output_path/${sample}.no_prophages"
        # Mark full pipeline as successfully completed
        echo 'done' > /mfs/bens/rdt_plus_data/xx_pipeline_tracking/${sample}.pipeline.done
        conda deactivate
        exit 0
    fi

    # Normal successful case (≥1 prophage)
    echo 'done' > "$output_path/${sample}.checkv.done"
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished CheckV for sample: ${sample}\033[0m"
fi

conda deactivate
echo


# --------------------------------------------------- #
        # Filter out Not-determined phage contigs
# --------------------------------------------------- #
conda activate seqkit_env
phage_fasta=/mfs/bens/rdt_plus_data/04_genomad/$sample/${sample}_filtered_1.5kb_find_proviruses/${sample}_filtered_1.5kb_provirus.fna
filtered_phage_names=/mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}_phages_filtered.txt
filtered_phages=/mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}_filtered_phages.fna
input_file=/mfs/bens/rdt_plus_data/05_checkv/$sample/quality_summary.tsv

if [[ -f "$output_path/${sample}.na-filter.done" && -s "$filtered_phages" ]]; then
    echo "na-filter - Skipping $sample — already completed."
else
    if [[ ! -f "/mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}.checkv.done" ]]; then
        echo -e "\033[1;31mERROR: na-filter cannot run because CheckV failed for $sample\033[0m"
        echo "Expected file not found: /mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}.checkv.done"
        exit 1
    fi
    # First check that the input quality_summary.tsv file exists
    if [[ -s "$input_file" ]]; then
        echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Filtering out Not-determined quality phages for sample: ${sample}\033[0m"
        echo "Filtering out phages for: $sample"

        # Create a .txt file with all the names of phage contigs, excluding "Not-determined"
        awk -F '\t' 'NR>1 && $8 != "Not-determined" {print $1}' "$input_file" > "$filtered_phage_names"


        # Now filter all_phages.fna using the list of passed phage contig names
        seqkit grep -f "$filtered_phage_names" "$phage_fasta" > "$filtered_phages" && echo 'done' > $output_path/${sample}.na-filter.done
        echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished na-filter for sample: ${sample}\033[0m"
    else
        # Missing checkv output error handling
        echo "quality_summary.tsv not found for sample $sample"
    fi
fi
conda deactivate
echo


# --------------------------------------------------- #
        # iPHoP: Phage-Host prediction
# --------------------------------------------------- #
conda activate iphop_env
prophages=/mfs/bens/rdt_plus_data/05_checkv/${sample}/${sample}_filtered_phages.fna
output_path=/mfs/bens/rdt_plus_data/06_iphop/$sample
iphop_database=/mfs/databases/iphop_db/Jun_2025_pub_rw

if [[ -f "$output_path/${sample}.iphop.done" && -s "$output_path/Host_prediction_to_genome_m90.csv" ]]; then
    echo "iPHoP - Skipping $sample — already completed."
else
    if [[ ! -f "/mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}.na-filter.done" ]]; then
        echo -e "\033[1;31mERROR: iPHoP cannot run because na-filter failed for $sample\033[0m"
        echo "Expected file not found: /mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}.na-filter.done"
        exit 1
    fi
    rm -rf "$output_path"

    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running iPHoP for sample: ${sample}\033[0m"
    iphop predict --fa_file "$prophages" --db_dir "$iphop_database" --out_dir "$output_path"  && echo 'done' > $output_path/${sample}.iphop.done
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished iPHoP for sample: ${sample}\033[0m"
fi
conda deactivate
echo


# --------------------------------------------------- #
                   # END OF PIPELINE
# --------------------------------------------------- #
echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] |------------------------------------------------------------------------------------------|\033[0m"
echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`]      Checking full pipeline for successful completion of each step for sample: ${sample}\033[0m"
echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] |------------------------------------------------------------------------------------------|\033[0m"

# Each of these .done files should exist for the pipeline to have been successfully completed
required_done_files=(
    "/mfs/bens/rdt_plus_data/01_chopper/${sample}.chopper.done"
    "/mfs/bens/rdt_plus_data/02_sra-human-scrubber/${sample}.sra-human-scrubber.done"
    "/mfs/bens/rdt_plus_data/03_flye/$sample/${sample}.flye.done"
    "/mfs/bens/rdt_plus_data/03_flye/$sample/${sample}.seqkit.done"
    "/mfs/bens/rdt_plus_data/04_genomad/$sample/${sample}.genomad.done"
    "/mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}.checkv.done"
    "/mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}.na-filter.done"
    "/mfs/bens/rdt_plus_data/06_iphop/$sample/${sample}.iphop.done"
)

# Will append to this list any missing .done files
missing_files=()

# Check that each file exists
for file in "${required_done_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        missing_files+=("$file")
    fi
done

# If any missing: print error and exit non-zero
if [[ ${#missing_files[@]} -gt 0 ]]; then
    echo -e "\033[1;31mERROR: Pipeline did NOT finish successfully for $sample.\033[0m"
    echo "Missing .done files:"
    for mf in "${missing_files[@]}"; do
        echo "  - $mf"
    done
    exit 1
# Otherwise, the pipeline has successfully ran. Create a final pipeline.done file for easy tracking
else
    echo -e "\033[1;32m[`date +"%Y-%m-%d %H:%M:%S"`] FULL PIPELINE SUCCESSFULLY COMPLETED for sample: ${sample}\033[0m"
    echo 'done' > /mfs/bens/rdt_plus_data/xx_pipeline_tracking/${sample}.pipeline.done
fi

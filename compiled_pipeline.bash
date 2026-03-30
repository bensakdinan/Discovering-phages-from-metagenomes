# This script is an alternate way of running the pipeline
# This pipeline takes in a single metagenomic reads fastq file as input and runs the entire pipeline on that sample

# --------------------------------------------------- #
                      # Minimap2
# --------------------------------------------------- #
sample=$1

echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Starting P:H ratio pipeline for sample: ${sample}\033[0m"
if [[ -f "/mfs/bens/rdt_plus_data/05_checkv/$sample/${sample}.no_prophages" ]]; then
   echo "Skipping $sample no prophages"
   exit 0
fi

conda activate ratio_env
reads_file=/mfs/bens/rdt_plus_data/00_reads/${sample}.fastq.gz
output_path=/mfs/bens/rdt_plus_data/07_minimap/$sample
assembly_file=/mfs/bens/rdt_plus_data/03_flye/$sample/assembly.fasta

mkdir -p $output_path

if [[ -f "$output_path/${sample}.minimap2.done" && -s "$output_path/${sample}.bam" ]]; then
    echo "Minimap2 - Skipping $sample — already completed."
else
    rm $output_path/${sample}.bam
    rm $output_path/${sample}.bam.bai

    echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Running Minimap2 for sample: ${sample}\033[0m"
    minimap2 -ax map-ont $assembly_file $reads_file | samtools sort -o $output_path/${sample}.bam
    samtools index $output_path/${sample}.bam && echo "done" > "$output_path/${sample}.minimap2.done"
    echo -e "\033[1;32m[`date +\"%Y-%m-%d %H:%M:%S\"`] DONE - Finished Minimap2 for sample: ${sample}\033[0m"
fi

# --------------------------------------------------- #
                # Make prophage BED files 
# --------------------------------------------------- #

genomad_tsv=/mfs/bens/rdt_plus_data/04_genomad/${sample}/${sample}_filtered_1.5kb_find_proviruses/${sample}_filtered_1.5kb_provirus.tsv
qc_list=/mfs/bens/rdt_plus_data/05_checkv/${sample}/${sample}_phages_filtered.txt
assembly=/mfs/bens/rdt_plus_data/03_flye/${sample}/assembly.fasta
output_path=/mfs/bens/rdt_plus_data/08_bed_files/${sample}

mkdir -p "${output_path}"

awk '
  NR==FNR { keep[$1]; next }
  ($1 in keep)
' "${qc_list}" "${genomad_tsv}" \
| awk 'BEGIN{OFS="\t"} {print $2, $3-1, $4, $1}' \
| sort -k1,1 -k2,2n \
> "${output_path}/${sample}_prophages.bed"

# --------------------------------------------------- #
        # Get host complement BED coordinates 
# --------------------------------------------------- #

prophage_bed_dir=/mfs/bens/rdt_plus_data/08_bed_files/$sample
assembly_file=/mfs/bens/rdt_plus_data/03_flye/$sample/assembly.fasta
output_path=/mfs/bens/rdt_plus_data/08_bed_files/$sample

# Make .fai genome file and sort prophages bed 
samtools faidx $assembly_file
cut -f1,2 $assembly_file.fai > $assembly_file.genome
sort -k1,1 -k2,2n $prophage_bed_dir/${sample}_prophages.bed \
    > $prophage_bed_dir/${sample}_prophages.sorted.bed

# Get the complement regions of genomes (host regions)
bedtools complement \
  -i $prophage_bed_dir/${sample}_prophages.sorted.bed \
  -g $assembly_file.genome \
> $output_path/${sample}_host_regions.bed

# --------------------------------------------------- #
                # Calculate P:H Ratios 
# --------------------------------------------------- #
minimap_path=/mfs/bens/rdt_plus_data/07_minimap/$sample
output_path=/mfs/bens/rdt_plus_data/09_ratios/$sample
prophage_bed_dir=/mfs/bens/rdt_plus_data/08_bed_files/$sample

if [[ -f "$output_path/${sample}.ratios.done" && -s "$output_path/${sample}_prophage_host_ratios.tsv" ]]; then
    echo "Ratio calculations - Skipping $sample — already completed."
else
  rm -rf "${output_path}"
  echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Calculating prophage-host ratios for sample: ${sample}\033[0m"

  mkdir -p "${output_path}"

  # Prophage depth
  bedtools coverage \
    -a $prophage_bed_dir/${sample}_prophages.sorted.bed \
    -b $minimap_path/${sample}.bam \
    -mean \
  > $output_path/${sample}_prophage_depth.tsv

  # Host depth
  bedtools coverage \
    -a $prophage_bed_dir/${sample}_host_regions.bed \
    -b $minimap_path/${sample}.bam \
    -mean \
  > $output_path/${sample}_host_depth.tsv

  # One value per contig
  awk 'BEGIN{OFS="\t"}
  {
    sum[$1] += $4
    n[$1]++
  }
  END{
    for (c in sum)
      print c, sum[c]/n[c]
  }' $output_path/${sample}_host_depth.tsv \
  > $output_path/${sample}_host_mean_per_contig.tsv

  # Per prophage
  awk -v sample="$sample" '
  BEGIN{
    FS=OFS="\t"
    print "sample_id","prophage_name","host_contig","prophage_host_ratio"
  }
  NR==FNR {
    host_cov[$1]=$2
    next
  }
  {
    contig=$1
    prophage=$4
    prophage_cov=$5

    if (contig in host_cov && host_cov[contig] > 0)
      ratio = prophage_cov / host_cov[contig]
    else
      ratio = "NA"

    print sample, prophage, contig, ratio
  }' \
  $output_path/${sample}_host_mean_per_contig.tsv \
  $output_path/${sample}_prophage_depth.tsv \
  > $output_path/${sample}_prophage_host_ratios.tsv \
  && \
  [ -s $output_path/${sample}_prophage_host_ratios.tsv ] \
  && \
  echo "done" > "$output_path/${sample}.ratios.done"

  echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Finished calculating prophage-host ratios for sample: ${sample}\033[0m"
fi 

echo -e "\033[1;34m[`date +\"%Y-%m-%d %H:%M:%S\"`] Finished P:H ratio pipeline for sample: ${sample}\033[0m"
conda deactivate

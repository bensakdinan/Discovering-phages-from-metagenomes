# Identifying phages from metagenomes pipeline

## Getting started 
You will need the following software (downloading within conda environments is recommended):
- **Seqkit** for various fasta handling purposes (https://anaconda.org/bioconda/seqkit)
- **Chopper** for trimming and quality controlling (QC) Nanopore reads (https://github.com/wdecoster/chopper) 
- **Flye** for metagenome assembly ([https://anaconda.org/bioconda/flye](https://github.com/mikolmogorov/Flye)
- **GeNomad** for identification of mobile genetic elements, phages, prophages (https://github.com/apcamargo/genomad)
- **CheckV** for assessing viral contig completeness and contiguity (https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-running-checkv)
- **iPHoP** for estimating phage-host pairings (https://bitbucket.org/srouxjgi/iphop/src)

Conda is used throughout this pipeline to manage software and packages. Visit the links pasted above for download instruction for each tool. `conda activate [environment_name]` is used at the top of each script to activate the conda environments for each tool. 

Additionally, you will need to run these commands in the background as they can take multiple hours/days. Running these scripts directly in the shell will most likely result in them getting killed as soon as your computer goes to sleep. You can do this using `nohup ./script.bash > script.01.out &`. This will run your script in the background and will direct all standard out to a log file named "script.01.out". See [scripts](https://github.com/bensakdinan/Discovering-phages-from-metagenomes/tree/b6ab6b8b9b18fe3db6c5176e9afa8fd9afe10ba9/scripts) for `nohup` commands to run each file.

**NOTICE**: If you are working on a cluster that does not support `conda` (ie. SLURM job submission), your scripts will look a little different. To parallelize jobs, you would likely not be using GNU parallel as I am using here.

## Trimming and cleaning reads 
The first step is to clean our reads with **chopper**. This involves removing sequencing adapters, filtering out low-score and short reads.
`--trim-approach trim-by-quality` trims low quality bases from ends of reads using quality score `cutoff` of 20. You should consult the Q-score distribution of your sequencing run to ensure not too much data is lost, but Q20 is a good standard. 
`--minlength 500` filter out reads less than 500bp
Adapter sequences are likely already removed during basecalling with Dorado, but if not, consult `--headcrop <HEADCROP>`

```bash
#!/bin/bash
conda activate chopper_env

# Iterates through a .txt file of your sample names. Adjust n of parallel jobs accordingly
cat /path/to/sample_names.txt | parallel -j 10 '
  # This variable stores the current sample name
  sample={} 

  # Directories
  reads_file=/path/to/reads/${sample}.fastq.gz
  output_path=/path/to/output/00_chopper

  # && will create an out file if the previous command ran and exited successfully
  chopper --trim-approach trim-by-quality --cutoff 20 --minlength 500 -i $reads_file | gzip > $output_path/${sample}_cleaned.fastq.gz && echo "done" > "${sample}.chopper.done"  
'
```

## Metagenome assembly
Now that we have cleaned and trimmed fastq.gz reads, we can assemble our metagenome using **Flye**. There are two main methods for metagenome assembly: reference based, and *de novo*. Reference based approaches are more straightforward, it involves mapping reads to known genomes. Consdering that vast majority of phages do not have reference genomes, we instead take a *de novo* approach. **Flye** accomplishes this by  overlapping reads to construct what is known as an **assembly graph**. 

`--meta` flag for metagenome assembly
`--nano-hq` as our input since we previously cleaned reads with minmum Q20
`-t 8` 8 threads, but adjust accordingly

```bash
#!/bin/bash
conda activate flye_env

cat /path/to/sample_names.txt | parallel -j 10 '
  sample={}

  # Directories
  reads_file=/path/to/00_chopper/${sample}_cleaned.fastq.gz
  output_path=/path/to/output/01_flye/$sample

  # Make output directory for each sample (-p if it does not already exist)
  mkdir -p $output_path

  flye --nano-hq $input_file --out-dir $output_path --meta -t 8 && echo 'done' > $output_path/${sample}.flye.done
'
```

## Filter out contigs <1.5kb
Short contigs can skew taxonomic assignments, add noise that obscures abundance estimation, potentially represent lower quality contigs, or have missing genes. We can filter out contigs below 1.5kb as a decent standard of practice.

```bash
#!/bin/bash
conda activate seqkit_env

cat /path/to/sample_names.txt | parallel -j 10 '
  sample={}

  # Directories
  metagenome_file=/path/to/01_flye/$sample/assembly.fasta
  output_path=/path/to/01_flye/$sample # We will store our filtered contigs in the same directory, but under a new name

  seqkit seq -m 1500 $metagenome_file > ${output_path}/${sample}_filtered_1.5kb.fna
'
```

## Identify viral contigs
Using our assembled and 1.5kb filtered metagenome, we can use **geNomad** to identify mobile genetic elements. This includes, plasmids, viruses + phages, and prophages. GeNomad functions by using a database of over 200,000 viral proteins to annotate and taxonomically classify contigs. Additionally, it can identfy prophages by identifying flanking bacterial DNA and common prophage genes (ie. Integrase). 

`end-to-end` will run the entire geNomad pipeline
`--splits 8` split resources 
`--cleanup` delete intermediate files

```bash
#!/bin/bash
conda activate genomad_env

cat /path/to/sample_names.txt | parallel -j 10 '
  sample={}

  metagenome_file=/path/to/01_flye/$sample/${sample}_cleaned.fastq.gz
  output_path=path/to/output/02_genomad/$sample
  genomad_database=/path/to/genomad_db_v

  genomad end-to-end --cleanup --splits 8 "$metagenome_file" "$output_path" "$genomad_database" && echo 'done' > $output_path/${sample}.genomad.done
'
```

## Assess quality and completeness of viral contigs
Now that we have identified our viral contigs, we can assess their quality and completeness to filter out ones with low scores. First step though, is to decide what phages we want to work with. If you only want to look at extracellular phages, those will be found in `/path/to/genomad_output/$sample/final.contigs_summary/viruses.fna`. If you also want to look at prophages, those are found in `/path/to/genomad_output/$sample/final.contigs_find_proviruses/$sample/final.contigs_provirus.fna`. If you only want to look at specific taxa of phages, you can filter by taxonomy using this .tsv file `/path/to/genomad_output/$sample/final.contigs_summary/final.contigs_virus_summary.tsv`. 

If you want to combine both extracellular phages and prophages, you can concatenate both .fna files like so (iterate for all samples if you have many samples): 
```bash
cat /path/to/genomad_output/$sample/final.contigs_summary/viruses.fna /path/to/genomad_output/$sample/final.contigs_find_proviruses/$sample/final.contigs_provirus.fna > /path/to/genomad_output/$sample/all_phages.fna 
```

Now that you have decided which phages you want to work with (I'll be using all phages for the rest of this tutorial), we have to assess the quality and completeness of their genomes. To do so, we will use CheckV. CheckV also functions by comparing the contiguity and completeness of putative viral genomes to a database of fully assembled genomes. 

`end_to_end` to execute the entire checkV pipeline
`-t 16` 16 threads for computational efficiency, adjusted as needed

```bash
#!/bin/bash
conda activate checkv_env

cat /path/to/sample_names.txt | parallel -j 10 '
  sample={}

  input_phages=/path/to/genomad_output/$sample/all_phages.fna
  output_path=/path/to/output/03_checkv/$sample
  export CHECKVDB=/path/to/checkv-db-v1.5   # CheckV requires database directory to be exported

  checkv end_to_end $input_phages $output_path -t 16 && echo 'done' > $output_path/${sample}.checkv.done
  # You can also add `-d [database_path]` in the command if the export statement does not work
'
```

Now we have calculated the quality and completeness of all of our viral contigs. These stats can be found in the output file `quality_summary.tsv`. Here is an example output:
```
contig_id	contig_length	provirus	proviral_length	gene_count	viral_genes	host_genes	checkv_quality	miuvig_quality	completeness	completeness_method	contamination	kmer_freq	warnings
k141_71124|provirus_4909_72203  67295   No      NA      94      25      3       High-quality    High-quality    100.0   AAI-based (high-confidence)     0.0     1.0 
k141_97103|provirus_92502_148357	55856	No	NA	68	9	6	Medium-quality	Genome-fragment	70.32	HMM-based (lower-bound)	0.0	1.0	
k141_111526|provirus_760_49311	48552	Yes	7888	46	3	20	Low-quality	Genome-fragment	16.04	HMM-based (lower-bound)	83.75	1.0	
k141_1791|provirus_1_9293       9293    No      NA      13      0       1       Not-determined  Genome-fragment NA      NA      0.0     1.0     no viral genes detected

```

We can filter out phage contigs that are`"Not-determined` or `Low-quality`. The script below :
```bash
#!/bin/bash
conda activate seqkit_env

for input_file in /path/to/03_checkv/*/quality_summary.tsv; do
    # Gets sample name from the directory name
    sample=$(basename "$(dirname "$input_file")")

    # This is the concatenated fasta of all our phages 
    phage_fasta=/path/to/genomad_output/${sample}/all_phages.fna
    filtered_phage_names=/path/to/03_checkv/${sample}/${sample}_phages_filtered.txt
    filtered_phages=/path/to/03_checkv/${sample}/${sample}_filtered_phages.fna

    # First check that the input quality_summary.tsv file exists
    if [[ -f "$input_file" ]]; then
        echo "Filtering out phages for: $sample"

        # Create a .txt file with all the names of phage contigs, excluding those with Not-determined of Low-quality 
        awk -F '\t' 'NR>1 && $8!="Not-determined" && $8!="Low-quality" {print $1}' \
            "$input_file" > "$filtered_phage_names"

        # Now filter all_phages.fna using the list of passed phage contig names
        seqkit grep -f "$filtered_phage_names" "$phage_fasta" > "$filtered_phages"
    else
        # Missing checkv output error handling
        echo "quality_summary.tsv not found for sample $sample"
    fi
done
```
After this step is complete, we now have successfully generated a multifasta (.fna) file containing all of our phage contigs `/path/to/genomad_output/${sample}/${sample}_filtered_phages.fna`. From here you can do more analyses, this is completely dependent on your project now! One useful thing to do next is estimate the bacterial hosts for each of the phages you have now identified. 

## Estimate phage hosts
For this we will use iPHoP. This tool works by running 6 different host prediction methods:
- Host-based alignment methods: BLAST, CRISPR, SpacePharer
- Host-based alignment-independent: Propkaryote Host Predictor (PHP), VirHostMatcher
- Host-independent: Random Forest Assignment of Hosts (RaFAH)
iPHoP then aggregates the results of these 6 methods to compute a weighted score for each phage host, then declaring the most confident prediction.

```bash
#!/bin/bash
conda activate iphop_env

cat /mfs/bens/aim2/yaffe_data/day1_srs_ids.txt | parallel -j 10 '
  sample={}

  prophages=/path/to/03_checkv/${sample}/${sample}_filtered_phages.fna
  iphop_database=/path/to/iphop/database
  output_dir=/path/to/output/04_iphop/$sample
  
  iphop predict --fa_file "$prophages" --db_dir "$iphop_database" --out_dir "$output_dir"  && echo 'done' > /path/to/04_iphop/${sample}.iphop.done
'
```

The main output from iPHoP will be `Host_prediction_to_genome_m90.csv`, which looks like this:
```
Virus,Host genome,Host taxonomy,Main method,Confidence score,Additional methods
k141_100819|provirus_20279_27577,RS_GCF_025150245.1,d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Dorea;s__Dorea formicigenerans,blast,91.20,iPHoP-RF;85.40
k141_100819|provirus_20279_27577,RS_GCF_025149785.1,d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Bariatricus;s__Bariatricus comes,blast,90.30,iPHoP-RF;75.20
k141_100819|provirus_20279_27577,RS_GCF_014287705.1,d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Dorea_D;s__Dorea_D hominis,blast,90.20,iPHoP-RF;80.80
k141_105171|provirus_22169_37874,GB_GCA_963588145.1,d__Bacteria;p__Bacillota;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Gemmiger;s__Gemmiger sp963588145,blast,95.00,iPHoP-RF;94.10
k141_105171|provirus_22169_37874,GB_GCA_025758045.1,d__Bacteria;p__Bacillota;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Gemmiger;s__Gemmiger sp900539695,iPHoP-RF,94.10,blast;92.20
```

## Further directions
There are many analyses to run, I recommend also clustering these phages into vOTUs using vcontact2 (https://bitbucket.org/MAVERICLab/vcontact2)

## Thank you for taking a look through this pipeline, good luck and happy phage-hunting!

# Identifying-phages-from-metagenomes

**"What about the phages?"** is a question I often see raised by microbiologists during the question period following a seminar. It's easy to neglect the phage fraction in a given biological sample, and it's all too common to forget about these microbial parasites. Phages are uniquely difficult to study due to their lack of a single marker gene. This greatly complicates sequencing since there is no universal PCR target, nor is there an easy genetic identification for phages. The result is that phages are often ignored or forgotten about. 

**I don't want phages to be ignored.** 

To overcome this obstacle, metagenomic sequencing is used to capture phage DNA from samples. A metagenome is the collective genetic material of any and all organisms within a given sample (bacteria, viruses, fungi, eukaryotes, etc). Put simply, by sequencing *everything* you *will* capture phage DNA. Metagenomic sequencing generates immense volumes of data which certainly provides immense potential, but this high-throughput output also significantly complicates analysis. I beilieve this to be a major barrier of entry to metagenomics.

Metagenomics can be daunting for microbiologists to start exploring, especially for those with little to no experience in bioinformatics. The immense volume of data in metagenomics can be intimidating and I often hear from fellow microbiologists that they don't exactly know where to begin. 

However, identifying phages from metagenomes does not need to be complicated. Here I want to provide **[this basic pipeline](https://github.com/bensakdinan/Identifying-phages-from-metagenomes/blob/57223ad5c71596f36c5ff666f38b2118843c35de/pipeline.md)** for extracting quality-controlled phage genomes from metagenomic read datasets. This pipeline uses current benchmarked tools for phage identification. It includes commented out code, brief theoretical explanations of each step, and definitions for various command line options. There are certainly steps that can be taken to further optimize this pipeline and it is not perfect, but I aim here to present a basic foundation. 

I hope this helps make phage discovery more accessible for microbiologists getting started in bioinformatics!

## Info before we start
- I am working under the assumption that you have basic skills for navigating in command line.
- I am also using Oxford Nanopore long reads in this pipeline. Long reads help assemble higher quality genomes.
    - Read QC/trimming and metagenome assembly will use different tools if you are using short reads (ie. FastP for QC and MEGAHIT/MetaSPAdes for assembly). The remainder of this pipeline will still work with short read datasets. 
- I am also working under the assumption that your reads have already been trimmed (if using Oxford Nanopore, their basecaller, Dorado, automatically trims adapter sequences)
    - Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic) can be used to trim short reads
- This pipeline can also be applied to identify prophage sequences and plasmids. Although some downstream analysis does not apply for plasmids (ie. phage-host prediction)

## References
Check out my preprint that uses this pipeline to investigate prophage induction by antibiotics in human gut metagenomes!
https://doi.org/10.64898/2026.03.11.711154 

- Wang, H., Sun, C., Li, Y. et al. Complementary insights into gut viral genomes: a comparative benchmark of short- and long-read metagenomes using diverse assemblers and binners. Microbiome 12, 260 (2024). https://doi.org/10.1186/s40168-024-01981-z
- Shen, Wei, Botond Sipos, and Liuyang Zhao. 2024. “ SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing.” iMeta 3, e191. https://doi.org/10.1002/imt2.191
- Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
- Wouter De Coster, Rosa Rademakers, NanoPack2: population-scale evaluation of long-read sequencing data, Bioinformatics, Volume 39, Issue 5, May 2023, btad311, https://doi.org/10.1093/bioinformatics/btad311
- Mikhail Kolmogorov, Derek M. Bickhart, Bahar Behsaz, Alexey Gurevich, Mikhail Rayko, Sung Bong Shin, Kristen Kuhn, Jeffrey Yuan, Evgeny Polevikov, Timothy P. L. Smith and Pavel A. Pevzner "metaFlye: scalable long-read metagenome assembly using repeat graphs", Nature Methods, 2020 doi:10.1038/s41592-020-00971-x
- Nayfach, S., Camargo, A.P., Schulz, F. et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat Biotechnol 39, 578–585 (2021). https://doi.org/10.1038/s41587-020-00774-7
- Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. — Nature Biotechnology (2023), DOI: 10.1038/s41587-023-01953-y.
- Bin Jang, H., Bolduc, B., Zablocki, O., Kuhn, J. H., Roux, S., Adriaenssens, E. M., … Sullivan, M. B. (2019). Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks. Nature Biotechnology. https://doi.org/10.1038/s41587-019-0100-8
- Sakdinan, Ben. et al. BioRxiv. "Species-specific prophage induction by ciprofloxacin in human gut metagenomes".

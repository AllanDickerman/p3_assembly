# Genome Assembly Service

## Overview

The bacterial Genome Assembly Service allows single or multiple assemblers to be invoked to compare results. Several assembly workflows or "strategies" are available that have been tuned to fit certain data types or desired analysis criteria such as throughput or rigor. 
Short reads (Illumina) or long reads (Nanopore or PacBio), or a mixture, can be submitted.
Available assembler programs include:
* SPAdes
: General purpose assembler for Illumina reads with additional modes for special cases (metagenomes, plasmids, single-cell)
* Unicycler
: Uses multiple SPAdes runs with varying parameters and adds steps to extend and merge contigs. It is particularly effective if both Illumina and long-reads are submitted.
* Flye
: Modern long-read assembler.
* Canu
: Older long-read assembler.

An "auto" option lets the pipeline select the assembler for the job:
* Unicycler if just Illumina reads or Illumina plus long reads are submitted.
* Flye if just long reads are submitted.

Prior to assembly, the following two manipulations of Illumina read data are optional:
* Trimming low-quality ends by TrimGalore.
* Normalizing by BBNorm (downsampling to an average depth of 100X attempting to smooth variation in coverage).

Once the assembly process has started by clicking the Assemble button, the genome is queued as a "job" for the Assembly Service to process, and will increment the count in the Jobs information box on the bottom right of the page. Once the assembly job has successfully completed, the output file will appear in the workspace, available for use in the BV-BRC comparative tools and downloaded if desired.

The principal output is the assembly in a single file in fasta format.
A report is provided in HTML format describing the input data, any transformations, the assembly process and the output. A Bandage plot and a Quast statistical summary are incorporated.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication:
* [GenomeAssembly2](app_specs/GenomeAssembly2.md)


## See also

* [Genome Assembly Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/genome_assembly_service.html)
* [Genome Assembly Service](https://www.bv-brc.org/app/Assembly2)
* [Genome Assembly Service Tutorial](https://www.bv-brc.org/docs//tutorial/genome_assembly/assembly.html)



## References

1.	Wick, R.R., et al., Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS computational biology, 2017. 13(6): p. e1005595.
2.	Bankevich, A., et al., SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology, 2012. 19(5): p. 455-477.
3.	Koren, S., et al., Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome research, 2017. 27(5): p. 722-736.
4.	Nurk, S., et al., metaSPAdes: a new versatile metagenomic assembler. Genome research, 2017. 27(5): p. 824-834.
5.	Antipov, D., et al., plasmidSPAdes: assembling plasmids from whole genome sequencing data. bioRxiv, 2016: p. 048942.
6.	Krueger, F., Trim Galore: a wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries. URL http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ (Date of access: 28/04/2016), 2012.
7.	Wick, R.R., et al., Bandage: interactive visualization of de novo genome assemblies. Bioinformatics, 2015. 31(20): p. 3350-3352.
8.	Gurevich, A., et al., QUAST: quality assessment tool for genome assemblies. Bioinformatics, 2013. 29(8): p. 1072-1075.
9. Kolmogorov, M. et al., "Assembly of Long Error-Prone Reads Using Repeat Graphs", Nature Biotechnology, 2019 doi:10.1038/s41587-019-0072-8
10. BBNorm: https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/44377-introducing-bbnorm-a-read-normalization-and-error-correction-tool


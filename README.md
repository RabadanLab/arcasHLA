# arcasHLA: high resolution HLA typing from RNA seq #

## Installation ##
To install arcasHLA, download or clone this repository.

### Requirements ###

- [Samtools](http://www.htslib.org/)
- [bedtools](http://bedtools.readthedocs.io/)
- [pigz](https://zlib.net/pigz/)
- [Kallisto v0.44.0](https://pachterlab.github.io/kallisto/)
- Python 3.6
    - NumPy
    - Pandas
    - [Biopython](https://biopython.org/wiki/Download)

## Usage ##

To see the list of available tools, simply enter `arcasHLA`. To view the required and optional arguments for any of the tools enter `arcasHLA [command] -h`.

### Extract reads ###
arcasHLA takes sorted bam files and extracts chromosome 6 reads and related HLA sequences. If the bam file is not indexed, this tool will run samtools index before extracting reads.

    arcasHLA extract [options][--paired] /path/to/sample.bam 
    
Output: `sample.1.fq.gz`, `sample.2.fq.gz`

### Genotype - complete ###

To predict the most likely genotype (complete alleles), input the fastqs produced by `extract`.

    arcasHLA genotype [options] /path/to/sample.1.fq.gz /path/to/sample.2.fq.gz
  
Output: `sample.json`, `sample.em.json`, `sample.genotype.json`
  
If you have previously run `genotype` on a sample, you can run `genotype` again directly from `sample.json` to retype without aligning with Kallisto again. This is useful if you want to try different populations, genes and other parameters.

    arcasHLA genotype [options] /path/to/sample.json

### Genotype - partial ###
Following genotyping, partial alleles can be predicted. This requires aligning the reads to an alternate, partial allele reference. The `sample.genotype.json` file from the previous step is required.

    arcasHLA partial [options] -G /path/to/sample.genotype.json /path/to/sample.1.fq.gz /path/to/sample.2.fq.gz
    
Output: `sample.partial.json`, `sample.genotype.partial.json`

### Merge jsons ###
To make analysis easier, this command will merge all jsons produced by genotyping. Simply provide the input directory and the run name.

    arcasHLA table [options] -i /path/to/indir --run run_name

### Change reference ###
To update the reference to the latest IMGT/HLA version, run

    arcasHLA reference --update
    
If you are running multiple tools to type HLAs, it can be helpful to use the same version of IMGT/HLA. You can select the version you like using the commithash from [the IMGT/HLA github](https://github.com/ANHIG/IMGTHLA/commits/Latest).

    arcasHLA reference --version [commithash]
# arcasHLA: high resolution HLA typing from RNA seq #

arcasHLA performs high resolution genotyping for HLA class I and class II genes from RNA sequencing, supporting both paired and single-end samples.

### Dependencies ###
arcasHLA requires the following utilities:
- coreutils

Make sure the following programs are in your `PATH`:
- [Samtools v1.19](http://www.htslib.org/)
- [bedtools v2.27.1](http://bedtools.readthedocs.io/)
- [pigz v2.3.1](https://zlib.net/pigz/)
- [Kallisto v0.44.0](https://pachterlab.github.io/kallisto/)
- Python 3.6

arcasHLA requires the following Python modules:
- [Biopython v1.77 (or lower)](https://biopython.org/wiki/Download)
- NumPy
- SciPy
- Pandas

### Test ###
In order to test arcasHLA partial typing, we need to roll back the reference to an earlier version. First, fetch IMGT/HLA database version 3.24.0:
```
./arcasHLA reference --version 3.24.0
```
Extract reads:
```
./arcasHLA extract test/test.bam -o test/output -t 8 -v
```
Genotyping (no partial alleles):
```
./arcasHLA genotype test/output/test.extracted.1.fq.gz test/output/test.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o test/output -t 8 -v
```
Expected output in `test/output/test.genotype.json`:
```
{"A": ["A*01:01:01", "A*03:01:01"], 
 "B": ["B*39:01:01", "B*07:02:01"], 
 "C": ["C*08:01:01", "C*01:02:01"], 
 "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"], 
 "DQA1": ["DQA1*02:01:01", "DQA1*05:03"], 
 "DQB1": ["DQB1*02:02:01", "DQB1*06:09:01"], 
 "DRB1": ["DRB1*10:01:01", "DRB1*14:02:01"]}
```
Partial typing:
```
./arcasHLA partial test/output/test.extracted.1.fq.gz test/output/test.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -G test/output/test.genotype.json -o test/output -t 8 -v
```
Expected output in `test/output/test.partial_genotype.json`:
```
{"A": ["A*01:01:01", "A*03:01:01"], 
 "B": ["B*07:02:01", "B*39:39:01"],
 "C": ["C*08:01:01", "C*01:02:01"], 
 "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"], 
 "DQA1": ["DQA1*02:01:01", "DQA1*05:03"],
 "DQB1": ["DQB1*06:04:01", "DQB1*02:02:01"],
 "DRB1": ["DRB1*03:02:01", "DRB1*14:02:01"]}
```
Remember to update the HLA reference using the following command.
```
./arcasHLA reference --update
```

### Usage ###

To see the list of available tools, simply enter `arcasHLA`. To view the required and optional arguments for any of the tools enter `arcasHLA [command] -h`.

- `extract` : Extracts reads mapped to chromosome 6 and any HLA decoys or chromosome 6 alternates.
- `genotype` : Genotypes HLA alleles from extracted reads (no partial alleles).
- `partial` : Genotypes partial HLA alleles from extracted reads and output from `genotype` (optional).
- `reference` : Update, specify version or force rebuilding of HLA reference.
- `merge` : merge genotyping output for multiple samples into a single json file.

### Extract reads ###

arcasHLA takes sorted BAM files and extracts chromosome 6 reads and related HLA sequences. If the BAM file is not indexed, this tool will run samtools index before extracting reads. By default, `extract` outputs a single FASTQ file; use the `--paired` flag for paired-end samples.

    arcasHLA extract [options] /path/to/sample.bam 
    
Output: `sample.1.fq.gz`, `sample.2.fq.gz`

#### Options: ####
- `--single`          : single-end reads (default: False)                                                                             
- `--unmapped`        : include unmapped reads, recommended if the aligner used marks multimapping reads as unmapped (default: False) 
- `--log FILE`        : log file for run summary (default: sample.extract.log)                                                        
- `--o, --outdir DIR` : output directory (default: `.`)                                                                               
- `--temp DIR`        : temp directory (default: `/tmp`)                                                                              
- `--keep_files`      : keep intermediate files (default: False)                                                                      
- `-t, --threads INT` : number of threads (default: 1)                                                                                
- `-v, --verbose`     : verbosity (default: False)                

### Genotype ###

#### From FASTQs ####
To predict the most likely genotype (no partial alleles), input the FASTQs produced by `extract` or the original FASTQs with all reads (experimental - use with caution).

```
arcasHLA genotype [options] /path/to/sample.1.fq.gz /path/to/sample.2.fq.gz
```

Output: `sample.alignment.p`, `sample.em.json`, `sample.genotype.json`

#### From intermediate alignment file ####  
If you have previously run `genotype` on a sample, you can run `genotype` again directly from `sample.alignment.p` to retype without aligning with Kallisto again. This is useful if you want to try different populations, genes and other parameters.
```
arcasHLA genotype [options] /path/to/sample.alignment.p
``` 
#### Example `.genotype.json` ####

```
{'A': ['A*01:01:01', 'A*29:02:01'],
 'B': ['B*08:01:01', 'B*44:03:01'],
 'C': ['C*07:01:01', 'C*16:01:01'],
 'DQA1': ['DQA1*02:01:01', 'DQA1*05:01:01'],
 'DQB1': ['DQB1*02:01:01', 'DQB1*02:02:01'],
 'DRB1': ['DRB1*03:01:01', 'DRB1*07:01:01']}
```

#### Options ####
- `-g, --genes GENES`       : comma separated list of HLA genes (ex. A,B,C,DQA1,DQB1,DRB1)
- `-p, --population POPULATION`  : sample population, options are asian_pacific_islander, black, caucasian, hispanic, native_american and prior (default: Prior)
- `--min_count INT`   : minimum gene read count required for genotyping (default: 75)
- `--tolerance FLOAT` : convergence tolerance for transcript quantification (default: 10e-7)
- `--max_iterations INT` : maximmum number of iterations for transcript quantification (default: 1000)
- `--drop_iterations INT` : number of iterations before dropping low support alleles, a lower number of iterations is recommended for single-end and low read count samples (default: paired - 10, single - 4)
- `--drop_threshold FLOAT` : proportion of maximum abundance an allele needs to not be dropped (default: 0.1)
- `--zygosity_threshold FLOAT` : threshold for ratio of minor to major allele nonshared count to determine zygosity (default: 0.15)
- `--log FILE`        : log file for run summary (default: `sample.genotype.log`)                                                        
- `--o, --outdir DIR` : output directory (default: `.`)                                                                               
- `--temp DIR`        : temp directory (default: `/tmp`)                                                                              
- `--keep_files`      : keep intermediate files (default: False)                                                                      
- `-t, --threads INT` : number of threads (default: 1)                                                                                
- `-v, --verbose`     : verbosity (default: False)
- `--single`          : Include flag to indicate if single-end FASTQs (paired-end if missing)
- `-l, --avg`         : Estimated average fragment length for single-end reads (default: 200)
- `-s, --std`         : Estimated standard deviation of fragment length (default: 20)


### Genotype - partial (optional) ###
Following genotyping, partial alleles can be predicted. This requires aligning the reads to an alternate, partial allele reference. The `sample.genotype.json` file from the previous step is required.

```
arcasHLA partial [options] -G /path/to/sample.genotype.json /path/to/sample.1.fq.gz /path/to/sample.2.fq.gz
```
   
Output: `sample.partial_alignment.p`, `sample.partial_genotype.json`

The options for partial typing are the same as genotype. Partial typing can be run from the intermediate alignment file.
 
### Merge jsons ###
To make analysis easier, this command will merge all jsons produced by genotyping into a single table. All `.genotype.json` files will be merged into a single `run.genotypes.tsv` file and all `.partial_genotype.json` files will be merged into `run.partial_genotypes.tsv`. In addition, HLA locus read counts and relative abundance produced by alignment will be merged into a single tsv file.
```
arcasHLA merge [options]
```
#### Options ####
- `--run RUN` : run name
- `--i, --indir DIR` : input directory (default: `.`)     
- `--o, --outdir DIR` : output directory (default: `.`)                                                                  
- `-v, --verbose`     : toggle verbosity

### Convert HLA nomenclature ###
arcasHLA convert changes alleles in a tsv file from its input form to a specified grouped nomenclature (P-group or G-group) or a specified number of fields (i.e. 1, 2 or 3 fields in resolution). This file can be produced by arcasHLA merge or any tsv following the same structure:

| subject      	| A1         	| A2         	| B1         	| B2         	| C1         	| C2         	|
|--------------	|------------	|------------	|------------	|------------	|------------	|------------	|
| subject_name 	| A*01:01:01 	| A*01:01:01 	| B*07:02:01 	| B*07:02:01 	| C*04:01:01 	| C*04:01:01 	|

P-group (alleles sharing the same amino acid sequence in the antigen-binding region) and G-group (alleles sharing the same base sequence in the antigen-binding region) can only be reduced to 1-field resolution as alleles with differing 2nd fields can be in the same group. By the same reasoning, P-group cannot be converted into G-group.

```
arcasHLA convert --resolution [resolution] genotypes.tsv
```
#### Options ####
- `-r, --resolution RESOLUTION` : output resolution (1, 2, 3) or grouping (g-group, p-group)
- `-o, --outfile FILE` : output file (default: `./run.resolution.tsv`)
- `-f, --force` : force conversion for grouped alleles even if it results in loss of resolution
- `-v, --verbose`     : toggle verbosity

### Change reference ###
To update the reference to the latest IMGT/HLA version, run

```
arcasHLA reference --update
```
If you are running multiple tools to type HLAs, it can be helpful to use the same version of IMGT/HLA. You can select the version you like using the commithash from the [IMGT/HLA Github](https://github.com/ANHIG/IMGTHLA/commits/Latest).

```
arcasHLA reference --version [commithash]
```

If you suspect there is an issue  with the reference files, rebuild the reference with the following command
```
arcasHLA reference --rebuild
```

Note: if your reference was built with arcasHLA version <= 0.1.1 and you wish to change your reference to versions >= 3.35.0, it may be necessary to remove the IMGTHLA folder due to the need for Git Large File Storage to properly download hla.dat.

```
rm -rf dat/IMGTHLA
arcasHLA reference --update
```

#### Options ####
- `--update` : update to latest IMGT/HLA version
- `--version` : checkout IMGT/HLA version using commithash
- `--rebuild` : rebuild HLA database
- `-v, --verbose`     : verbosity (default: False)   

## Citations ##
Orenbuch R, Filip I, Comito D, et al (2019) arcasHLA: high resolution HLA typing from RNA seq. Bioinformatics doi:[10.1093/bioinformatics/btz474](http://dx.doi.org/10.1093/bioinformatics/btz474)

Orenbuch R, Filip I, Rabadan R (2020) HLA Typing from RNA Sequencing and Applications to Cancer. Methods Mol. Biol. doi: 10.1007/978-1-0716-0327-7_5 (https://link.springer.com/protocol/10.1007%2F978-1-0716-0327-7_5)


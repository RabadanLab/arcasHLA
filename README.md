### Dependencies ###

Install `arcasHLA` through bioconda with:
```
conda install arcas-hla -c bioconda -c conda-forge
conda activate arcas-hla
```
**Important**: Please include channels `bioconda` and `conda-forge` as above.

`arcasHLA` can also be installed through the [environment.yml](./environment.yml) file in this repo:
```
conda env create -f environment.yml
conda activate arcas-hla
```

### Test ###

**(Update 2023-09-29)**: The below tests are now implemented as a pytest [suite](./test/test_arcas_hla.py). You can run this locally by building the docker environment and running pytest. From the current directory:

```
docker build -t <image-name> -f Docker/Dockerfile .
docker run --rm -v /path/to/repo:/app <image-name> pytest
```
-----

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

arcasHLA takes sorted BAM files and extracts chromosome 6 reads and related HLA sequences. If the BAM file is not indexed, this tool will run samtools index before extracting reads. By default, `extract` outputs paired FASTQ files; use the `--single` flag for single-end samples.

    arcasHLA extract [options] /path/to/sample.bam 
    
Output: `sample.extracted.1.fq.gz`, `sample.extracted.2.fq.gz`

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

## Build Customized References ##

#### Input: arcasHLA genotypes.json ####
Customized references can be built from arcasHLA genotype outputs.
```
./arcasHLA customize genotypes.json -o ~/ref
```
#### Input: HLA tsv ####

Customized references can be built from a tab-separated file with the following structure:

| subject | A1      | A2      | B1      | B2      | C1      | C2      |
|---------|---------|---------|---------|---------|---------|---------|
| Example | A*01:01 | A*02:01 | B*07:01 | B*52:01 | C*04:01 | C*18:01 |

```
./arcasHLA customize hla.tsv -o ~/ref
```
#### Options: ####
```
usage: arcasHLA customize [options]

optional arguments:
  -h, --help            show this help message and exit

  -G , --genotype       comma-separated list of HLA alleles (e.g. A*01:01,A*11:01,...)
                        arcasHLA output genotype.json or genotypes.json
                        or tsv with format specified in README.md
  -s , --subject        subject name, only required for list of alleles
  -g , --genes          comma separated list of HLA genes
                        default: all
                        options: A, B, C, DMA, DMB, DOA, DOB, DPA1, DPB1, DQA1,
                        DQB1, DRA, DRB1, DRB3, DRB5, E, F, G, H, J, K, L

  --transcriptome TRANSCRIPTOME
                        transcripts to include besides input HLAs
                         options: full, chr6, none
                          default: full

  --resolution RESOLUTION
                        genotype resolution, only use >2 when typing performed with assay or Sanger sequencing
                          default: 2

  --grouping GROUPING   type/number of transcripts to include per allele
                         single - one 3-field resolution transcript per allele (e.g. A*01:01:01)
                        g-group - all transcripts with identical binding regions
                          default: protein group - all transcripts with identical protein types (2 fields the same)

  -o , --outdir         out directory

  --temp                temp directory

  --keep_files          keep intermediate files

  -t , --threads
  -v, --verbose
```

## Quantification ##
Note: if the reference was built with the `--chr6` flag, you should run `quant` with extracted chromosome 6 FASTQs (see `extract`).

```
./arcasHLA quant --ref /path/to/ref/sample FASTQ
```

Example:
```
./arcasHLA quant --ref ~/ref/Pt23 -t 8 -o /Volumes/quant/ /Volumes/fastq/Pt23_pre.1.fq.gz /Volumes/fastq/Pt23_pre.2.fq.gz
```

#### Options: ####
```
usage: arcasHLA quant [options] FASTQs

positional arguments:
  file               list of fastq files

optional arguments:
  -h, --help         show this help message and exit

  --sample SAMPLE    sample name
  --ref              arcasHLA quant_ref path (e.g. "/path/to/ref/sample")

  -o , --outdir      out directory

  --temp             temp directory

  --keep_files       keep intermediate files

  -l AVG, --avg AVG  Estimated average fragment length for single-end reads
                       default: 200

  -s STD, --std STD  Estimated standard deviation of fragment length for single-end reads
                       default: 20

  --single           Include flag if single-end reads. Default is paired-end.

  -t , --threads
  -v, --verbose
```
 

## Citations ##
* Orenbuch R, Filip I, Comito D, et al (2019) arcasHLA: high resolution HLA typing from RNA seq. Bioinformatics doi:[10.1093/bioinformatics/btz474](http://dx.doi.org/10.1093/bioinformatics/btz474)
* Orenbuch R, Filip I, Rabadan R (2020) HLA Typing from RNA Sequencing and Applications to Cancer. Methods Mol. Biol. doi: 10.1007/978-1-0716-0327-7_5 (https://link.springer.com/protocol/10.1007%2F978-1-0716-0327-7_5)
* Filip, I., Wang, A., Kravets, O. et al. Pervasiveness of HLA allele-specific expression loss across tumor types. Genome Med 15, 8 (2023). https://doi.org/10.1186/s13073-023-01154-x

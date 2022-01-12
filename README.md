# arcasHLA: high resolution HLA typing from RNA seq #

arcasHLA performs high resolution genotyping for HLA class I and class II genes from RNA sequencing, supporting both paired and single-end samples.

# arcasHLA quantification #

## Setup ##
```
git clone https://github.com/roseorenbuch/arcasHLA.git
cd arcasHLA
```
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
./arcasHLA quant --ref ~/ref/Pt23 -t 8 -o /Volumes/quant/ /Volumes/fastq/Pt23_pre.1.fq.gz /Volumes/f
astq/Pt23_pre.2.fq.gz
```

#### Options: ####
```
usage: arcasHLA quant [options] FASTQs

positional arguments:
  file             list of fastq files

optional arguments:
  -h, --help       show this help message and exit

  --sample SAMPLE  sample name
  --ref            arcasHLA quant_ref path (e.g. "/path/to/ref/sample")

  -o , --outdir    out directory

  --temp           temp directory

  --keep_files     keep intermediate files

  -t , --threads
  -v, --verbose
```

## Merge ##
Merge will create a tsv file containing all the quantification results ("run.quant.tsv").
```
./arcasHLA merge -i /Volumes/quant/ --run test -o ./
```
Merge will also now create a tsv file containing all gene counts when supplied a folder with ".genotype.log" files.

### Dependencies ###
Make sure the following programs are in your `PATH`:
- [Samtools](http://www.htslib.org/)
- [bedtools](http://bedtools.readthedocs.io/)
- [pigz](https://zlib.net/pigz/)
- [Kallisto v0.44.0](https://pachterlab.github.io/kallisto/)
- Python 3.6
- GNU Parallel

arcasHLA requires the following Python modules:
- [Biopython](https://biopython.org/wiki/Download)
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
./arcasHLA extract test/test.bam -o test/output --paired -t 8 -v
```
Complete typing:
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
Before further usage, remember to update to the 3.34.0.
```
./arcasHLA reference --version 3.34.0
```
At this time, cloning the latest version of IMGTHLA (3.35.0) results in a corrupted database due to an issue with GitHub's Large File Storage. To update the arcasHLA's reference to the current version, run the following commands.
```
curl https://media.githubusercontent.com/media/ANHIG/IMGTHLA/Latest/hla.dat > dat/IMGTHLA/hla.dat
./arcasHLA reference --rebuild --v
```

### Usage ###

To see the list of available tools, simply enter `arcasHLA`. To view the required and optional arguments for any of the tools enter `arcasHLA [command] -h`.

- `extract` : Extracts reads mapped to chromosome 6 and any HLA decoys or chromosome 6 alternates.
- `genotype` : Genotypes complete HLA alleles from extracted reads.
- `partial` : Genotypes partial HLA alleles from extracted reads and output from `genotype`.
- `reference` : Update, specify version or force rebuilding of HLA reference.
- `merge` : merge genotyping output for multiple samples into a single json file.

### Extract reads ###

arcasHLA takes sorted BAM files and extracts chromosome 6 reads and related HLA sequences. If the BAM file is not indexed, this tool will run samtools index before extracting reads. By default, `extract` outputs a single FASTQ file; use the `--paired` flag for paired-end samples.

    arcasHLA extract [options] /path/to/sample.bam 
    
Output: `sample.1.fq.gz`, `sample.2.fq.gz`

#### Options: ####
- `--paired`          : paired-end reads (default: False)                                                                             
- `--unmapped`        : include unmapped reads, recommended if the aligner used marks multimapping reads as unmapped (default: False) 
- `--log FILE`        : log file for run summary (default: sample.extract.log)                                                        
- `--o, --outdir DIR` : output directory (default: `.`)                                                                               
- `--temp DIR`        : temp directory (default: `/tmp`)                                                                              
- `--keep_files`      : keep intermediate files (default: False)                                                                      
- `-t, --threads INT` : number of threads (default: 1)                                                                                
- `-v, --verbose`     : verbosity (default: False)                

### Genotype - complete ###

#### From FASTQs ####
To predict the most likely genotype (complete alleles), input the FASTQs produced by `extract`.

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
- `--tolerance FLOAT` : convergence tolerance for transcript quantification (default: 10e-7)
- `--max_iterations INT` : maximmum number of iterations for transcript quantification (default: 1000)
- `--drop_iterations INT` : number of iterations before dropping low support alleles, a lower number of iterations is recommended for single-end and low read couunt samples (default: paired - 10, single - 4)
- `--drop_threshold FLOAT` : proportion of maximum abundance an allele needs to not be dropped (default: 0.1)
- `--zygosity_threshold FLOAT` : threshold for ratio of minor to major allele nonshared count to determine zygosity (default: 0.15)
- `--log FILE`        : log file for run summary (default: sample.genotype.log)                                                        
- `--o, --outdir DIR` : output directory (default: `.`)                                                                               
- `--temp DIR`        : temp directory (default: `/tmp`)                                                                              
- `--keep_files`      : keep intermediate files (default: False)                                                                      
- `-t, --threads INT` : number of threads (default: 1)                                                                                
- `-v, --verbose`     : verbosity (default: False)   

### Genotype - partial ###
Following genotyping, partial alleles can be predicted. This requires aligning the reads to an alternate, partial allele reference. The `sample.genotype.json` file from the previous step is required.

```
arcasHLA partial [options] -G /path/to/sample.genotype.json /path/to/sample.1.fq.gz /path/to/sample.2.fq.gz
```
   
Output: `sample.partial_alignment.p`, `sample.partial_genotype.json`

The options for partial typing are the same as complete. Partial typing, like complete, can be run from the intermediate alignment file.
 
### Merge jsons ###
To make analysis easier, this command will merge all jsons produced by genotyping. All `.genotype.json` files will be merged into a single `run.genotypes.json` file and all `.partial_genotype.json` files will be merged into `run.partial_genotypes.json`.
```
arcasHLA merge [options]
```
#### Options ####
- `--run RUN` : run name
- `--i, --indir DIR` : input directory (default: `.`)     
- `--o, --outdir DIR` : output directory (default: `.`)                                                                  
- `-v, --verbose`     : verbosity (default: False)   

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
#### Options ####
- `--update` : update to latest IMGT/HLA version
- `--version` : checkout IMGT/HLA version using commithash
- `--rebuild` : rebuild HLA database
- `-v, --verbose`     : verbosity (default: False)   

## Citation ##
R. Orenbuch, I. Filip, D. Comito, J. Shaman, I. Pe’er, and R. Rabadan. arcasHLA:
high resolution HLA typing from RNA seq. bioRxiv doi: [10.1101/479824](https://doi.org/10.1101/479824)

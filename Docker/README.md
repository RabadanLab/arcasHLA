# Container for arcasHLA #

Installs up-to-date versions of arcasHLA and all dependencies:

- arcasHLA
- [Git Large File Storage](https://github.com/git-lfs/git-lfs/wiki/Installation)
- coreutils
- [Samtools v1.19](http://www.htslib.org/)
- [bedtools v2.29.1](http://bedtools.readthedocs.io/)
- [pigz v2.3.1](https://zlib.net/pigz/)
- [Kallisto v0.44.0](https://pachterlab.github.io/kallisto/)
- python 3.6
- [Biopython v1.77](https://biopython.org/wiki/Download)
- numpy
- scipy
- pandas
- pytest

### Build ###
In order to use this arcasHLA container, install Docker and build in this directory:
```
docker build -t <image_name> .
```
### Run ###
Interactively ("image_name" is as above):
```
docker run -it --entrypoint bash -v <path/to/files>:<docker/path/to/files> <image_name>
```
Noninteractively ("image_name" is as above), e.g. 'arcasHLA extract':
```
docker run \
	--rm -v <path/to/files>:<docker/path/to/files> \
	<image_name> \
	arcasHLA extract --o docker/path/to/files/out_dir \
	[other options] docker/path/to/files/sample.bam & 

```

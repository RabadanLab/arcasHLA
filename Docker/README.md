# Container for arcasHLA #

Installs up-to-date versions of arcasHLA and all dependencies from `environment.yml`.

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

FROM continuumio/miniconda3:latest

WORKDIR /app

# Install Git LFS
RUN apt-get update && \
    apt-get install -y git curl
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install -y git-lfs

RUN git clone -b buendia/ref_patch https://github.com/RabadanLab/arcasHLA.git
RUN conda env create -f /app/arcasHLA/environment.yml
RUN echo "source activate arcas-hla" > ~/.bashrc
ENV PATH /opt/conda/envs/arcas-hla/bin:$PATH

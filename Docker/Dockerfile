FROM continuumio/miniconda3:latest

WORKDIR /app

# Install Git LFS
RUN apt-get update && \
    apt-get install -y git curl pigz
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install -y git-lfs

RUN git clone https://github.com/RabadanLab/arcasHLA.git
RUN conda env create -f /app/arcasHLA/environment.yml
RUN echo "source activate arcas-hla" > ~/.bashrc
ENV PATH /opt/conda/envs/arcas-hla/bin:$PATH
ENV PATH="${PATH}:/app/arcasHLA/"

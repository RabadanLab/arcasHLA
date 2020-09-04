# work from latest LTS ubuntu release
FROM ubuntu:18.04

# set the environment variables
ENV kallisto_version 0.44.0
ENV samtools_version 1.9
ENV bedtools_version 2.29.2
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# run update and install necessary tools from package manager
RUN apt-get update -y && apt-get install -y \
    build-essential \
    cmake \
    zlib1g-dev \
    libhdf5-dev \
    libnss-sss \
    curl \
    autoconf \
    bzip2 \
    python3-dev \
    python3-pip \
    python3-biopython \
    pigz \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    unzip

# install python libraries
RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install pandas

# install kallisto
RUN mkdir -p /usr/bin/kallisto \
    && curl -SL https://github.com/pachterlab/kallisto/archive/v${kallisto_version}.tar.gz \
    | tar -zxvC /usr/bin/kallisto

RUN mkdir -p /usr/bin/kallisto/kallisto-${kallisto_version}/build
RUN cd /usr/bin/kallisto/kallisto-${kallisto_version}/build && cmake ..
RUN cd /usr/bin/kallisto/kallisto-${kallisto_version}/ext/htslib && autoreconf
RUN cd /usr/bin/kallisto/kallisto-${kallisto_version}/build && make
RUN cd /usr/bin/kallisto/kallisto-${kallisto_version}/build && make install

# install samtools
WORKDIR /usr/bin/
RUN curl -SL https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    > samtools-${samtools_version}.tar.bz2
RUN tar -xjvf samtools-${samtools_version}.tar.bz2
WORKDIR /usr/bin/samtools-${samtools_version}
RUN ./configure
RUN make
RUN make install

# install bedtools
WORKDIR /usr/bin
RUN curl -SL https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools-${bedtools_version}.tar.gz \
    > bedtools-${bedtools_version}.tar.gz
RUN tar -xzvf bedtools-${bedtools_version}.tar.gz
WORKDIR /usr/bin/bedtools2
RUN make
RUN ln -s /usr/bin/bedtools2/bin/bedtools /usr/bin/bedtools

# git lfs
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install -y git-lfs
RUN git lfs install --system --skip-repo

# install arcasHLA
WORKDIR /home/
RUN git clone https://github.com/roseorenbuch/arcasHLA-quant.git arcasHLA-quant

ENV PATH="${PATH}:/home/arcasHLA-quant/"

RUN arcasHLA reference --update

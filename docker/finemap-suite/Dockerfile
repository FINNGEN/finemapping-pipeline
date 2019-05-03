FROM google/cloud-sdk:slim

LABEL maintainer="mkanai@broadinstitute.org"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y wget make zlib1g-dev g++ libopenblas-dev gfortran bedtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# HTSLIB (for bgzip)
ENV HTSLIB_VERSION 1.9
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --disable-bz2 --disable-lzma && \
    make && make install && \
    rm -r ~/htslib-${HTSLIB_VERSION}*

# LDSTORE
ENV LDSTORE_VERSION v2.0b
ADD ./docker/finemap-suite/ldstore_v2.0b_x86_64 /usr/bin/ldstore

# FINEMAP
ENV FINEMAP_VERSION v1.3.1
RUN wget --quiet http://www.christianbenner.com/finemap_${FINEMAP_VERSION}_x86_64.tgz -P ~ && \
    tar xf ~/finemap_${FINEMAP_VERSION}_x86_64.tgz -C ~ && \
    cp ~/finemap_${FINEMAP_VERSION}_x86_64/finemap_${FINEMAP_VERSION}_x86_64 /usr/bin/finemap && \
    rm -r ~/finemap_${FINEMAP_VERSION}_x86_64*

# Add scripts
RUN pip install --no-cache numpy scipy pandas pybedtools
ENV PYTHONPATH /home:${PYTHONPATH}
COPY python/make_finemap_inputs.py /usr/local/bin/
COPY python/utils.py /home/

VOLUME [ "/root/.config/gcloud" ]
CMD [ "/bin/bash" ]

FROM ubuntu:22.04
#FROM google/cloud-sdk:slim
LABEL maintainer="mkanai@broadinstitute.org"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y curl gnupg2 lsb-release
# Install gcsfuse
RUN export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s` && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    echo "deb https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update --fix-missing && \
    apt-get install -y wget gawk make zlib1g-dev g++ libopenblas-dev gfortran bedtools python3-pip gcsfuse google-cloud-cli && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# HTSLIB (for bgzip)
ENV HTSLIB_VERSION 1.10.2
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --disable-bz2 --disable-lzma && \
    make && make install && \
    rm -r ~/htslib-${HTSLIB_VERSION}*

# LDSTORE
# finemap-suite:0.1 -- v2.0b (on the private website)
# finemap-suite:0.3 -- v2.0b_0504
# finemap-suite:0.4 -- v2.0b_0508
# finemap-suite:0.5 -- v2.0b_0510
# finemap-suite:1.2 -- v2.0 (public)
ENV LDSTORE_VERSION v2.0
# ADD ./docker/finemap-suite/ldstore_${LDSTORE_VERSION}_x86_64 /usr/bin/ldstore
RUN wget --quiet http://www.christianbenner.com/ldstore_${LDSTORE_VERSION}_x86_64.tgz -P ~ && \
    tar xf ~/ldstore_${LDSTORE_VERSION}_x86_64.tgz -C ~ && \
    cp ~/ldstore_${LDSTORE_VERSION}_x86_64/ldstore_${LDSTORE_VERSION}_x86_64 /usr/bin/ldstore && \
    rm -r ~/ldstore_${LDSTORE_VERSION}_x86_64*

# FINEMAP
# finemap-suite:0.1 -- v1.3.1 (on the website)
# finemap-suite:0.2 -- v1.4_0503
# finemap-suite:0.4 -- v1.4_0508
# finemap-suite:0.5 -- v1.4b_0510
# finemap-suite:1.1 -- v1.4 (public)
ENV FINEMAP_VERSION v1.4
# ADD ./docker/finemap-suite/finemap_${FINEMAP_VERSION}_x86_64 /usr/bin/finemap
RUN wget --quiet http://www.christianbenner.com/finemap_${FINEMAP_VERSION}_x86_64.tgz -P ~ && \
    tar xf ~/finemap_${FINEMAP_VERSION}_x86_64.tgz -C ~ && \
    cp ~/finemap_${FINEMAP_VERSION}_x86_64/finemap_${FINEMAP_VERSION}_x86_64 /usr/bin/finemap && \
    rm -r ~/finemap_${FINEMAP_VERSION}_x86_64*

# for pysam
ENV HTSLIB_CONFIGURE_OPTIONS "--disable-bz2 --disable-lzma"
# Add scripts
RUN python3 -m pip install --no-cache cython && \
    python3 -m pip install --no-cache numpy scipy pandas pybedtools==0.8.2
ENV PYTHONPATH /home:${PYTHONPATH}
COPY python/make_finemap_inputs.py python/filter_and_summarize.py /usr/local/bin/
COPY python/utils.py /home/
RUN chmod +x /usr/local/bin/make_finemap_inputs.py

VOLUME [ "/root/.config/gcloud" ]
CMD [ "/bin/bash" ]

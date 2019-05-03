FROM r-base:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Modified FROM google/cloud-sdk:slim
ARG CLOUD_SDK_VERSION=244.0.0
ENV CLOUD_SDK_VERSION=$CLOUD_SDK_VERSION

ARG INSTALL_COMPONENTS
ENV PATH "$PATH:/opt/google-cloud-sdk/bin/"
RUN apt-get update -qqy && apt-get install -qqy \
    curl \
    gcc \
    python-dev \
    python-pip \
    apt-transport-https \
    lsb-release \
    openssh-client \
    git \
    gnupg \
    && pip install -U crcmod && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=${CLOUD_SDK_VERSION}-0 $INSTALL_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud --version

# Install additional components
RUN apt-get update --fix-missing && \
    apt-get install -y libopenblas-dev libssl-dev libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# SUSIE
# dev version: 0.8.1.0521
ENV SUSIE_VERSION e81a5ca
RUN pip install --no-cache numpy pandas && \
    Rscript -e 'install.packages(c("devtools", "argparse", "reticulate", "dplyr", "data.table", "here", "R.utils"), noCache = TRUE, repos="https://cloud.r-project.org")' && \
    Rscript -e "Sys.setenv(TAR = \"/bin/tar\"); devtools::install_github(\"stephenslab/susieR@${SUSIE_VERSION}\")"

# Add scripts
COPY R/run_susieR.R /usr/local/bin/

VOLUME [ "/root/.config/gcloud" ]
CMD [ "/bin/bash" ]

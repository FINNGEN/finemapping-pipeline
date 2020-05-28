FROM r-base:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install additional components
RUN apt-get update --fix-missing && \
    apt-get install -y python-dev python-pip libxml2-dev libopenblas-dev libssl-dev libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# SUSIE
# dev version: 0.9.1.0.894ba2f
ENV SUSIE_VERSION 894ba2f
RUN Rscript -e 'install.packages(c("devtools", "argparse", "dplyr", "stringr", "data.table", "here", "R.utils"), noCache = TRUE, repos="https://cloud.r-project.org")' && \
    Rscript -e "Sys.setenv(TAR = \"/bin/tar\"); devtools::install_github(\"stephenslab/susieR@${SUSIE_VERSION}\")"

# Add scripts
COPY ./R/run_susieR.R /usr/local/bin/
RUN chmod +x /usr/local/bin/run_susieR.R

VOLUME [ "/root/.config/gcloud" ]
CMD [ "/bin/bash" ]

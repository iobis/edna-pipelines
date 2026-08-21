FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
        pandoc \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
        rmarkdown \
        knitr \
        DT \
        dplyr \
        tidyr \
        tibble \
        purrr \
        plotly

WORKDIR /work

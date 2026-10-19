# Computational environment for the full pipeline: R and package versions from
# renv.lock, plus Quarto for the reports. The image holds the environment; the
# repository is mounted at run time, so data and results land on the host.
#
# Build: docker build -t csf-proteomics .
# Run:   docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
#          -v "$PWD":/project csf-proteomics

FROM rocker/r-ver:4.6.1

ARG QUARTO_VERSION=1.10.18

# Package library outside /project, so mounting the repository does not hide it.
# Binary CRAN packages from Posit Package Manager (Ubuntu noble) speed up restore.
ENV RENV_PATHS_LIBRARY=/opt/renv/library \
    RENV_PATHS_CACHE=/opt/renv/cache \
    RENV_CONFIG_REPOS_OVERRIDE=https://p3m.dev/cran/__linux__/noble/latest

# System libraries required by CRAN binaries (fs, stringi, graphics devices)
# and by Bioconductor packages compiled from source
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates unzip git \
      libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
      libuv1-dev libicu-dev \
      libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libjpeg-dev libtiff-dev libwebp-dev \
    && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL -o /tmp/quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && dpkg -i /tmp/quarto.deb \
    && rm /tmp/quarto.deb

WORKDIR /project

# Restore packages first, so code changes do not invalidate this layer
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/
RUN R -e "renv::restore(prompt = FALSE)" \
    && chmod -R a+rX /opt/renv

COPY . .

CMD ["bash", "run_pipeline.sh"]

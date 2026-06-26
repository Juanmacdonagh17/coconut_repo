# ---------------------------------
# Builder
# ---------------------------------
FROM continuumio/miniconda3:25.3.1-1 AS builder

LABEL org.opencontainers.image.title="Coconut Builder Image" \
      org.opencontainers.image.description="Builder stage for Coconut application with Conda environment" \
      org.opencontainers.image.authors="Juan Mac Donagh (SBG-UNQ), Gustavo Parisi (SBG-UNQ)" \
      org.opencontainers.image.source="https://github.com/Juanmacdonagh17/coconut_repo" \
      org.opencontainers.image.vendor="SBG-UNQ" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.created="2025-06-10" \
      maintainer.email="macjuan17@gmail.com" \
      stage.purpose="Build environment for Coconut application" \
      build.dependencies="miniconda3,gcc,libcurl" \
      build.target="coconut executable"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN conda install -c conda-forge -y \
        libcurl \
        curl \
        gcc_linux-64 \
        gxx_linux-64 && \
    conda clean -afy

ENV SCRIPT_DIR="/usr/src/coconut" \
    COCONUT_SOURCE="/usr/src/coconut/coconut.c" \
    CONDA_GCC="/opt/conda/bin/x86_64-conda-linux-gnu-gcc" \
    CONDA_PREFIX="/opt/conda"

WORKDIR $SCRIPT_DIR
COPY *.c *.h ./

RUN conda run $CONDA_GCC -o coconut "$COCONUT_SOURCE" \
    -I"$CONDA_PREFIX/include" \
    -L"$CONDA_PREFIX/lib" \
    -Wl,-rpath,"$CONDA_PREFIX/lib" \
    -lcurl \
    -lm

# ---------------------------------
# Runner
# ---------------------------------
FROM debian:bookworm AS runner

LABEL org.opencontainers.image.title="Coconut" \
      org.opencontainers.image.description="A suite for transcripts, codon usage and protein structure analysis" \
      org.opencontainers.image.url="https://github.com/Juanmacdonagh17/coconut_repo" \
      org.opencontainers.image.documentation="https://github.com/Juanmacdonagh17/coconut_repo" \
      org.opencontainers.image.authors="Juan Mac Donagh (SBG-UNQ), Gustavo Parisi (SBG-UNQ)" \
      org.opencontainers.image.source="https://github.com/Juanmacdonagh17/coconut_repo" \
      org.opencontainers.image.vendor="SBG-UNQ" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.created="2025-06-10" \
      maintainer.email="macjuan17@gmail.com" \
      stage.purpose="Production runtime for Coconut application" \
      runtime.dependencies="libcurl4" \
      security.features="non-root-user" \
      base.image="debian:bookworm"

# Create non-root user and install dependencies in a single layer
RUN useradd -ms /bin/bash -u 10001 cocouser && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        libcurl4 \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /home/cocouser

# Copy only the necessary files
COPY --from=builder --chown=cocouser:cocouser /usr/src/coconut/coconut /home/cocouser/coconut

# Switch to non-root user
USER cocouser:cocouser

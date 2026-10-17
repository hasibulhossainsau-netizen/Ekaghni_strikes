# syntax=docker/dockerfile:1

FROM python:3.10-slim-bookworm AS python-dependencies

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        libffi-dev \
        libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY requirements-docker.txt .
RUN python -m pip wheel \
    --wheel-dir /wheels \
    --extra-index-url https://download.pytorch.org/whl/cpu \
    -r requirements-docker.txt

FROM python:3.10-slim-bookworm AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    HOME=/tmp \
    PATH="/usr/local/bin:${PATH}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        bash \
        coreutils \
        gawk \
        libgomp1 \
        make \
        mafft \
        ncbi-blast+ \
        tini \
    && rm -rf /var/lib/apt/lists/*

COPY --from=python-dependencies /wheels /wheels
COPY requirements-docker.txt /tmp/requirements-docker.txt
RUN python -m pip install \
        --no-index \
        --find-links=/wheels \
        -r /tmp/requirements-docker.txt \
    && rm -rf /wheels /tmp/requirements-docker.txt

WORKDIR /workspace
COPY Makefile README.md Snakefile config.yaml ./
COPY config/ ./config/
COPY envs/ ./envs/
COPY scripts/ ./scripts/
COPY data/fasta/test_sample.fasta ./data/fasta/test_sample.fasta

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["make", "smoke-test", "mock-test"]

FROM ghcr.io/astral-sh/uv:python3.13-bookworm

# Install the project into `/app`
WORKDIR /app

RUN apt-get update && apt-get install -y \
    g++ \
    gcc \
    zlib1g-dev \
    ca-certificates \
    tar \
    bedtools

COPY . .

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []
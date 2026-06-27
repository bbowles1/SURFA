# Example from https://github.com/astral-sh/uv-docker-example/blob/main/Dockerfile
FROM ghcr.io/astral-sh/uv:python3.13-bookworm

# Install the project into `/app`
WORKDIR /app

RUN apt-get update && apt-get install -y \
    g++ \
    gcc \
    zlib1g-dev \
    ca-certificates \
    tar \
    bedtools \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Install dependencies first (without the project itself) for optimal layer caching.
# Mounts pyproject.toml and uv.lock without copying them into the image at this stage.
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-install-project --no-dev

# Copy source and install the project itself.
# Kept as a separate layer so dependency installation above is only
# re-run when pyproject.toml / uv.lock change, not on every code edit.
COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:$PATH"

# Make test script(s) executable
RUN chmod +x tests/integration/build_test_db.sh

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []
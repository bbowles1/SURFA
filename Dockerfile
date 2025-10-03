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
    # dev dependencies below
    less \
    vim \
    # npm dependencies below
    curl \
    && curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Install the project's dependencies using the lockfile and settings
# essentially uses uv sync with explicit pyproject, lock files
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-install-project --no-dev

# Then, add the rest of the project source code and install it
# Installing separately from its dependencies allows optimal layer caching
ADD . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:$PATH"

# make test script(s) executable
RUN chmod +x tests/build_test_db.sh

# Install NPM dependencies
RUN npm install

# expose flask port
EXPOSE 5000

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []

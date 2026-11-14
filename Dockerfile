# abscore in a container.
#
# Why this exists: ANARCI assigns antibody numbering by shelling out to HMMER's
# `hmmscan`, a compiled C program. A Python lockfile cannot install it, and no
# native Windows build of HMMER exists, so on Windows the tool has to run
# inside Linux -- either WSL2 or a container.
#
# On macOS or Linux the container is optional; `brew install hmmer` or
# `apt-get install hmmer` plus `uv sync --frozen` is simpler. See the README.
#
#   docker build -t abscore .                     # CLI, the default stage
#   docker build -t abscore-app --target app .    # the Streamlit app
#   docker build --target test .                  # runs the suite on Linux, fails the
#                                                 # build if anything is red
#
#   docker run --rm abscore score --seq EVQLVESGG...
#   docker run --rm abscore validate
#   cat chain.fasta | docker run --rm -i abscore score
#   docker run --rm -v "$PWD/out:/out" abscore score --seq EVQ... --save /out/report.html
#   docker run --rm -p 8501:8501 abscore-app
#
# Size. The default image carries the CLI only. The Streamlit extra pulls
# pyarrow, pydeck, altair and pillow -- around 250 MB of packages that no score
# depends on -- so the app lives in its own stage and nobody scoring a sequence
# pays for it. The build is also split: dependencies are resolved in a builder
# stage and the finished virtualenv is copied into a clean runtime, so the
# image holds one copy of it rather than one per `uv sync`.

# ---------------------------------------------------------------- builder ---
FROM python:3.11-slim AS builder

# A pinned uv, taken from its own published image rather than a curl-to-shell
# script: the version is then visible here and reproducible.
COPY --from=ghcr.io/astral-sh/uv:0.12.5 /uv /usr/local/bin/uv

ENV UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_DOWNLOADS=never

WORKDIR /app

# Dependencies first, in their own layer. They change when the lockfile
# changes; the source changes constantly. Splitting the two means editing a
# module does not re-resolve the stack.
COPY pyproject.toml uv.lock README.md ./
RUN uv sync --frozen --no-dev --no-install-project

COPY src ./src
COPY panels ./panels
COPY examples ./examples
COPY app.py ./
RUN uv sync --frozen --no-dev

# ------------------------------------------------------------------- base ---
# Shared by both final stages. It is deliberately NOT last in this file: with
# no --target, `docker build` builds the LAST stage, so the CLI sits at the
# bottom and the app is opt-in. Getting that order wrong silently ships the
# heavier image, which is how the first version of this file came to be 2.1 GB.
FROM python:3.11-slim AS runtime-base

# HMMER is the reason this file exists.
RUN apt-get update \
    && apt-get install -y --no-install-recommends hmmer \
    && rm -rf /var/lib/apt/lists/* \
    && hmmscan -h > /dev/null

ENV PATH="/opt/venv/bin:$PATH" \
    PYTHONDONTWRITEBYTECODE=1

COPY --from=builder /opt/venv /opt/venv
COPY --from=builder /app /app
WORKDIR /app

# Fail the build rather than ship an image whose numbering is broken. This is
# the check that HMMER, ANARCI and the committed germline data agree inside the
# image, not merely on the machine that wrote this file.
RUN abscore score --quiet --no-epitopes --format summary \
    --seq EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS \
    | grep -q 'chain=H'

# -------------------------------------------------------------------- app ---
# The base plus Streamlit. Its extra pulls pyarrow, pydeck, altair and pillow --
# roughly 250 MB of packages that no score depends on -- so anyone scoring a
# sequence does not carry them.
FROM runtime-base AS app

COPY --from=ghcr.io/astral-sh/uv:0.12.5 /uv /usr/local/bin/uv
ENV UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_DOWNLOADS=never
RUN uv sync --frozen --no-dev --extra app && rm /usr/local/bin/uv

EXPOSE 8501
ENTRYPOINT ["streamlit", "run", "app.py", "--server.address", "0.0.0.0"]
CMD []

# ------------------------------------------------------------------- test ---
# Not part of any shipped image: building this target runs the whole suite
# inside the container and fails if anything is red. It is how the claim "the
# scores are identical on Linux" is checked rather than assumed, and it needs
# the dev group and the test data, neither of which belongs in the CLI image.
FROM runtime-base AS test

COPY --from=ghcr.io/astral-sh/uv:0.12.5 /uv /usr/local/bin/uv
ENV UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_DOWNLOADS=never
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --all-extras
COPY tests ./tests
RUN pytest -q

# -------------------------------------------------------------------- cli ---
# Last, so it is what a bare `docker build` produces.
FROM runtime-base AS cli

ENTRYPOINT ["abscore"]
CMD ["--help"]

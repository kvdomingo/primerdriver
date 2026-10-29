# PrimerDriver: Automated design of mutagenic PCR primers

![PrimerDriver](https://res.cloudinary.com/kdphotography-assets/image/upload/v1587460290/primerdriver/PrimerDriver_logo.png)

![GitHub](https://img.shields.io/github/license/kvdomingo/primerdriver)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/kvdomingo/primerdriver?include_prereleases)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fkvdomingo%2Fprimerdriver%2Fmain%2Fpyproject.toml)

## Introduction

_PrimerDriver_ is a user-friendly bioinformatics platform specifically designed
to generate primers for site-directed mutagenesis experiments into workflows
with an easy-to-use command-line interface. “Traversing” PrimerDriver lanes,
users can generate possible mutagenic primers upon input of a DNA sequence
(`DNA`), design primers for an array of species codon expression systems through
direct mutation of the amino acid (`PRO`), characterize and report user-designed
primers (`CHAR`).

> [!WARNING]
> `PRO` mode back-translates the input protein into the selected organism's
> preferred codons and designs primers against that reconstruction. The primers
> match a real template only if the gene was synthesized with exactly those
> codons; for a cloned or cDNA-derived gene, use `DNA` mode with the actual
> coding sequence. See the
> [documentation](https://primerdriver-docs.kvd.studio/software/#choosing-primer-design-modes).

## Usage

### A. Running standalone

You can access and download the CLI from the
[releases page](https://github.com/kvdomingo/primerdriver/releases). Currently,
we only have prebuilt binaries for 64-bit Windows and Linux. For other
OS/architectures, see the section on Building from Source, under Developing
Locally.

Run the program in a terminal using

```shell
primerdriver -h
```

This will run the help program. For first-time users, the program can be run in
interactive mode by passing the `-i` flag:

```shell
primerdriver -i
```

This will walk you through each option step-by-step.
Batch design can be performed by including
[`primerdriver`](primerdriver/__main__.py) as part of a shell script.

### B. Running as Python module

#### Prerequisites

- [Git](https://git-scm.com/downloads)
- [Mise](https://mise.jdx.dev)

#### Setup

1. Clone the repo to your local device and `cd` into it

    ```shell
    git clone https://github.com/kvdomingo/primerdriver.git
    cd primerdriver
    ```

2. Install dependencies:

    ```shell
    mise install
    uv sync
    ```

Run the program as follows:

```shell
uv run python -m primerdriver -h
```

### C. Running from Docker

Run as follows:

```shell
docker run -it --entrypoint python kvdomingo/primerdriver -m primerdriver -h
```

## Web application

For a more interactive experience, visit the
[web application](https://primerdriver.kvd.studio).

## Documentation

The documentation is available at <https://primerdriver-docs.kvd.studio>.

## Contributing

Open a PR or raise an
[issue](https://github.com/kvdomingo/primerdriver/issues).
You may also email Nomer or Kenneth, depending on the nature of the issue.

## Developing locally

### Prerequisites

- [Docker](https://www.docker.com/get-started)
- [Mise](https://mise.jdx.dev)

### Installing

A step by step series of examples that tell you how to get a
development environment running

1. Clone the repo.
2. Install prerequisites:

    ```shell
    mise install
    prek install
    ```

3. Install backend dependencies:

    ```shell
    uv sync
    ```

4. Run the development servers:

    ```shell
    task
    ```

    This runs Docker Compose in the foreground with `--watch`; containers
    restart when `uv.lock` or `ui/bun.lockb` change. Stop with `Ctrl+C` or
    `task stop`.

Wait a few minutes for all the containers to start, then access the
local servers in your browser at:

- Web app: <http://localhost:3000>
- API: <http://localhost:5000> (docs at <http://localhost:5000/api/docs>)
- Docs: <http://localhost:8000>

### Web app

The UI in [`ui/`](ui) is React + TypeScript, built with Vite and Bun, using
Zustand for shared state, TanStack Router for routing and TanStack Form with zod
schemas for the design forms. API calls go through openapi-react-query, typed by
`ui/src/api/openapi.gen.ts`; regenerate it after changing the API:

```shell
docker compose exec frontend bun run gen:api
```

The `frontend` container installs dependencies into the mounted `ui/`, so add
packages through it to avoid mixing file ownership in `ui/node_modules`:

```shell
docker compose exec frontend bun add <package>
```

### Tests and checks

```shell
uv run pytest            # backend tests
uv run ruff check .      # backend lint
cd ui
bun run build            # typecheck (tsc) and production build
bunx biome check src     # frontend lint and format
```

`prek` runs ruff, Biome and the other pre-commit hooks on each commit.

### Building from source

Build the CLI binary (PyInstaller, Linux x64) or the Docker images:

```shell
task build-cli    # standalone binary
task build-web    # web app + API image
task build-docs   # documentation image
```

### Deployment

```shell
git add .
git commit -m "DESCRIPTIVE_COMMIT_MESSAGE"
git push origin your_feature_branch
```

where `your_feature_branch` should summarize the changes you are implementing
following the
[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) format
(e.g., `feat/xxxx`, `fix/yyyy`).

## Authors

- **Numeriano Amer "Nomer" E. Gutierrez** - Project Lead, Molecular Biologist -
  [Email](mailto:ngutierrez@evc.pshs.edu.ph) |
  [GitHub](https://github.com/nomgutierrez)
- **Kenneth V. Domingo** - Technical Lead -
  [Email](mailto:kvdomingo@alum.up.edu.ph) | [Website](https://kvd.studio) |
  [GitHub](https://github.com/kvdomingo)
- **Shebna Rose D. Fabilloren** - Technical Consultant -
  [Email](mailto:sdfabilloren@up.edu.ph)
- **Carlo M. Lapid** - Project Adviser - [Email](mailto:cmlapid@up.edu.ph)

## Versioning

This project complies with [SemVer](https://semver.org) for versioning. For
all available versions, see
[tags](https://github.com/kvdomingo/primerdriver/tags).

## License

This project is licensed under the [GPLv3 License](./LICENSE).

Each pull request gets a preview at `https://primerdriver-pr-<n>.lab.kvd.studio`.

# PrimerDriver: Automated design of mutagenic PCR primers

![PrimerDriver](https://res.cloudinary.com/kdphotography-assets/image/upload/v1587460290/primerdriver/PrimerDriver_logo.png)

![GitHub](https://img.shields.io/github/license/kvdomingo/primerdriver)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/kvdomingo/primerdriver?include_prereleases)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/primerdriver)

> [!IMPORTANT]
> PrimerDriver is no longer under active development.
> It may occasionally receive minor updates, but these will not add new functionality.
> The web app will be kept running on a best-effort basis, but note that we do not guarantee any SLAs.
> This project remains open-source, so anyone may fork and deploy their own instance.

## Introduction

_PrimerDriver_ is a user-friendly bioinformatics platform specifically designed to generate primers for site-directed
mutagenesis experiments into workflows with an easy-to-use command-line interface. “Traversing” PrimerDriver lanes,
users can generate possible mutagenic primers upon input of a DNA sequence (`DNA`), design primers for an array of
species codon expression systems through direct mutation of the amino acid (`PRO`), characterize and report
user-designed primers (`CHAR`).

## Usage

### A. Running standalone

You can access and download the CLI from the
[releases page](https://github.com/kvdomingo/primerdriver/releases). Currently, we only have prebuilt binaries for
64-bit Windows and Linux. For other OS/architectures, see the section on Building from Source, under Developing Locally.

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
    ```

Run the program as follows:

```shell
poetry run python -m primerdriver -h
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

The documentation is available at https://primerdriver-docs.kvd.studio.

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
    pip install -U pre-commit
    pre-commit install
    ```
3. Install backend dependencies:
    ```shell
    poetry install --no-root --with dev
    ```
4. Run the development servers:
    ```shell
    task
    ```

Wait a few minutes for all the containers to start, then access the
local servers in your browser at:

- Web app: http://localhost:3000
- Docs: http://localhost:8000

### Building from source

Run the script:

```shell
# On a Linux machine, will build for Linux x64 only
# On a Windows machine, will build for Windows and Linux x64
task build
```

### Deployment

```shell
git add .
git commit -m "DESCRIPTIVE_COMMIT_MESSAGE"
git push origin your_feature_branch
```

where `your_feature_branch` should summarize the changes you are implementing following
the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) format
(e.g., `feat/xxxx`, `fix/yyyy`).

## Authors

- **Numeriano Amer "Nomer" E. Gutierrez** - Project Lead, Molecular
  Biologist - [Email](mailto:ngutierrez@evc.pshs.edu.ph) | [GitHub](https://github.com/nomgutierrez)
- **Kenneth V. Domingo** - Technical
  Lead - [Email](mailto:kvdomingo@alum.up.edu.ph) | [Website](https://kvd.studio) | [GitHub](https://github.com/kvdomingo)
- **Shebna Rose D. Fabilloren** - Technical Consultant - [Email](mailto:sdfabilloren@up.edu.ph)
- **Carlo M. Lapid** - Project Adviser - [Email](mailto:cmlapid@up.edu.ph)

## Versioning

This project complies with [SemVer](https://semver.org) for versioning. For
all available versions, see
[tags](https://github.com/kvdomingo/primerdriver/tags).

## License

This project is licensed under the [GPLv3 License](./LICENSE).

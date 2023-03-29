# PrimerDriver: Automated design of mutagenic PCR primers

![PrimerDriver](https://res.cloudinary.com/kdphotography-assets/image/upload/v1587460290/primerdriver/PrimerDriver_logo.png)

![GitHub](https://img.shields.io/github/license/kvdomingo/primerdriver)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/kvdomingo/primerdriver?include_prereleases)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/django)

## CLI

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

- [Python 3.10+](https://www.python.org/downloads/)
- [Poetry](https://python-poetry.org/docs/#installation)
- [Git](https://git-scm.com/downloads)

#### Setup

1. Clone the repo to your local device
    ```shell
    git clone https://github.com/kvdomingo/primerdriver.git
    ```

2. Install dependencies:
    ```shell
    poetry install
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
[web application](https://primerdriver.kvdstudio.app).

## Documentation

The documentation is available at https://kvdomingo.github.io/primerdriver/.

## Contributing

Open a PR or raise an
[issue](https://github.com/kvdomingo/primerdriver/issues).
You may also email Nomer or Kenneth, depending on the nature of the issue.

## Developing locally

### Prerequisites

- [Docker](https://www.docker.com/get-started)
- [Task](https://taskfile.dev/#/installation)

### Installing

A step by step series of examples that tell you how to get a
development environment running

1. Clone the repo.
2. Install backend dependencies:

```shell
poetry install
```

### Running local server

Setup the Docker containers:

```shell
task
```

Wait a few minutes for all the containers to start, then access the
local server in your browser at http://localhost:8000.

### Building from source

Run the script:

```shell
# Build for Windows on a Windows machine / for Linux on a Linux machine 
task build-all
```

### Deployment

```shell
git add .
git commit -m "DESCRIPTIVE_COMMIT_MESSAGE"
git push origin your_feature_branch
```

where `your_feature_branch` should summarize the changes you are implementing
(e.g., `feature/implementing-xxxx-feature`, `bugfix/crush-critical-yyyy-bug`).

## Authors

- **Numeriano Amer "Nomer" E. Gutierrez** - Project Lead, Molecular
  Biologist - [Email](mailto:ngutierrez@evc.pshs.edu.ph) | [GitHub](https://github.com/nomgutierrez)
- **Kenneth V. Domingo** - Lead Developer, Technical
  Consultant - [Email](mailto:kvdomingo@alum.up.edu.ph) | [Website](https://kvdomingo.xyz) | [GitHub](https://github.com/kvdomingo)
- **Shebna Rose D. Fabilloren** - Technical Consultant - [Email](mailto:sdfabilloren@up.edu.ph)
- **Carlo M. Lapid** - Project Adviser - [Email](mailto:cmlapid@up.edu.ph)

## Versioning

This project complies with [SemVer](https://semver.org) for versioning. For
all available versions, see
[tags](https://github.com/kvdomingo/primerdriver/tags).

## License

This project is licensed under the [GPLv3 License](./LICENSE).

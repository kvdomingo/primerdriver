# PrimerDriver: Automated design of mutagenic PCR primers
![PrimerDriver](https://raw.githubusercontent.com/kvdomingo/primerdriver-api/master/sdm/static/sdm/media/private/PrimerDriver_logo.png)

![GitHub](https://img.shields.io/github/license/kvdomingo/primerdriver)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/kvdomingo/primerdriver?include_prereleases)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/django)

## CLI

You can access and download the CLI from the 
[releases page](https://github.com/kvdomingo/primerdriver/releases). 
Extract the files to your local machine and run the program via
```shell
python -m primerdriver
```

This will automatically run the help program, which can also be accessed 
by passing a `-h` flag. For first-time users, the program can be run in 
interactive mode by passing the `-i` flag. This will walk you through 
each option step-by-step. Batch design can be performed by 
including [`primerdriver.py`](pdcli/primerdriver.py) as part of a shell script.

## Web application
For a more interactive experience, the updated web application can be 
accessed via https://primerdriver.vercel.app.

## Documentation
The documentation is available at https://kvdomingo.github.io/primerdriver/.

## Contributing
Open a PR or raise an 
[issue](https://github.com/kvdomingo/primerdriver/issues). 
You may also email Nomer or Kenneth, depending on the nature of the issue.

## Developing locally

### Prerequisites
- Python 3.9 or above
- Node.js LTS 14 or above
- Docker

### Installing
A step by step series of examples that tell you how to get a 
development environment running
   
1. Clone and extract the repo.
2. Install backend dependencies:
```shell
pip install -r requirements.dev.txt
```

### Running local server

Setup the Docker containers:
```shell
docker compose up --build
```

Wait a few minutes for all the containers to start, then access the 
local server in your browser at http://localhost:8000.

### Deployment
```bash
> git add .
> git commit -m <DESCRIPTIVE_COMMIT_MESSAGE>
> git push origin <GITHUB_BRANCH>
```

where `GITHUB_BRANCH` should summarize the changes you are implementing 
(e.g., `feature/implementing-xxxx-feature`, `bugfix/crush-critical-yyyy-bug`).


## Authors
- **Numeriano Amer "Nomer" E. Gutierrez** - Project Lead, Molecular Biologist - [Email](mailto:ngutierrez@evc.pshs.edu.ph) | [GitHub](https://github.com/nomgutierrez)
- **Kenneth V. Domingo** - Lead Developer, Technical Consultant - [Email](mailto:kvdomingo@up.edu.ph) | [Website](https://kvdomingo.xyz) | [GitHub](https://github.com/kvdomingo)
- **Shebna Rose D. Fabilloren** - Technical Consultant - [Email](mailto:sdfabilloren@up.edu.ph)
- **Carlo M. Lapid** - Project Adviser - [Email](mailto:cmlapid@up.edu.ph)

## Versioning
This project complies with [SemVer](https://semver.org) for versioning. Fo
all available versions, see 
[tags](https://github.com/kvdomingo/primerdriver/tags).

## License
This project is licensed under the [GPLv3 License](./LICENSE).

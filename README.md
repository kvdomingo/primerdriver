# PrimerDriver: Automated design of mutagenic PCR primers
![PrimerDriver](https://raw.githubusercontent.com/kvdomingo/primerdriver/master/sdm/static/sdm/media/private/PrimerDriver_logo.png)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT?style=for-the-badge)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/kvdomingo/primerdriver?include_prereleases)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/django)
![node-lts](https://img.shields.io/node/v-lts/react)

## CLI

### Prerequisites
- Python 3.6.8 and above
- Node.js LTS 10 and above

You can access and download the CLI from the [releases page](https://github.com/kvdomingo/primerdriver/releases). Extract the files to your local machine and run the program via
```bash
> pipenv run python primerdriver.py
```

This will automatically run the help program, which can also be accessed by passing a `-h` flag. For first-time users, the program can be run in interactive mode by passing the `-i` flag. This will walk you through each option step-by-step. Batch design can be performed by including [`primerdriver.py`](./primerdriver.py) as part of a shell script.

## Web application
For a more interactive experience, the updated web application can be accessed via https://primerdriver.herokuapp.com/.

## Documentation
The documentation is available at https://kvdomingo.github.io/primerdriver/.

## Contributing
Open a PR or raise an [issue](https://github.com/kvdomingo/primerdriver/issues). You can also email us.

## Developing locally

### Installing
A step by step series of examples that tell you how to get a development env running

1. Install Git, Python, and Node.js
1. Install/update Python package manager (`pipenv`)
```bash
> python -m pip install -U pipenv
```
1. Clone and extract repository to your machine
```bash
> git clone https://github.com/kvdomingo/primerdriver.git
```
1. Checkout a new `develop` branch. **DO NOT** make any modifications directly in the `master` branch. Similarly, **DO NOT** push directly to the `master` branch.
```bash
> git checkout -b develop
```
1. Install remaining prerequisites
```bash
> pipenv install --dev
> npm install -D
```

### Running local server

1. Open a terminal and run
```bash
> pipenv run python manage.py runserver
```
1. Open another terminal and run
```bash
> npm run watch
```

Access the local server in your browser at http://localhost:8000.

### Deployment
```bash
> git add .
> git commit -m <DESCRIPTIVE_COMMIT_MESSAGE>
> git push origin <GITHUB_BRANCH>
```

where `GITHUB_BRANCH` should summarize the changes you are implementing (e.g., `feature/implementing-xxxx-feature`, `bugfix/crush-critical-yyyy-bug`).


## Authors
- **Kenneth V. Domingo** - [Email](mailto:kvdomingo@up.edu.ph) | [Website](https://kvdomingo.xyz) | [GitHub](https://github.com/kvdomingo)
- **Numeriano Amer E. Gutierrez** - [Email](mailto:ngutierrez@evc.pshs.edu.ph) | [GitHub](https://github.com/nomgutierrez)
- **Shebna Rose D. Fabilloren** - [Email](mailto:sdfabilloren@up.edu.ph)
- **Carlo M. Lapid** - [Email](mailto:cmlapid@up.edu.ph)

## Versioning
This project complies with [SemVer](https://semver.org) for versioning. For all available versions, see [tags](https://github.com/kvdomingo/primerdriver/tags).

## License
This project is licensed under the [MIT License](./LICENSE).

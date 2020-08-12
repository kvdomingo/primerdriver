# PrimerDriver: Automated design of mutagenic PCR primers
![PrimerDriver](https://raw.githubusercontent.com/kvdomingo/primerdriver/master/sdm/static/sdm/media/private/PrimerDriver_logo.png)

[![CircleCI](https://circleci.com/gh/kvdomingo/primerdriver.svg?style=svg)](https://circleci.com/gh/kvdomingo/primerdriver?style=for-the-badge)
[![This project is using Percy.io for visual regression testing.](https://percy.io/static/images/percy-badge.svg)](https://percy.io/Kenneth-V-Domingo/primerdriver?style=for-the-badge)
[![codecov](https://codecov.io/gh/kvdomingo/primerdriver/branch/master/graph/badge.svg)](https://codecov.io/gh/kvdomingo/primerdriver?style=for-the-badge)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT?style=for-the-badge)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/kvdomingo/primerdriver?include_prereleases&style=for-the-badge)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/django?style=for-the-badge)

## CLI

### Prerequisites
- Python 3.6.8 and above
- Node.js LTS 10 and above

You can access and download the CLI from the [releases page](https://github.com/kvdomingo/primerdriver/releases). Extract the files to your local machine and run the program via

```shell
> python pdcli.py
```

This will automatically run the help program, which can also be accessed by passing a `-h` flag. For first-time users, the program can be run in interactive mode by passing the `-i` flag. This will walk you through each option step-by-step. Batch design can be performed by including [`pdcli.py`](./pdcli.py) as part of a shell script.

## Web application
For a more interactive experience, the updated web application can be accessed via https://primerdriver.herokuapp.com/.

## Documentation
The documentation is available at https://kvdomingo.github.io/primerdriver/.

## Contributing
Open a PR or raise an [issue](https://github.com/kvdomingo/primerdriver/issues). You can also email us.

## Authors
- **Kenneth V. Domingo** - [Email](mailto:kvdomingo@up.edu.ph) | [Website](https://kvdomingo.xyz) | [GitHub](https://github.com/kvdomingo)
- **Numeriano Amer E. Gutierrez** - [Email](mailto:ngutierrez@evc.pshs.edu.ph) | [GitHub](https://github.com/nomgutierrez)
- **Shebna Rose D. Fabilloren** - [Email](mailto:sdfabilloren@up.edu.ph)
- **Carlo M. Lapid** - [Email](mailto:cmlapid@up.edu.ph)

## Versioning
This project complies with [SemVer](https://semver.org) for versioning. For all available versions, see [tags](https://github.com/kvdomingo/primerdriver/tags).

## License
This project is licensed under the [MIT License](./LICENSE).

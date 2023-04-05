.. _getting_started:

Getting Started
===============

.. _requirements:

Requirements
------------

To use *PrimerDriver*, you should have:
    - `Python 3.8.1 or above <https://python.org/downloads>`_
    - `Poetry <https://python-poetry.org/docs/#installation>`_

Additional packages can be installed via the ``pyproject.toml`` file, which will be covered in the next section. `poetry` will automatically handle creating an isolated virtual environment for the project.


.. _installaton:

Installation
------------

First, download the latest release from `GitHub <https://github.com/kvdomingo/primerdriver/releases>`_ and extract it to your local machine. Open a terminal and ``cd`` to where you extracted the files. Install the required packages via

.. code-block:: console

    $ poetry install

Now you are ready to design mutagenic primers.


.. _quickstart:

Quickstart
----------

To make sure everything is in working order, run the following commands

.. code-block:: console

    $ poetry run python --version
    $ poetry run python -m primerdriver -h

If no errors occur, then you are good to go. The PrimerDriver CLI can cater to first-time users by providing an interactive mode which allows the user to go through the settings step-by-step. This is triggered by passing the ``-i`` flag when running the program with no additional arguments

.. code-block:: console

    $ poetry run python -m primerdriver -i

The user is greeted with a customized header and the first question::

                    ---.   .------------.
                    ||||\ /||||||||||||||\
                  Primer · Driver v1.3.2
            \|||||||||||/ \|||||||
             `---------`   `------

    (c) 2020 Kenneth Domingo & Nomer Gutierrez

    Enter primer mode [dna/pro/char]: _

Here, the user is asked for the desired mode, from DNA-based design (``dna``), protein-based design (``pro``), or characterization of a user-provided sequence (``char``). Try going through all the questions with some example answers as follows::

                    ---.   .------------.
                    ||||\ /||||||||||||||\
                  Primer · Driver v1.3.2
            \|||||||||||/ \|||||||
             `---------`   `------

    (c) 2020 Kenneth Domingo & Nomer Gutierrez

    Enter primer mode [dna/pro/char]: char
    Enter primer sequence: CGATCGTACGGACGCAGCTCGTAGCTACGATCGATCGATCGATCGTACGTACGTACGATCGTACGATCGATCGTACG
    Enter mutation type [s/i/d]: s
    Enter number of mismatched bases: 1

(*Note*: All inputs are case-insensitive). The resulting characterization is::

    |                  | Primer 1                                                                      |
    |------------------+-------------------------------------------------------------------------------|
    | Forward          | CGATCGTACGGACGCAGCTCGTAGCTACGATCGATCGATCGATCGTACGTACGTACGATCGTACGATCGATCGTACG |
    | Reverse          | CGTACGATCGATCGTACGATCGTACGTACGTACGATCGATCGATCGATCGTAGCTACGAGCTGCGTCCGTACGATCG |
    | Fwd length       | 77 bp                                                                         |
    | Rev length       | 77 bp                                                                         |
    | Fwd GC content   | 54.55%                                                                        |
    | Rev GC content   | 54.55%                                                                        |
    | Fwd melting temp | 94.76 C                                                                       |
    | Rev melting temp | 94.76 C                                                                       |
    | Fwd mol. weight  | 10083.12 g/mol                                                                |
    | Rev mol. weight  | 10074.11 g/mol                                                                |
    | Fwd mismatch     | 1.30%                                                                         |
    | Rev mismatch     | 1.30%                                                                         |
    | Fwd ends in G/C  | True                                                                          |
    | Rev ends in G/C  | True                                                                          |

    Save? [y/n] _

In the last line, the user is asked whether or not the characterization should be saved. If the answer is no (``n``), the program terminates. Otherwise, a final prompt appears::

    Enter filename: _

In this case, the program will automatically detect the file extension and save to that format. Supported file formats are ``.csv``, ``.txt``, ``.html``, and ``.fasta``. (*Note*: Primer characterizations cannot be saved as FASTA). Aside from manually entering a DNA sequence at the prompt, the user can instead provide the filename of a FASTA file containing the sequence, provided that the correct relative location is entered. By default, PrimerDriver will look in the same directory as ``primerdriver.py``.

Let's try another example. Suppose we have a FASTA file named ``read.fasta`` in the root directory. The content of the file is as follows::

    >Seq1
    CGATCGTACGGACGCAGCTCGTAGCTACGATCGATCGATCGATCGTACGTACGTACGATCGTACGATCGATCGTACG

Enter the following at the interactive CLI::


                    ---.   .------------.
                    ||||\ /||||||||||||||\
                  Primer · Driver v1.3.2
            \|||||||||||/ \|||||||
             `---------`   `------

    (c) 2020 Kenneth Domingo & Nomer Gutierrez

    Enter primer mode [dna/pro/char]: dna
    Enter primer sequence: read.fasta
    Enter mutation type [s/i/d]: s
    Enter target base: C
    Enter replacement for target base: G
    Enter position of target: 25

The output should yield::

    |                  | Primer 1                        |
    |------------------+---------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC     |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 27 bp                           |
    | Rev length       | 31 bp                           |
    | Fwd GC content   | 55.56%                          |
    | Rev GC content   | 58.06%                          |
    | Fwd melting temp | 76.05 C                         |
    | Rev melting temp | 77.28 C                         |
    | Fwd mol. weight  | 3554.21 g/mol                   |
    | Rev mol. weight  | 4102.70 g/mol                   |
    | Fwd mismatch     | 3.70%                           |
    | Rev mismatch     | 3.23%                           |
    | Fwd ends in G/C  | True                            |
    | Rev ends in G/C  | True                            |

    |                  | Primer 2                         |
    |------------------+----------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC      |
    | Reverse          | CGTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 27 bp                            |
    | Rev length       | 32 bp                            |
    | Fwd GC content   | 55.56%                           |
    | Rev GC content   | 59.38%                           |
    | Fwd melting temp | 76.05 C                          |
    | Rev melting temp | 77.69 C                          |
    | Fwd mol. weight  | 3554.21 g/mol                    |
    | Rev mol. weight  | 4213.80 g/mol                    |
    | Fwd mismatch     | 3.70%                            |
    | Rev mismatch     | 3.12%                            |
    | Fwd ends in G/C  | True                             |
    | Rev ends in G/C  | True                             |

    |                  | Primer 3                            |
    |------------------+-------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC         |
    | Reverse          | GATCGTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 27 bp                               |
    | Rev length       | 35 bp                               |
    | Fwd GC content   | 55.56%                              |
    | Rev GC content   | 57.14%                              |
    | Fwd melting temp | 76.05 C                             |
    | Rev melting temp | 77.87 C                             |
    | Fwd mol. weight  | 3554.21 g/mol                       |
    | Rev mol. weight  | 4626.18 g/mol                       |
    | Fwd mismatch     | 3.70%                               |
    | Rev mismatch     | 2.86%                               |
    | Fwd ends in G/C  | True                                |
    | Rev ends in G/C  | True                                |

    |                  | Primer 4                             |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC          |
    | Reverse          | CGATCGTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 27 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 55.56%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 76.05 C                              |
    | Rev melting temp | 78.28 C                              |
    | Fwd mol. weight  | 3554.21 g/mol                        |
    | Rev mol. weight  | 4737.28 g/mol                        |
    | Fwd mismatch     | 3.70%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 5                         |
    |------------------+----------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC      |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATCG |
    | Fwd length       | 27 bp                            |
    | Rev length       | 32 bp                            |
    | Fwd GC content   | 55.56%                           |
    | Rev GC content   | 59.38%                           |
    | Fwd melting temp | 76.05 C                          |
    | Rev melting temp | 77.69 C                          |
    | Fwd mol. weight  | 3554.21 g/mol                    |
    | Rev mol. weight  | 4253.83 g/mol                    |
    | Fwd mismatch     | 3.70%                            |
    | Rev mismatch     | 3.12%                            |
    | Fwd ends in G/C  | True                             |
    | Rev ends in G/C  | True                             |

    |                  | Primer 6                             |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC          |
    | Reverse          | GATCGTACGGACGCAGCTCGTAGGTACGATCGATCG |
    | Fwd length       | 27 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 55.56%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 76.05 C                              |
    | Rev melting temp | 78.28 C                              |
    | Fwd mol. weight  | 3554.21 g/mol                        |
    | Rev mol. weight  | 4777.31 g/mol                        |
    | Fwd mismatch     | 3.70%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 7                         |
    |------------------+----------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC      |
    | Reverse          | CGGACGCAGCTCGTAGGTACGATCGATCGATC |
    | Fwd length       | 27 bp                            |
    | Rev length       | 32 bp                            |
    | Fwd GC content   | 55.56%                           |
    | Rev GC content   | 59.38%                           |
    | Fwd melting temp | 76.05 C                          |
    | Rev melting temp | 77.69 C                          |
    | Fwd mol. weight  | 3554.21 g/mol                    |
    | Rev mol. weight  | 4213.80 g/mol                    |
    | Fwd mismatch     | 3.70%                            |
    | Rev mismatch     | 3.12%                            |
    | Fwd ends in G/C  | True                             |
    | Rev ends in G/C  | True                             |

    |                  | Primer 8                            |
    |------------------+-------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC         |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATCGATC |
    | Fwd length       | 27 bp                               |
    | Rev length       | 35 bp                               |
    | Fwd GC content   | 55.56%                              |
    | Rev GC content   | 57.14%                              |
    | Fwd melting temp | 76.05 C                             |
    | Rev melting temp | 77.87 C                             |
    | Fwd mol. weight  | 3554.21 g/mol                       |
    | Rev mol. weight  | 4626.18 g/mol                       |
    | Fwd mismatch     | 3.70%                               |
    | Rev mismatch     | 2.86%                               |
    | Fwd ends in G/C  | True                                |
    | Rev ends in G/C  | True                                |

    |                  | Primer 9                             |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC          |
    | Reverse          | CGTACGGACGCAGCTCGTAGGTACGATCGATCGATC |
    | Fwd length       | 27 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 55.56%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 76.05 C                              |
    | Rev melting temp | 78.28 C                              |
    | Fwd mol. weight  | 3554.21 g/mol                        |
    | Rev mol. weight  | 4737.28 g/mol                        |
    | Fwd mismatch     | 3.70%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 10                            |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC          |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATCGATCG |
    | Fwd length       | 27 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 55.56%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 76.05 C                              |
    | Rev melting temp | 78.28 C                              |
    | Fwd mol. weight  | 3554.21 g/mol                        |
    | Rev mol. weight  | 4777.31 g/mol                        |
    | Fwd mismatch     | 3.70%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 11                            |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATC          |
    | Reverse          | CGGACGCAGCTCGTAGGTACGATCGATCGATCGATC |
    | Fwd length       | 27 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 55.56%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 76.05 C                              |
    | Rev melting temp | 78.28 C                              |
    | Fwd mol. weight  | 3554.21 g/mol                        |
    | Rev mol. weight  | 4737.28 g/mol                        |
    | Fwd mismatch     | 3.70%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 12                       |
    |------------------+---------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG    |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 28 bp                           |
    | Rev length       | 31 bp                           |
    | Fwd GC content   | 57.14%                          |
    | Rev GC content   | 58.06%                          |
    | Fwd melting temp | 77.76 C                         |
    | Rev melting temp | 78.17 C                         |
    | Fwd mol. weight  | 3705.34 g/mol                   |
    | Rev mol. weight  | 4102.70 g/mol                   |
    | Fwd mismatch     | 3.57%                           |
    | Rev mismatch     | 3.23%                           |
    | Fwd ends in G/C  | True                            |
    | Rev ends in G/C  | True                            |

    |                  | Primer 13                        |
    |------------------+----------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG     |
    | Reverse          | CGTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 28 bp                            |
    | Rev length       | 32 bp                            |
    | Fwd GC content   | 57.14%                           |
    | Rev GC content   | 59.38%                           |
    | Fwd melting temp | 77.76 C                          |
    | Rev melting temp | 78.58 C                          |
    | Fwd mol. weight  | 3705.34 g/mol                    |
    | Rev mol. weight  | 4213.80 g/mol                    |
    | Fwd mismatch     | 3.57%                            |
    | Rev mismatch     | 3.12%                            |
    | Fwd ends in G/C  | True                             |
    | Rev ends in G/C  | True                             |

    |                  | Primer 14                           |
    |------------------+-------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG        |
    | Reverse          | GATCGTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 28 bp                               |
    | Rev length       | 35 bp                               |
    | Fwd GC content   | 57.14%                              |
    | Rev GC content   | 57.14%                              |
    | Fwd melting temp | 77.76 C                             |
    | Rev melting temp | 78.76 C                             |
    | Fwd mol. weight  | 3705.34 g/mol                       |
    | Rev mol. weight  | 4626.18 g/mol                       |
    | Fwd mismatch     | 3.57%                               |
    | Rev mismatch     | 2.86%                               |
    | Fwd ends in G/C  | True                                |
    | Rev ends in G/C  | True                                |

    |                  | Primer 15                            |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG         |
    | Reverse          | CGATCGTACGGACGCAGCTCGTAGGTACGATCGATC |
    | Fwd length       | 28 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 57.14%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 77.76 C                              |
    | Rev melting temp | 79.17 C                              |
    | Fwd mol. weight  | 3705.34 g/mol                        |
    | Rev mol. weight  | 4737.28 g/mol                        |
    | Fwd mismatch     | 3.57%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 16                        |
    |------------------+----------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG     |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATCG |
    | Fwd length       | 28 bp                            |
    | Rev length       | 32 bp                            |
    | Fwd GC content   | 57.14%                           |
    | Rev GC content   | 59.38%                           |
    | Fwd melting temp | 77.76 C                          |
    | Rev melting temp | 78.58 C                          |
    | Fwd mol. weight  | 3705.34 g/mol                    |
    | Rev mol. weight  | 4253.83 g/mol                    |
    | Fwd mismatch     | 3.57%                            |
    | Rev mismatch     | 3.12%                            |
    | Fwd ends in G/C  | True                             |
    | Rev ends in G/C  | True                             |

    |                  | Primer 17                            |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG         |
    | Reverse          | GATCGTACGGACGCAGCTCGTAGGTACGATCGATCG |
    | Fwd length       | 28 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 57.14%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 77.76 C                              |
    | Rev melting temp | 79.17 C                              |
    | Fwd mol. weight  | 3705.34 g/mol                        |
    | Rev mol. weight  | 4777.31 g/mol                        |
    | Fwd mismatch     | 3.57%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 18                           |
    |------------------+-------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG        |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATCGATC |
    | Fwd length       | 28 bp                               |
    | Rev length       | 35 bp                               |
    | Fwd GC content   | 57.14%                              |
    | Rev GC content   | 57.14%                              |
    | Fwd melting temp | 77.76 C                             |
    | Rev melting temp | 78.76 C                             |
    | Fwd mol. weight  | 3705.34 g/mol                       |
    | Rev mol. weight  | 4626.18 g/mol                       |
    | Fwd mismatch     | 3.57%                               |
    | Rev mismatch     | 2.86%                               |
    | Fwd ends in G/C  | True                                |
    | Rev ends in G/C  | True                                |

    |                  | Primer 19                            |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG         |
    | Reverse          | CGTACGGACGCAGCTCGTAGGTACGATCGATCGATC |
    | Fwd length       | 28 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 57.14%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 77.76 C                              |
    | Rev melting temp | 79.17 C                              |
    | Fwd mol. weight  | 3705.34 g/mol                        |
    | Rev mol. weight  | 4737.28 g/mol                        |
    | Fwd mismatch     | 3.57%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |

    |                  | Primer 20                            |
    |------------------+--------------------------------------|
    | Forward          | GCAGCTCGTAGGTACGATCGATCGATCG         |
    | Reverse          | GTACGGACGCAGCTCGTAGGTACGATCGATCGATCG |
    | Fwd length       | 28 bp                                |
    | Rev length       | 36 bp                                |
    | Fwd GC content   | 57.14%                               |
    | Rev GC content   | 58.33%                               |
    | Fwd melting temp | 77.76 C                              |
    | Rev melting temp | 79.17 C                              |
    | Fwd mol. weight  | 3705.34 g/mol                        |
    | Rev mol. weight  | 4777.31 g/mol                        |
    | Fwd mismatch     | 3.57%                                |
    | Rev mismatch     | 2.78%                                |
    | Fwd ends in G/C  | True                                 |
    | Rev ends in G/C  | True                                 |
    
    Too many results; truncating output...

    Save? [y/n] _

Notice that the output is silenced if the number of potential primers exceed 20. In this case, the exceeding primers are still stored and can be accessed by saving the results. The above example can be run in single command mode via

.. code-block:: console

    $ poetry run python -m primerdriver --mode dna --sequence read.fasta --mutation-type sub --target C --position 25 --replacement G --save primers.fasta

or in shorthand form via

.. code-block:: console

    $ poetry run python -m primerdriver -M dna -s read.fasta -m sub -t C -p 25 -r G --save primers.fasta

The ``--save`` argument is optional and can be omitted if the first 20 primers suffice for the user. As you can see, this can become a powerful tool especially when batch designing primers, by including it as part of a shell script.

Protein-based design works similarly and uses *Homo sapiens* expression system by default

.. code-block:: console

    $ poetry run python -m primerdriver --mode pro --sequence CAISBVAIVBAIVBCAICBASCBAVQVFEWQEPFQEHVSDBVSKZDBNCSD --mutation-type sub -t Q -p 26 -r R

.. _software_features:

*****************
Software Features
*****************


.. _parameter_setting:

Parameter Setting
=================

The user parameter setting is possible for the inputs in PrimerDriver. The tool designs the desired mutagenic primer pairs according to user input and performs site-directed mutagenesis from single- or multiple-bases within the sequence.

The user needs to open the ``primerdriver/settings.json`` file to define the range of allowable values in most of the parameters. The ``json`` file with the default settings appears as follows: ::
    
    {
        "Tm_range_min": 75,
        "Tm_range_max": 85,
        "gc_range_min": 40,
        "gc_range_max": 60,
        "length_min": 25,
        "length_max": 45,
        "flank5_range_min": 11,
        "flank5_range_max": 21,
        "flank3_range_min": 11,
        "flank3_range_max": 21,
        "forward_overlap5": 9,
        "forward_overlap3": 9,
        "terminate_gc" : 1,
        "center_mutation": 0,
        "primer_mode": "overlapping",
        "expression_system": "Homo sapiens"
    }

The following definitions and considerations should be made for designing mutagenic primers:

- ``Tm_range``
    This is where the user sets the minimum and maximum allowable melting temperature (:math:`T_m`) of primers in a PCR reaction. Depending on the type of mutation, the melting temperature is calculated differently. Estimation of the Tm of primers for base substitution is defined by

    .. math::

        T_m = 81.5 + 0.41(\%\textrm{GC}) - \frac{675}{N} - \%\textrm{mismatch}

    where :math:`N` is the primer length in base pairs, and :math:`\%\textrm{GC}, \%\textrm{mismatch}` are whole numbers. A modified formula is used for calculating :math:`T_m` that intends to introduce deletions or insertions:

    .. math::

        T_m = 81.5 + 0.41(\%\textrm{GC}) - \frac{675}{N}

- ``gc_range``
    The GC content is calculated as
    
    .. math::

        \%\textrm{GC} = \frac{\textrm{# of G + # of C}}{\textrm{total length of primer}}

    and the range of GC in a primer is to be set by the user. The user can follow the percentage as prescribed by a protocol or the user can just set the range from 0 to 100.

- ``length``
    The total length of the desired mutagenic primer

- ``flank5_range``
    This is the length of the region that is to the left of the mutation when read in standard sequence notation.

- ``flank3_range``
    This refers to the length of the region that is to the right of the mutation when read in 5’ to 3’.

- ``forward_overlap5``
    When generating overlapping primers, this is where the user sets the minimum number of overlapping bases of the 5’ end of the forward primer to the reverse primer.

- ``forward_overlap3``
    This sets the minimum number of base overlaps from the 3’ end of the forward primer to the reverse primer.

- ``terminate_gc``
    Set in boolean values, this refers to having a GC clamp at both ends of the primers.

- ``center_mutation``
    This sets the mutation is centered within the primer. Flanking regions differ by at most one base pair. These are set by boolean values.

- ``primer_mode``
    This gives the user the choice of having overlapping primers or complementary primers. Both modes use protocol from QuikChange™ Site-Directed Mutagenesis Kit by Stratagene®.

    - Complementary
        The values for the flanking regions will just be applied to the forward primer. The reverse primers will just be the reverse complements of the forward primers.
    
    - Overlapping
        Both the forward and reverse primers will be screened for the flanking region values. Forward overlap values will also be applied when generating the primer pairs.

- ``expression_system``
    When using Protein based primer design, the user can choose from a number of commonly used model organisms used in molecular and cell biology. This enables the tool to generate primers that use the optimum codon for the chosen organism.


.. _input_format:

Input format
============

In PrimerDriver, sequences can be added either by copying/pasting into the input cursor. The tool also accepts one or multiple sequences when uploaded in FASTA format. The input sequence must match the DNA sequence that will serve as template for PCR. Before processing the sequences, a basic input emending (FASTA headers, unsupported characters) is performed.


.. _interface:

Choosing an interface
=====================

As a powerful and useful feature in DNA engineering, the Web server offers the possibility of designing primers for site-directed mutagenesis. The tool can accommodate both DNA and protein sequences to incorporate base pair insertions, deletions, and substitutions as specified by the user. This can cater to an array of primer designs for site-directed single mutagenesis. PrimerDriver also lets you choose from two different command-line modes:

- Interactive mode
    Guides the user through a step-by-step input prompt. Primarily aimed towards demonstration purposes and first-time users of the program.

- Single-command mode
    Requires the user to include all the inputs and arguments needed to perform the task in a single command. Primarily aimed towards more advanced tasks such as batch primer design.


.. _design_mode:

Choosing primer design modes
============================

PrimerDriver asks you to choose from three different modes when using the command-line interface:

    - Primer characterization

    - DNA-based primer design

    - Protein-based primer design

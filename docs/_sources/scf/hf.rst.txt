.. _hf:

*************
Hartree--Fock
*************

===========
Description
===========

SCF can be run by specifying the following values to the keyword ``title``:

* Restricted HF for closed-shell systems: ``rhf`` or ``hf``
* Restricted open-shell HF: ``rohf``
* Unrestricted HF: ``uhf``
* Two-component HF using ECP basis sets spin-orbit SCF: ``soscf``

RHF can be run with the fast multipole method (FMM). For RHF-FMM, ``"cfmm" : "true"``
has to be specified in :ref:`molecule-toc`.

========
Keywords
========
The default values are recommended unless mentioned otherwise.

.. topic:: ``thresh``

   | **Description**: SCF convergence threshold for the root-mean-square of the error vector.
   | **Datatype**: double
   | **Default**: :math:`1.0\times 10^{-8}`

.. topic:: ``maxiter``/``maxiter_scf``

   | **Description**: number of iterations and number of SCF interations, after which the program will terminate if convergence is not reached.
   | **Datatype**: int
   | **Default**: :math:`100`

.. topic:: ``diis_start``

   | **Description**: after the specified iteration, we will begin using the DIIS algorithm to accelerate the convergence.
   | **Datatype**: int
   | **Default**: :math:`1`


.. topic:: ``thresh_overlap``

   | **Description**: Overlap threshold used to identify linear dependancies in the atomic basis set.
   | **Datatype**: double
   | **Default**: :math:`1.0\times 10^{-8}`

.. topic:: ``multipole``

   | **Description**: rank of Cartesian multipole moments to be printed out.
   | **Datatype**: int
   | **Values** : :math:`1, 2`
   | **Default** : :math:`1` (dipoles)

.. topic:: ``dma``

   | **Description**: options to print out multipole moments from distributed multipole analysis.
   | **Datatype**: int
   | **Values** : :math:`0, 1, 2, 3`
   | **Default** : :math:`0` (not print out)


.. topic:: ``charge``

   | **Description**: molecular charge
   | **Datatype**: int
   | **Default** : :math:`0`

.. topic:: ``nopen``

   | **Description**: number of unpaired electrons in high-spin unrestricted Hartree--Fock
   | **Datatype** : int
   | **Default** : (number of electrons - charge) % 2

.. topic:: ``restart``

   | **Description**: save an archive in each iteration to allow for restarting the calculation
   | **Datatype**: bool
   | **Default**: false

Keywords for RHF-FMM
====================

.. topic:: ``ns``

   | **Description**: level of descritization which controls the number of lowest-level boxes in one dimension for FMM
   | **Datatype**: int
   | **Default**: :math:`4`

.. topic:: ``ws``

   | **Description**: well-separatedness index, which is the number of boxes that must separate
                      two collections of charges before they are considered distant
                      and can interact through multipole expansions
   | **Datatype**: int
   | **Default**: :math:`0`

.. topic:: ``lmax``

   | **Description**: order of the multipole expansions in FMM-J
   | **Datatype**: int
   | **Default**: :math:`10`

.. topic:: ``exchange``

   | **Description**: whether to include far-field exchange using occ-RI-FMM
   | **Datatype**: bool 
   | **Default**: true 

.. topic:: ``lmax_exchange``

   | **Description**: order of the multipole expansions in FMM-K
   | **Datatype**: int
   | **Default**: :math:`2`

========
Examples
========
Below are some examples for SCF calculations using RHF, ROHF, UHF, SOSCF, and RHF-FMM.

RHF
===

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : "false",
     "geometry" : [
       { "atom" : "F",  "xyz" : [ -0.000000,     -0.000000,      2.720616]},
       { "atom" : "H",  "xyz" : [ -0.000000,     -0.000000,      0.305956]}
     ]
   },

   {
     "title" : "hf",
     "thresh" : 1.0e-8
   }

   ]}

The converged SCF energy is :math:`-99.84772354` after :math:`11` iterations.

ROHF
====
.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : "false",
     "geometry" : [
       { "atom" : "C",  "xyz" : [   -0.000000,     -0.000000,      3.000000] },
       { "atom" : "H",  "xyz" : [    0.000000,      0.000000,      0.000000] }
     ]
   },

   {
     "title" : "rohf",
     "nopen" : 1
   }

   ]}

The converged SCF energy is :math:`-38.16810629` after :math:`10` iterations.

UHF
===
.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : false,
     "geometry" : [
       { "atom" : "O",  "xyz" : [  -0.000000,     -0.000000,      1.500000]},
       { "atom" : "H",  "xyz" : [  -0.000000,     -0.000000,      0.000000]}
     ]
   },

   {
     "title" : "uhf",
     "nopen" : 1
   }

   ]}

The converged SCF energy is :math:`-75.28410147` after :math:`12` iterations. The expectation value of :math:`S^2` is :math:`0.7536`.

SOSCF
=====

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "ecp28mdf",
     "df_basis" : "tzvpp-jkfit",
     "angstrom" : "true",
     "geometry" : [
       { "atom" : "Br",  "xyz" : [  0.000000,      0.000000,      0.000000]},
       { "atom" :  "H",  "xyz" : [  0.000000,      1.420000,      0.000000],
                        "basis" : "sto-3g"}
     ]
   },

   {
     "title" : "soscf"
   }

   ]}

RHF-FMM
=======
.. figure:: hf-graphene.png
    :width: 30 % 
    :align: center
    :alt: alternate text
    :figclass: align-center


.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "3-21g",
     "angstrom" : "true",
     "cfmm" : "true",
     "schwarz_thresh" : "1.0e-8",
     "geometry" : [
       { "atom" : "C", "xyz" : [     -0.710000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      1.420000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -0.710000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -1.420000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.810000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.100000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.810000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.810000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.100000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.810000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.100000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.100000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.100000000,   -4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -3.550000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -2.840000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [      3.550000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -4.970000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -5.680000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [      4.970000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [     -3.550000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -2.840000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      3.550000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      2.840000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -4.970000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -5.680000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      4.970000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      5.680000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -3.550000000,    6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [     -2.840000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [      3.550000000,   -6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [      2.840000000,   -4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -4.970000000,    6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [     -5.680000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [      4.970000000,   -6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [      5.680000000,   -4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -2.840000000,    7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [     -2.840000000,   -7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      1.420000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -0.710000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -1.420000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -0.710000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -1.420000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      1.420000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,    6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [      1.420000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -0.710000000,   -6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [     -1.420000000,   -4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -0.710000000,    6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [     -1.420000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,   -6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [      1.420000000,   -4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,    8.608292514,    0.000000000] },
       { "atom" : "C", "xyz" : [      1.420000000,    7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [     -0.710000000,   -8.608292514,    0.000000000] },
       { "atom" : "C", "xyz" : [     -1.420000000,   -7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [     -0.710000000,    8.608292514,    0.000000000] },
       { "atom" : "C", "xyz" : [     -1.420000000,    7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [      0.710000000,   -8.608292514,    0.000000000] },
       { "atom" : "C", "xyz" : [      1.420000000,   -7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [      4.970000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      5.680000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [     -4.970000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      3.550000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      2.840000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [     -3.550000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      4.970000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      5.680000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -4.970000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -5.680000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      3.550000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      2.840000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -3.550000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -2.840000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      4.970000000,    6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [      5.680000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -4.970000000,   -6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [     -5.680000000,   -4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [      3.550000000,    6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [      2.840000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [     -3.550000000,   -6.148780367,    0.000000000] },
       { "atom" : "C", "xyz" : [     -2.840000000,   -4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [      2.840000000,    7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [      2.840000000,   -7.378536440,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.810000000,    1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.100000000,    0.000000000,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.810000000,   -1.229756073,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.810000000,    3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.100000000,    2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.810000000,   -3.689268220,    0.000000000] },
       { "atom" : "C", "xyz" : [     -7.100000000,   -2.459512147,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.100000000,    4.919024293,    0.000000000] },
       { "atom" : "C", "xyz" : [      7.100000000,   -4.919024293,    0.000000000] },
       { "atom" : "H", "xyz" : [      1.250000000,    9.543599950,    0.000000000] },
       { "atom" : "H", "xyz" : [     -1.250000000,   -9.543599950,    0.000000000] },
       { "atom" : "H", "xyz" : [      5.510000000,    7.084087803,    0.000000000] },
       { "atom" : "H", "xyz" : [     -5.510000000,   -7.084087803,    0.000000000] },
       { "atom" : "H", "xyz" : [      3.380000000,    8.313843876,    0.000000000] },
       { "atom" : "H", "xyz" : [      3.380000000,   -8.313843876,    0.000000000] },
       { "atom" : "H", "xyz" : [      7.640000000,    5.854331730,    0.000000000] },
       { "atom" : "H", "xyz" : [      7.640000000,   -5.854331730,    0.000000000] },
       { "atom" : "H", "xyz" : [     -7.640000000,    5.854331730,    0.000000000] },
       { "atom" : "H", "xyz" : [     -7.640000000,   -5.854331730,    0.000000000] },
       { "atom" : "H", "xyz" : [     -5.510000000,    7.084087803,    0.000000000] },
       { "atom" : "H", "xyz" : [      5.510000000,   -7.084087803,    0.000000000] },
       { "atom" : "H", "xyz" : [     -3.380000000,    8.313843876,    0.000000000] },
       { "atom" : "H", "xyz" : [     -3.380000000,   -8.313843876,    0.000000000] },
       { "atom" : "H", "xyz" : [     -1.250000000,    9.543599950,    0.000000000] },
       { "atom" : "H", "xyz" : [      1.250000000,   -9.543599950,    0.000000000] },
       { "atom" : "H", "xyz" : [      8.890000000,    1.229756073,    0.000000000] },
       { "atom" : "H", "xyz" : [     -8.890000000,   -1.229756073,    0.000000000] },
       { "atom" : "H", "xyz" : [      8.890000000,    3.689268220,    0.000000000] },
       { "atom" : "H", "xyz" : [     -8.890000000,   -3.689268220,    0.000000000] },
       { "atom" : "H", "xyz" : [     -8.890000000,    1.229756073,    0.000000000] },
       { "atom" : "H", "xyz" : [      8.890000000,   -1.229756073,    0.000000000] },
       { "atom" : "H", "xyz" : [     -8.890000000,    3.689268220,    0.000000000] },
       { "atom" : "H", "xyz" : [      8.890000000,   -3.689268220,    0.000000000] }
     ]
   },

   {
     "df" : false,
     "ns" : "5",
     "lmax" : "10",
     "ws" : "0.0",
     "exchange" : true,
     "lmax_exchange" : "2",
     "title" : "hf",
     "thresh" : 1.0e-6
   }

   ]}

The converged SCF energy is :math:`-3629.48403990` after :math:`11` iterations, and the far-field
exchange contribution is :math:`7.736` milli-hartree.

==========
References
==========

BAGEL References
================

+-----------------------------------------------+----------------------------------------------------------------------------+
|          Description of Reference             |                               Reference                                    |
+===============================================+============================================================================+
| Exact exchange evaluation using occ-RI-FMM    | H\.-A. Le and T. Shiozaki, J. Chem. Theory Comput. **14**, 1228 (2018)     |
+-----------------------------------------------+----------------------------------------------------------------------------+

General References
==================
+-----------------------------------------------+----------------------------------------------------------------------------------+
|          Description of Reference             |                               Reference                                          |
+===============================================+==================================================================================+
| General text on electronic structure theory   | A\. Szabo and N. S. Ostlund,                                                     |
|                                               | *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory* |
|                                               | (McGraw-Hill, New York, 1989).                                                   |
+-----------------------------------------------+----------------------------------------------------------------------------------+
| References for fast multipole method in       | C\. A. White, B. G. Johnson, P. M. W. Gill, and M. Head-Gordon,                  |
| quantum chemistry                             | Chem. Phys. Lett. **230**, 8 (1994).                                             |
+-----------------------------------------------+----------------------------------------------------------------------------------+
|                                               | M\. C. Strain, G. E. Scuseria, and M. J. Frisch, Science **271**, 51 (1996).     |
+-----------------------------------------------+----------------------------------------------------------------------------------+


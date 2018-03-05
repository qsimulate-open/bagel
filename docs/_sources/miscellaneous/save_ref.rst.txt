.. _save_ref:

******************************************
Saving and loading MOs and CI coefficients
******************************************

===========
Description
===========
A binary file can be produced to store the geometry, converged molecular orbitals, and for multi-configurational methods, CI coefficients.
This file can be read back to provide a reference for a future calculation.

For some non-relativistic methods, this feature is redundant with the restart capability provided by writing and reading Molden files, but
the binary archives used in this module are compatible with every method in BAGEL.

Binary archives generated using a different version of BAGEL might not be readable.

Commands: ``save_ref`` and ``load_ref``

========
Keywords
========

.. topic:: ``file``

   | **Description:** Name of the file to be generated.  The extension .archive will be appended.
   |      A complete path can be supplied, or the name of a file in the current working directory.
   | **Datatype:** string
   | **Default:**  "reference"

.. topic:: ``continue_geom``

   | **Description:**  This option is only used with load_ref.
   |      By default, the geometry is read from the binary archive file and used without modification.
   |      If this option is set to false, the molecular orbitals from the archive will be projected to a new geometry or
   |      basis set provided earlier in the input file.
   | **Datatype:** bool
   | **Default:**  true

=======
Example
=======

The following input files show how to use this restart feature.
The example given is for relativistic CASSCF calculations of the ground triplet for the oxygen dimer.
The initial calculation obtains a wavefunction for a 4-orbital active space comprised of the :math:`\pi` and :math:`\pi^\ast`
orbitals and saves the molecular orbital coefficients to a binary archive.

The second file loads this archive and uses it to give the initial guess for a calculation with a larger 6-orbital active space
which also includes the :math:`\sigma` and :math:`\sigma^\ast` orbitals.

The third file loads the same archive, projects to a slightly different geometry, and runs CASSCF at the new geometry.


Sample input
------------

Write binary archive:

.. code-block:: javascript

  { "bagel" : [

    {
      "title" : "molecule",
      "basis" : "tzvpp",
      "df_basis" : "tzvpp-jkfit",
      "angstrom" : true,
      "geometry" : [
        { "atom" : "O",  "xyz" : [   -0.000000,     -0.000000,      1.210000]},
        { "atom" : "O",  "xyz" : [   -0.000000,     -0.000000,      0.000000]}
      ]
    },

    {
      "title" : "hf",
      "charge" : "+2"
    },

    {
      "title" : "dhf",
      "charge" : "+2"
    },

    {
      "title"  : "zcasscf",
      "algorithm" : "second",
      "charge" : "0",
      "state" : [0, 0, 1],
      "thresh" : 1.0e-7,
      "nclosed"  : 5,
      "nact"   : 4
    },

    {
      "title" : "save_ref",
      "file" : "/your/directory/filename"
    }

  ]}


Load a binary archive without changing the reference:

.. code-block:: javascript

  { "bagel" : [

    {
      "title" : "load_ref",
      "file" : "/your/directory/filename"
    },

    {
      "title"  : "zcasscf",
      "algorithm" : "second",
      "charge" : "0",
      "state" : [0, 0, 1],
      "thresh" : 1.0e-7,
      "nclosed"  : 3,
      "nact"   : 8
    }

  ]}


Load a binary archive and use its orbitals to generate an initial guess at a new geometry.
We could also have changed the basis set in the input, rather than changing the atomic coordinates.

.. code-block:: javascript

  { "bagel" : [

    {
      "title" : "molecule",
      "basis" : "tzvpp",
      "df_basis" : "tzvpp-jkfit",
      "angstrom" : true,
      "geometry" : [
        { "atom" : "O",  "xyz" : [   -0.000000,     -0.000000,      1.220000]},
        { "atom" : "O",  "xyz" : [   -0.000000,     -0.000000,      0.000000]}
      ]
    },

    {
      "title" : "load_ref",
      "file" : "/your/directory/filename",
      "continue_geom" : false
    },

    {
      "title"  : "zcasscf",
      "algorithm" : "second",
      "charge" : "0",
      "state" : [0, 0, 1],
      "thresh" : 1.0e-7,
      "nclosed"  : 5,
      "nact"   : 4
    }

  ]}


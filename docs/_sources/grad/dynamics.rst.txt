.. _dynamics:

***************************
Interface to dynamics codes
***************************

Description
===========
BAGEL can be used for performing QM/MM calculations and direct dynamics simulations, by interfacing it to external software.
The "electrostatic embedding" scheme is used in inserting the QM region in MM region. BAGEL exports the gradients of the energy with respect
to the displacements of the QM particles and MM particles in a form that can be directly read by the dynamics codes.

How It Works
------------

The communications between the quantum chemistry codes and the molecular dynamics codes are based on the text files.
For example, in molecular dynamics (MD) jobs, after the trajectory initiates, the positions of the QM atoms and MM particles,
and the charges of MM particles are passed into BAGEL, in a form of text input.
The quantum chemical calculation is performed to calculate any needed quantities
(such as energies, its gradients, and NACMEs for nonadiabatic dynamics). This is passed (exported) to the form that can be
read by the dynamics program, and the dynamics program propagates the trajectory according to these gradients. After a
one time integration step, the positions of QM atoms and MM particles are again passed into BAGEL.
This continues until the trajectory integration ends.

In summary, to perform the dynamics simulations with external programs one needs:

(1) quantum chemistry input writer in the dynamics program

(2) quantum chemistry information output writer in the quantum chemistry program

(3) quantum chemistry output parser in the dynamics program

(2) is included in BAGEL (with a keyword of ``export`` and ``export_single``); the users should insert (1) and (3) in the dynamics code.

Input Format
------------

The BAGEL input for QM/MM jobs is same to that used in the quantum chemistry, except for the external charges.
The external MM charge can be added in BAGEL calculation, as

.. code-block:: javascript

  { "atom" : "Q", "xyz" : [   %lf,   %lf,   %lf ], "charge" :    %lf }

where ``xyz`` has the positions, and ``charge`` has the charge on the MM particle (in the units of e).


Export Format
-------------

Some of dynamics software reads the informations in a fixed form. The gradient export format is as follows.
The name of the file exported is ``ENERGY.txt`` (energy), ``FORCE_%d.txt`` (gradient, %d is the number of the state)
and ``NACME_%d_%d.txt`` (derivative coupling). The gradient is written in ``FORCE.txt`` when ``export_single`` is ``true``.

.. code-block:: cpp

                           0                   1                   2
        %6d         %20.10lf            %20.10lf            %20.10lf
  //(atomno)   (X component)       (Y component)       (Z component)

Note that the QM and MM atoms is not distinguished to each other in the exported output.
The exported gradients and nonadiabatic couplings are in the units of Hartrees/bohr and 1/bohr.
The energy is also exported, as

.. code-block:: cpp

          %20.10lf
    //   (statewise energy)

The energies for the multiple states in ascending order are exported in the same file when one does multi-state
calculations, such as SA-CASSCF and XMS-CASPT2.


Files to be Modified in Dynamics Software
-----------------------------------------

Some source codes in the dynamics software should be modified. Many of the dynamics software in fact have QM/MM
interfaces to the other quantum chemistry software, and by modifying the existent interface, one can perform the
QM/MM dynamics or gas phase dynamics simulations. For instance, to perform QM/MM calculations in ``GROMACS``, one should modify

``src/gromacs/mdlib/qmmm.c``,

add a BAGEL interface code to

``src/gromacs/mdlib/qm_bagel.c``.

This should contain the code for writing appropriately formatted input for BAGEL.
Additionally, include BAGEL for QM package for QM/MM in the CMake installation option in

``CMakeLists.txt``

(as of version 5.1.4).

Keywords
========

.. topic:: ``export``

   | **Description:** This option will export the nuclear gradient to a text file.
   | **Datatype:** bool
   | **Values:**
   |    ``true``: Export file
   |    ``false``: Do not export file
   | **Default:** false
   | **Recommendation:** Use ``true``.

.. topic:: ``export_single``

   | **Description:** This option will export the nuclear gradient to a text file for a single state.
   | **Datatype:** bool
   | **Values:**
   |    ``true``: Export file
   |    ``false``: Do not export file
   | **Default:** false
   | **Recommendation:** Use ``true`` with single state dynamics / optimizations.

Example
=======

Sample input
------------
A sample input for HF molecule using CASSCF.
This input computes the nuclear gradient for states 0 and 1 as well as the derivative coupling vector between these two states:

.. code-block:: javascript

  { "bagel" : [

  {
    "title" : "molecule",
    "basis" : "svp",
    "df_basis" : "svp-jkfit",
    "angstrom" : false,
    "geometry" : [
      { "atom" : "H",  "xyz" : [   -0.000000,     -0.000000,      1.700000] },
      { "atom" : "F",  "xyz" : [   -0.000000,     -0.000000,      0.000000] }
    ]
  },

  {
    "title" : "forces",
    "grads" : [
      { "title" : "force", "target" : 0 },
      { "title" : "force", "target" : 1 },
      { "title" : "nacme", "target" : 0, "target2" : 1 }
    ],
    "export" : true,
    "method" : [ {
      "title" : "casscf",
      "nopen" : 0,
      "nact" : 2,
      "nclosed" : 4,
      "nstate" : 2
    } ]
  }

  ]}

Executing BAGEL with this input will generate the following text files:

ENERGY.out
----------

.. code-block:: none

  -99.9135619715
  -99.5339025461

FORCE_0.out
-----------

.. code-block:: none

                         0                   1                   2
     0       -0.0000000000       -0.0000000000       -0.0019342533
     1        0.0000000000        0.0000000000        0.0019342533

FORCE_1.out
-----------

.. code-block:: none

                         0                   1                   2
     0        0.0000000000        0.0000000000       -0.2535235791
     1       -0.0000000000       -0.0000000000        0.2535235791

NACME_0_1.out
-------------

.. code-block:: none

                         0                   1                   2
     0       -0.0355272749       -0.0991581135        0.0000000000
     1       -0.0285426596       -0.0796637587       -0.0000000000

References
==========

General References
------------------

+-----------------------------------------------+--------------------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                             |
+===============================================+================================================================================+
| Nonadiabatic dynamics (Surface hopping)       | M\. Barbatti, WIREs Comput. Mol. Sci. **1**, 620 (2011).                       |
+-----------------------------------------------+--------------------------------------------------------------------------------+
| Excited state QM/MM in biomolecules           | E\. Brunk and U. Rothlisburger, Chem. Rev. **115**, 6217 (2015).               |
+-----------------------------------------------+--------------------------------------------------------------------------------+


.. _dhf:

********************
Dirac--Hartree--Fock
********************

===========
Description
===========

The Dirac--Hartree--Fock method performs a self-consistent field orbital optimization and energy calculation within a four-component relativistic framework.
The Dirac--Coulomb, Dirac--Coulomb--Gaunt, or full Dirac--Coulomb--Breit Hamiltonian can be used.
Density fitting is used for the two-electron integrals, and 2-spinor basis functions are generated using restricted kinetic balance (RKB).
External magnetic fields can be applied, in which case the spinor basis functions are generated using restricted magnetic balance (RMB) instead.

**Dirac--Hartree--Fock (DHF) should not be run with an odd number of electrons** in the absence of an external magnetic field, due to the Kramers degeneracy.
For open-shell molecules, it is recommended to run relativistic complete active space self-consistent field (ZCASSCF).
Dirac HF can be used to generate guess orbitals by increasing the molecular charge to remove unpaired electrons.

Calculations using DHF can be done using the keyword ``"title" : "dhf"``.

========
Keywords
========

The default values are recommended unless mentioned otherwise.

.. topic:: ``gaunt``

   | **Description**: Turns on the Gaunt interaction in the Hamiltonian.
   | **Datatype**: bool
   | **Default**: false

.. topic:: ``breit``

   | **Description**: Turns on the full Breit interaction in the Hamiltonian. 
   | **Datatype**: bool
   | **Default**: value copied from "gaunt" (if gaunt is true, breit is true)
   | **Recommendation**: Usually the Breit contribution is not important for molecular properties.

.. topic:: ``robust``

   | **Description**:  Determines whether or not to explicitly symmetrize the exchange matrix for numerical stability.
   | **Datatype**: bool
   | **Default**: false

.. topic:: ``maxiter (or maxiter_scf)``

   | **Description**:  Maximum number of iterations, after which the program will terminate if convergence is not reached.
   | **Datatype**: int
   | **Default**: 100

.. topic:: ``conv_ignore``

   | **Description:**  If set to "true," BAGEL will continue running even if the maximum iterations is reached without convergence.
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``diis_start``

   | **Description**:  After the specified iteration, we will begin using the DIIS algorithm the accelerate the convergence. 
   | **Datatype**: int
   | **Default**: 1

.. topic:: ``thresh (or thresh_scf)``

   | **Description**:  Convergence threshold for the root mean square of the error vector.
   | **Datatype**: double
   | **Default**: 1.0e-8

.. topic:: ``thresh_overlap``

   | **Description**:  Overlap threshold used to identify linear dependancy in the atomic basis set.
   | **Datatype**: double
   | **Default**: 1.0e-8

.. topic:: ``charge``

   | **Description**:  Molecular charge.
   | **Datatype**: int
   | **Default**: 0

.. topic:: ``pop``

   | **Description**:  If set to true, population analysis of the molecular orbitals will be printed to a file named dhf.log.
   | **Datatype**: bool
   | **Default**: false

=======
Example
=======

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
    "title" : "hf"
  },

  {
    "title" : "dhf",
    "gaunt" : true,
    "breit" : true
  }

  ]}

The non-relativistic SCF calculation converges in 13 iterations to :math:`-99.84772354`, and the Dirac HF converges after 9 iterations
to :math:`-99.92755305`.

==========
References
==========

BAGEL references
================
+-----------------------------------------------+-------------------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                            |
+===============================================+===============================================================================+
| Density fitted Dirac--Hartree--Fock method    | M\. S. Kelley and T. Shiozaki, J. Chem. Phys. **138**, 204113 (2013).         |
+-----------------------------------------------+-------------------------------------------------------------------------------+
| GIAO extension                                | R\. D. Reynolds and T. Shiozaki, Phys. Chem. Chem. Phys. **17**, 14280 (2015).|
+-----------------------------------------------+-------------------------------------------------------------------------------+

General references
==================
+-----------------------------------------------+-----------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                    |
+===============================================+=======================================================================+
| General text on relativistic electronic       | M\. Reiher and A. Wolf, *Relativistic Quantum Chemistry* (Wiley-VCH,  |
| structure, including Dirac--Hartree--Fock.    | Weinheim, 2009).                                                      |
+-----------------------------------------------+-----------------------------------------------------------------------+


.. _casscf:

****************************************************
Complete active space self-consistent field (CASSCF)
****************************************************

===========
Description
===========

CASSCF simultaneously optimizes the orbital coefficients and the CI coefficients of all the possible determinants generated from the active space.
The state-averaging scheme can be used for calculating excited states. When reference wavefunctions are given as a input to a CASSCF calculation,
a Hartree--Fock calculation is performed by default prior to CASSCF. For the Hartree--Fock options, see the SCF section.

A two-step second-order algorithm is implemented in BAGEL. At each macroiteration,
the CI coefficient is optimized by FCI calculations and the orbitals are updated. Microiterations are used for iteratively solving augmented Hessian problems for orbital updates.

The CASSCF algorithm in BAGEL is very robust. If it fails to converge, it is highly likely that your active space is wrong, or your guess orbitals are wrong.

Title: ``casscf``

========
Keywords
========

.. topic:: ``nstate``

   | **Description:** Number of states included in the state averaging.
   | **Datatype:** int
   | **Default:** 1.

.. topic:: ``nact``

   | **Description:** Number of active orbitals.
   | **Datatype:** int
   | **Default:** 0.

.. topic:: ``nclosed``

   | **Description:** Number of closed orbitals.
   | **Datatype:** int
   | **Default:** Number of electrons / 2.

.. topic:: ``active``

   | **Description:** Specify active orbitals. Note that the orbital index starts from 1.
   | **Datatype:** vector<int>
   | **Default:** Nact / 2 orbitals lower and higher from the valence orbital.
   | **Example:**
   |    [36, 37, 39] : include 36th, 37th, and 39th orbitals.

.. topic:: ``charge``

   | **Description**:  Molecular charge.
   | **Datatype**: int
   | **Default**: 0

.. topic:: ``nspin``

   | **Description:** The number associated with the spin states: 0 for singlet, 1 for doublet, 2 for triplet, etc. 
   | **Datatype:** int
   | **Default:** 0

.. topic:: ``algorithm``

   | **Description:** Orbital optimization algorithm.
   | **Datatype:** string
   | **Values:**
   |    ``second``: second-order algorithm.
   |    ``noopt``: no orbital optimization.
   | **Default:** ``second``

.. topic:: ``fci_algorithm``

   | **Description:** FCI algorithm employed in each macroiteration.
   | **Datatype:** string
   | **Values:**
   |    ``knowles``, ``handy``, ``kh``: Knowles--Handy Algorithm.
   |    ``harrison``, ``zarrabian``, ``hz``: Harrison--Zarrabian Algorithm.
   |    ``parallel``, ``dist``: Parallel FCI algorithm.
   | **Default:** ``parallel`` (when the number of active orbital is larger than 9 and number of process is larger than 8), ``knowles`` (otherwise)

.. topic:: ``thresh``

   | **Description:** Convergence threshold in macroiteration.
   | **Datatype:** double precision
   | **Default:** 1.0e-8.

.. topic:: ``thresh_micro``

   | **Description:** Convergence threshold in microiteration.
   | **Datatype:** double precision
   | **Default:** 5.0e-6.

.. topic:: ``conv_ignore``

   | **Description:**  If set to "true," BAGEL will continue running even if the maximum iterations is reached without convergence.  Normally an error is thrown and the program terminates.  
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``restart_cas``

   | **Description:**  If set to "true", after each macroiteration the orbitals will be written to a binary archive with filename "casscf_<iter>.archive". 
         They can be read back in using the "load_ref" module.  
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``natocc``

   | **Description**: If set to "true," occupation numbers of natural orbitals within the active space will be printed to casscf.log after each macroiteration.
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``canonical``

   | **Description**: If set to "true," optimized orbitals will be transformed to semi-canonical orbitals.
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``maxiter``

   | **Description:** Maximum number of macroiteration.
   | **Datatype:** int
   | **Default:** 50.
   | **Recommendation:** Increase if convergence is not obtained.

.. topic:: ``maxiter_micro``

   | **Description:** Maximum number of microiteration.
   | **Datatype:** int
   | **Default:** 100.

.. topic:: ``maxiter_fci``

   | **Description**: Maximum number of iterations in CI coefficient optimization 
   | **Datatype**: int
   | **Default**: copied from ``maxiter``

=======
Example
=======
Two-state CASSCF calculation of benzene. The active space of (6e,6o), which comprises three :math:`\pi` and three :math:`\pi^*` orbitals, is used.

Sample input
============

.. code-block:: javascript

  { "bagel" : [

  {
    "title" : "molecule",
    "basis" : "svp",
    "df_basis" : "svp-jkfit",
    "geometry" : [
    { "atom" : "C", "xyz" : [     -0.079002,      2.543870,      0.000000 ] },
    { "atom" : "C", "xyz" : [      2.557470,      2.543870,      0.000000 ] },
    { "atom" : "C", "xyz" : [      3.875630,      4.826190,      0.000000 ] },
    { "atom" : "C", "xyz" : [      2.557250,      7.109950,     -0.002266 ] },
    { "atom" : "C", "xyz" : [     -0.078588,      7.109800,     -0.003171 ] },
    { "atom" : "C", "xyz" : [     -1.396870,      4.826620,     -0.001289 ] },
    { "atom" : "H", "xyz" : [     -1.117900,      0.744245,      0.000850 ] },
    { "atom" : "H", "xyz" : [      3.595900,      0.743875,      0.002485 ] },
    { "atom" : "H", "xyz" : [      5.953730,      4.826340,      0.001198 ] },
    { "atom" : "H", "xyz" : [      3.596980,      8.909240,     -0.002377 ] },
    { "atom" : "H", "xyz" : [     -1.118170,      8.909350,     -0.004972 ] },
    { "atom" : "H", "xyz" : [     -3.474820,      4.826960,     -0.001629 ] }
    ]
  },
  {
    "title" : "hf"
  },
  {
    "title" : "casscf",
    "nstate" : 2,
    "nact" : 6,
    "nclosed" : 18,
    "active" : [17, 20, 21, 22, 23, 30]
  }
  ]}

the out of which is shown below. Note that the specified active orbitals are printed in the output.

.. code-block:: javascript

  ---------------------------
      CASSCF calculation     
  ---------------------------

 
    ==== Active orbitals : ===== 
         Orbital 17
         Orbital 20
         Orbital 21
         Orbital 22
         Orbital 23
         Orbital 30
    ============================ 

    * nstate   :      2
    * nclosed  :     18
    * nact     :      6
    * nvirt    :     90
  === CASSCF iteration (svp) ===

    * Using the second-order algorithm

      0  0      -230.58939332     3.18e-03      0.06
      0  1      -230.39960500     0.00e+00      0.06

         res : 9.89e-02   lamb: 1.00e+00   eps : -3.06e-02   step: 3.40e-01    0.02
         res : 2.13e-02   lamb: 1.00e+00   eps : -3.21e-02   step: 3.56e-01    0.02
         res : 2.72e-03   lamb: 1.00e+00   eps : -3.23e-02   step: 3.48e-01    0.02
         res : 2.64e-04   lamb: 1.00e+00   eps : -3.23e-02   step: 3.49e-01    0.03
         res : 3.09e-05   lamb: 1.00e+00   eps : -3.23e-02   step: 3.49e-01    0.04

      1  0      -230.60443140     3.79e-04      0.37
      1  1      -230.42264578     0.00e+00      0.37

         res : 1.06e-02   lamb: 1.00e+00   eps : -3.62e-04   step: 3.57e-02    0.02
         res : 2.13e-03   lamb: 1.00e+00   eps : -3.81e-04   step: 3.76e-02    0.02
         res : 3.31e-04   lamb: 1.00e+00   eps : -3.82e-04   step: 3.69e-02    0.02
         res : 4.53e-05   lamb: 1.00e+00   eps : -3.82e-04   step: 3.69e-02    0.03
         res : 4.21e-06   lamb: 1.00e+00   eps : -3.82e-04   step: 3.69e-02    0.04
         res : 5.33e-07   lamb: 1.00e+00   eps : -3.82e-04   step: 3.69e-02    0.02

      2  0      -230.60501843     1.17e-05      0.36
      2  1      -230.42244692     0.00e+00      0.36

         res : 2.38e-04   lamb: 1.00e+00   eps : -1.11e-07   step: 2.86e-04    0.02
         res : 3.75e-05   lamb: 1.00e+00   eps : -1.21e-07   step: 3.32e-04    0.02
         res : 8.19e-06   lamb: 1.00e+00   eps : -1.21e-07   step: 3.41e-04    0.02
         res : 1.11e-06   lamb: 1.00e+00   eps : -1.21e-07   step: 3.41e-04    0.03
         res : 8.72e-08   lamb: 1.00e+00   eps : -1.21e-07   step: 3.41e-04    0.02
         res : 2.97e-08   lamb: 1.00e+00   eps : -1.21e-07   step: 3.41e-04    0.03

      3  0      -230.60502169     3.97e-07      0.36
      3  1      -230.42244379     0.00e+00      0.36

         res : 6.62e-06   lamb: 1.00e+00   eps : -1.91e-10   step: 1.92e-05    0.02
         res : 1.37e-06   lamb: 1.00e+00   eps : -1.98e-10   step: 1.91e-05    0.02
         res : 2.63e-07   lamb: 1.00e+00   eps : -1.99e-10   step: 1.82e-05    0.02
         res : 4.76e-08   lamb: 1.00e+00   eps : -1.99e-10   step: 1.82e-05    0.03
         res : 5.65e-09   lamb: 1.00e+00   eps : -1.99e-10   step: 1.82e-05    0.02
         res : 1.68e-09   lamb: 1.00e+00   eps : -1.99e-10   step: 1.82e-05    0.03

      4  0      -230.60502176     2.01e-08      0.36
      4  1      -230.42244372     0.00e+00      0.36

         res : 3.09e-07   lamb: 1.00e+00   eps : -6.58e-13   step: 1.33e-06    0.03
         res : 7.14e-08   lamb: 1.00e+00   eps : -6.74e-13   step: 1.33e-06    0.03
         res : 1.28e-08   lamb: 1.00e+00   eps : -6.76e-13   step: 1.28e-06    0.04
         res : 2.75e-09   lamb: 1.00e+00   eps : -6.76e-13   step: 1.28e-06    0.02

      5  0      -230.60502177     1.21e-09      0.33
      5  1      -230.42244372     0.00e+00      0.33

    * Second-order optimization converged. *   


==========
References
==========

The second-order orbital optimization is implemented with an assistance of Takeshi Yanai (Institute for Molecular Science, Japan).

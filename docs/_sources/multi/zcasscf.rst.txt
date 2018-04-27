.. _zcasscf:

********************************************************************
Relativistic complete active space self-consistent field (RelCASSCF)
********************************************************************

===========
Description
===========
The relativistic analogue of CASSCF is implemented in BAGEL. The second-order algorithm which is basically the same as its non-relativistic analogue is used by default.

Title: ``zcasscf``

========
Keywords
========
.. topic:: ``state``

   | **Description**: Number of states computed for each spin value.  All are included in the state-averaging procedure when orbitals are optimized.
   | **Datatype**: vector<int>
   | **Default**:  There is no default; this parameter must be supplied in the input.
   | **Note**:  An array of integers is supplied, where each one indicates the number of states for a given spin value.  For example,
          the input [ 1 ] gives a singlet ground state, while [ 3, 0, 1 ] gives three singlets and one triplet (6 states total).
          Be careful!  While the spin values you specified are used in generating guess CI coefficients, the spin sectors will mix, and the
          algorithm returns the *n* lowest eigenstates regardless of their spin expectation values.

.. topic:: ``nact``

   | **Description**: Number of active orbitals
   | **Datatype**: int
   | **Default**: 0

.. topic:: ``nclosed``

   | **Description**:  Number of closed orbitals
   | **Datatype**: int
   | **Default**: Number of electrons / 2. 

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

.. topic:: ``algorithm``

   | **Description:** Orbital optimization algorithm.
   | **Datatype:** string
   | **Values:**
   |    ``second``: second-order algorithm.
   |    ``noopt``: no orbital optimization.
   | **Default:** ``second``

.. topic:: ``gaunt``

   | **Description**: Turns on the Gaunt interaction in the Hamiltonian.
   | **Datatype**: bool
   | **Default**: false

.. topic:: ``breit``

   | **Description**: Turns on the full Breit interaction in the Hamiltonian. 
   | **Datatype**: bool
   | **Default**: value copied from "gaunt" (if gaunt is true, breit is true)
   | **Recommendation**: Usually the Breit contribution is not important for molecular properties.

.. topic:: ``only_electrons``

   | **Description**:  This option allows the user to freeze all positronic orbitals and optimize only for rotations between electronic orbitals.
   | **Datatype**: bool
   | **Default**:   false

.. topic:: ``natocc``

   | **Description**: If set to "true," occupation numbers of natural orbitals within the active space will be printed to casscf.log after each macroiteration.
   | **Datatype**: bool
   | **Default**: false

.. topic:: ``canonical``

   | **Description**: If set to "true," optimized orbitals will be transformed to semi-canonical orbitals.
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``hcore_guess``

   | **Description**:  If set to true, the one-electron Hamiltonian is diagonalized to generate initial guess orbitals.
   | **Datatype**: bool
   | **Default**: false

.. topic:: ``maxiter``

   | **Description**: Maximum number of macroiterations.
   | **Datatype**: int
   | **Default**: 100

.. topic:: ``maxiter_micro``

   | **Description**: Maximum number of microiterations.
   | **Datatype**: int
   | **Default**: 20 

.. topic:: ``maxiter_fci``

   | **Description**: Maximum number of iterations in CI coefficient optimization 
   | **Datatype**: int
   | **Default**: copied from ``maxiter``

.. topic:: ``thresh_fci``

   | **Description**: Convergence threshold for the CI coefficients
   | **Datatype**: double
   | **Default**: Value copied from ``thresh``

.. topic:: ``conv_ignore``

   | **Description:**  If set to "true," BAGEL will continue running even if the maximum iterations is reached without convergence.  Normally an error is thrown and the program terminates.  
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``restart_cas``

   | **Description:**  If set to "true", after each macroiteration the orbitals will be written to a binary archive with filename "zcasscf_<iter>.archive". 
         They can be read back in using the "load_ref" module.  
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``pop``

   | **Description**: If set to true, population analysis of the molecular orbitals will be printed to a file names dhf.log.
   | **Datatype**: bool
   | **Default**: false

.. topic:: ``davidson_subspace``

   | **Description**:  Number of vectors retained in the limited-memory Davidson algorithm.
   | **Datatype**: int
   | **Default**: 20
   | **Recommendation**: Reduce if an insufficient amount of memory is available (do not reduce to a value lower than 3). 

.. topic:: ``print_thresh``

   | **Description**: Threshold below which CI coefficients are not printed.  
   | **Datatype**: double
   | **Default**: 0.05

.. topic:: ``spin_adapt``

   | **Description**: Spin-adapt the starting guess. 
   | **Datatype**: bool
   | **Default**: true
   | **Recommendation**: Use false if the error "generate_guess produced an invalid determinant" is generated. 

.. topic:: ``aniso``

   | **Description**: Performs magnetic anisotropy analysis using a pseudospin Hamiltonian, to obtain g-factors and zero-field splitting parameters 
   | **Datatype**: input block with several paramaters, given below:
   |
   |
   |   ``nspin`` 
   |
   |   **Description**: 2 * the spin value of the pseudospin Hamiltonian
   |   **Datatype**: integer
   |   **Default**: Number of states in the CASSCF calculation minus one
   |
   |   ``ranks`` 
   |
   |   **Description**: Which ranks of the Extended Stevens Operators to include in the zero-field splitting Hamiltonian
   |   **Datatype**: Array of integers
   |   **Default**: Even integers from 2 to ``nspin``
   |   **Recommendation**: Normally use default.  The array [ 2, 4, 6 ] can be used to specify the commonly-used sixth-order Hamiltonian if using the giant spin approximation.  
   |   
   |   ``print_operators`` 
   |
   |   **Description**: Option to print additional information.  If true, the Extended Stevens Operators will be printed in matrix form, along with the magnetic moment, spin, and orbital angular momentum matrices.  
   |   **Datatype**: bool
   |   **Default**: false
   |   
   |   ``zaxis`` and ``xaxis``
   |
   |   **Description**: Can be used to specify the orientation of the axes along which we should quantize the spin when mapping the zero-field anisotropy
   |   **Datatype**: Arrays of 3 doubles
   |   **Default**: The primary g-anisotropy axes are used
   |   
   |   ``states`` 
   |
   |   **Description**: Which electronic states are to be mapped to the pseudospin Hamiltonian
   |   **Datatype**: Array of integers
   |   **Default**: The lowest-lying states are used
   |   
   |   ``center``
   |
   |   **Description**: Coordinates of the spin center (typically a metal atom for a mononuclear complex), given in Bohr radii
   |   **Datatype**: Arrays of 3 doubles
   |   **Default**: [  0.0,  0.0,  0.0 ]
   |   
   |   ``energies``
   |
   |   **Description**: Energies to be used in the anisotropy analysis.  This option can be used to replace the CASSCF energies with those generated using another method, such as CASPT2.  
   |   **Datatype**: Arrays of nspin+1 doubles
   |   **Default**: The relativistic CASSCF energies are used.  
   |   **Recommendation**: Use the default.  This option is deprecated now that the "aniso" module is interfaced directly to relativistic CASPT2 and MRCI.  
   |   
   |   ``diagop``
   |
   |   **Description**: Operator to be diagonalized in order to determine a mapping from ab initio electronic states to pseudospin states.  
   |   **Datatype**: string - can be "Mu" for magnetic moment, "J", "S", or "L"
   |   **Default**: Mu
   |   
   |   ``phases`` and ``phase_full``
   |
   |   **Description**: This is a debugging tool.  It allows the use to apply an arbitrary phase shift to the pseudospin states.  
   |   **Datatype**: Array of doubles; double
   |   **Default**: No phase shift is applied

=======
Example
=======

.. code-block:: javascript

   { "bagel" : [ 
   
   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : false,
     "geometry" : [ 
       { "atom" : "F",  "xyz" : [ 0.000000, 0.000000, 3.720616 ]},
       { "atom" : "H",  "xyz" : [ 0.000000, 0.000000, 0.305956 ]}
     ]
   },
   
   {
     "title"  : "zcasscf",
     "state" : [1],
     "thresh" : 5.0e-7,
     "nact"   : 2,
     "nclosed"  : 4 
   }
   
   ]}


from which one obtains

.. code-block:: javascript

  ---------------------------
      CASSCF calculation     
  ---------------------------

  *** Geometry (Relativistic) ***
       - 3-index ints post                         0.00
       - 3-index ints prep                         0.00
       - 3-index ints                              0.00
       - 3-index ints post                         0.00

       - Geometry relativistic (total)             0.01

    * nclosed  :      4
    * nact     :      2
    * nvirt    :     32
    * gaunt    : false
    * breit    : false
    * active space: 2 electrons in 2 orbitals
    * time-reversal symmetry will be assumed.
       - Coulomb: half trans                       0.01
       - Coulomb: metric multiply                  0.03
       - Coulomb: J operator                       0.00
       - Coulomb: K operator                       0.00
    * nstate   :      1

  === Dirac CASSCF iteration (svp) ===

   * Using the second-order algorithm

         0   0        -99.87309219     1.03e-02      0.07

         res : 1.24e-01   lamb: 1.00e+00   eps : -4.46e-02   step: 2.31e-01    0.06
         res : 2.29e-02   lamb: 1.00e+00   eps : -4.93e-02   step: 2.94e-01    0.06
         res : 2.90e-03   lamb: 1.00e+00   eps : -4.94e-02   step: 2.93e-01    0.05
         res : 8.08e-04   lamb: 1.00e+00   eps : -4.94e-02   step: 2.94e-01    0.05
         res : 1.96e-04   lamb: 1.00e+00   eps : -4.94e-02   step: 2.94e-01    0.05
         res : 3.17e-04   lamb: 1.00e+00   eps : -4.94e-02   step: 2.94e-01    0.05
         res : 1.13e-04   lamb: 1.00e+00   eps : -4.94e-02   step: 2.94e-01    0.06
         res : 1.91e-05   lamb: 1.00e+00   eps : -4.94e-02   step: 2.94e-01    0.06
         1   0        -99.90008454     9.20e-04      0.57

         res : 7.91e-03   lamb: 1.00e+00   eps : -2.52e-04   step: 1.27e-02    0.06
         res : 1.77e-03   lamb: 1.00e+00   eps : -2.72e-04   step: 1.61e-02    0.05
         res : 5.19e-04   lamb: 1.00e+00   eps : -2.75e-04   step: 1.68e-02    0.05
         res : 7.75e-04   lamb: 1.00e+00   eps : -2.75e-04   step: 1.75e-02    0.05
         res : 5.66e-04   lamb: 1.00e+00   eps : -2.76e-04   step: 1.84e-02    0.06
         res : 2.02e-04   lamb: 1.00e+00   eps : -2.76e-04   step: 1.89e-02    0.05
         res : 6.00e-05   lamb: 1.00e+00   eps : -2.76e-04   step: 1.90e-02    0.06
         res : 6.44e-06   lamb: 1.00e+00   eps : -2.76e-04   step: 1.90e-02    0.06
         res : 1.06e-06   lamb: 1.00e+00   eps : -2.76e-04   step: 1.90e-02    0.06
         2   0        -99.90024315     1.95e-04      0.62

         res : 1.42e-03   lamb: 1.00e+00   eps : -1.03e-05   step: 2.44e-03    0.06
         res : 3.04e-04   lamb: 1.00e+00   eps : -1.11e-05   step: 3.12e-03    0.05
         res : 9.86e-05   lamb: 1.00e+00   eps : -1.11e-05   step: 3.18e-03    0.05
         res : 4.00e-05   lamb: 1.00e+00   eps : -1.11e-05   step: 3.19e-03    0.05
         res : 4.83e-05   lamb: 1.00e+00   eps : -1.11e-05   step: 3.20e-03    0.06
         res : 3.89e-05   lamb: 1.00e+00   eps : -1.12e-05   step: 3.24e-03    0.06
         res : 1.11e-05   lamb: 1.00e+00   eps : -1.12e-05   step: 3.25e-03    0.06
         res : 2.28e-06   lamb: 1.00e+00   eps : -1.12e-05   step: 3.25e-03    0.05
         res : 4.11e-07   lamb: 1.00e+00   eps : -1.12e-05   step: 3.25e-03    0.05
         res : 9.17e-08   lamb: 1.00e+00   eps : -1.12e-05   step: 3.25e-03    0.05
         3   0        -99.90025026     5.44e-05      0.67

         res : 3.82e-04   lamb: 1.00e+00   eps : -7.76e-07   step: 6.56e-04    0.05
         res : 8.18e-05   lamb: 1.00e+00   eps : -8.31e-07   step: 8.43e-04    0.05
         res : 2.63e-05   lamb: 1.00e+00   eps : -8.36e-07   step: 8.59e-04    0.05
         res : 9.66e-06   lamb: 1.00e+00   eps : -8.36e-07   step: 8.61e-04    0.06
         res : 1.02e-05   lamb: 1.00e+00   eps : -8.36e-07   step: 8.61e-04    0.06
         res : 1.06e-05   lamb: 1.00e+00   eps : -8.37e-07   step: 8.71e-04    0.05
         res : 2.98e-06   lamb: 1.00e+00   eps : -8.37e-07   step: 8.73e-04    0.05
         res : 6.25e-07   lamb: 1.00e+00   eps : -8.37e-07   step: 8.73e-04    0.05
         res : 1.21e-07   lamb: 1.00e+00   eps : -8.37e-07   step: 8.73e-04    0.05
         4   0        -99.90025079     1.49e-05      0.60

         res : 1.04e-04   lamb: 1.00e+00   eps : -5.79e-08   step: 1.78e-04    0.05
         res : 2.24e-05   lamb: 1.00e+00   eps : -6.20e-08   step: 2.29e-04    0.05
         res : 7.19e-06   lamb: 1.00e+00   eps : -6.23e-08   step: 2.34e-04    0.05
         res : 2.64e-06   lamb: 1.00e+00   eps : -6.23e-08   step: 2.34e-04    0.06
         res : 2.77e-06   lamb: 1.00e+00   eps : -6.24e-08   step: 2.35e-04    0.05
         res : 2.90e-06   lamb: 1.00e+00   eps : -6.24e-08   step: 2.37e-04    0.05
         res : 8.21e-07   lamb: 1.00e+00   eps : -6.24e-08   step: 2.38e-04    0.05
         res : 1.72e-07   lamb: 1.00e+00   eps : -6.24e-08   step: 2.38e-04    0.05
         5   0        -99.90025083     4.07e-06      0.55

         res : 2.83e-05   lamb: 1.00e+00   eps : -4.30e-09   step: 4.86e-05    0.05
         res : 6.10e-06   lamb: 1.00e+00   eps : -4.61e-09   step: 6.25e-05    0.05
         res : 1.96e-06   lamb: 1.00e+00   eps : -4.63e-09   step: 6.37e-05    0.06
         res : 7.19e-07   lamb: 1.00e+00   eps : -4.64e-09   step: 6.39e-05    0.05
         res : 7.56e-07   lamb: 1.00e+00   eps : -4.64e-09   step: 6.39e-05    0.06
         res : 7.92e-07   lamb: 1.00e+00   eps : -4.64e-09   step: 6.46e-05    0.05
         res : 2.25e-07   lamb: 1.00e+00   eps : -4.64e-09   step: 6.48e-05    0.06
         6   0        -99.90025083     1.11e-06      0.51

         res : 7.70e-06   lamb: 1.00e+00   eps : -3.20e-10   step: 1.32e-05    0.06
         res : 1.66e-06   lamb: 1.00e+00   eps : -3.43e-10   step: 1.70e-05    0.05
         res : 5.35e-07   lamb: 1.00e+00   eps : -3.44e-10   step: 1.74e-05    0.06
         res : 2.00e-07   lamb: 1.00e+00   eps : -3.45e-10   step: 1.74e-05    0.05
         7   0        -99.90025083     3.03e-07      0.36

    * Second-order optimization converged. *   



==========
References
==========

BAGEL references
================
+-----------------------------------------------+-----------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                    |
+===============================================+=======================================================================+
| Relativistic CASSCF                           | J\. E. Bates and T. Shiozaki, J. Chem. Phys. **142**, 044112 (2015).  |
+-----------------------------------------------+-----------------------------------------------------------------------+
| Second-order Relativistic CASSCF              | R\. D. Reynolds, T. Yanai, and T. Shiozaki, arXiv:1804.06470 (2018).  |
+-----------------------------------------------+-----------------------------------------------------------------------+

General references
==================
+-----------------------------------------------+-----------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                    |
+===============================================+=======================================================================+
| General text on relativistic electronic       | M\. Reiher and A. Wolf, *Relativistic Quantum Chemistry* (Wiley-VCH,  |
| structure                                     | Weinheim, 2009).                                                      |
+-----------------------------------------------+-----------------------------------------------------------------------+

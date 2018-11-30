.. _force:

****************************************
Nuclear gradient and derivative coupling
****************************************

Description
===========
The force section can be used to compute the analytical gradient (force), the numerical gradient by finite difference, or the non-adiabatic coupling matrix elements (NACME). Analytical gradients are implemented for unrestricted Hartree–Fock (UHF), restricted open-shell Hartree–Fock (ROHF), restricted Hartree–Fock (RHF), Dirac–Hartree–Fock (DHF), Møller–Plesset perturbation theory (MP2), complete active space self consistent field (CASSCF), and multireference perturbation theory (CASPT2). For CASSCF and CASPT2, it is also possible to perform multiple evaluations of the gradients without repeating the energy calculations.
In addition, the relaxed density in MP2, SA-CASSCF, and CASPT2 calculations can be printed in the Gaussian Cube format.

Keywords
========

Required Keywords
-----------------

.. topic:: ``title``

   | **Description:** The title of the type of gradient calculation being performed.
   | **Datatype:** string
   | **Values:** (force, nacme, dgrad, forces)
   |    ``force``: Calculates the gradient (force)
   |    ``nacme``: Calculates the non-adiabatic coupling matrix elements
   |    ``dgrad``: Difference gradient (only available for CASSCF)
   |    ``forces``: Calculates many gradients [available for multi-state methods, such as SA-CASSCF and (X)MS-CASPT2]
   | **Default:** N/A

.. topic:: ``method``

   | **Description:** The method input block allows the user to specify one or more methods to be used in the gradient calculation. See section on input structure for more information.

.. topic:: ``grads``

   | **Description:** The grads input block allows the user to specify one or more gradients to be calculated. The keywords ``title``, ``target``, ``target2``, ``nacmtype``, ``maxziter`` can be specified for each gradient evaluation. See below for example.

.. topic:: ``nacmtype``

   | **Description:** Type of non-adiabatic coupling matrix element to be computed
   | **Datatype:** string
   | **Values:**
   |    ``full``: Full non-adiabatic coupling
   |    ``interstate``: Interstate coupling
   |    ``etf``: Non-adiabtic coupling with built-in electronic translational factor (ETF)
   |    ``noweight``: Interstate coupling without weighting it by energy gap
   | **Default:** full.

Optional Keywords
-----------------

.. topic:: ``numerical``

   | **Description:** The gradients will be computed by finite difference.
   | **Datatype:** bool
   | **Values:**
   |    ``true``: Uses finite difference
   |    ``false`` : Uses analytical gradient
   | **Default:** false
   | **Recommendation:** If available, use analytical gradient. If analytical gradient is not available, BAGEL automatically switches to numerical gradient.

.. topic:: ``dx``

   | **Description:** The step size used in the displacements in the finite difference calculations. The units are bohr.
   | **Datatype:** double precision
   | **Default:** 1.0e-3

.. topic:: ``target``

   | **Description:** The target state for the energy and gradient evaluation (e.g. which state in a state-averaged CASSCF calculation)
   | **Datatype:** int
   | **Values:**
   |    ``int``: ground state = 0
   | **Default:** 0

.. topic:: ``target2``

   | **Description:** In an NACME or DGRAD calculation, target2 designates the target state for the second state.
   | **Datatype:** int
   | **Values:**
   |    ``int``: first exited state = 1
   | **Default:** 1

.. topic:: ``dipole``

   | **Description:** Calculate unrelaxed (transition) dipole moments, available in SA-CASSCF and XMS-CASPT2.
   | **Datatype:** bool
   | **Default:** false

.. topic:: ``ciderivative``

   | **Description:** Evaluate the CI derivative and full gradient in XMS-CASPT2. If false, only unrelaxed density will be calculated.
   | **Datatype:** bool
   | **Default:** true

.. topic:: ``export``

   | **Description:** This option will export the nuclear gradient to a text file.
   | **Datatype:** bool
   | **Default:** false
   | **Recommendation:** This is used to interface with the QM/MM program. See section on non-adiabatic dynamics.

.. topic:: ``export_single``

   | **Description:** This option will export the nuclear gradient to a text file for a single state.
   | **Datatype:** bool
   | **Default:** false
   | **Recommendation:** This is used to interface with the QM/MM program. See section on non-adiabatic dynamics.

.. topic:: ``maxziter``

   | **Description:** Maximum number of Z-vector iterations for gradient evaluation. Applies to SA-CASSCF, CASPT2, and MP2 calculations.
   | **Datatype:** int
   | **Default:** 100
   | **Recommendation:** Increase the value when Z-vector equation does not converge.

.. topic:: ``nproc``

   | **Description:** The numerical gradient code is embarrassingly parallelized so that the displacements in the finite difference calculations can be run at the same time. The nproc keyword allows the user to specify the number of MPI processes to be used for each energy calculation.
   | **Datatype:** int
   | **Default:** 1

.. topic:: ``density_print``

   | **Description:** Print relaxed densities in the Gaussian Cube format. Applies to SA-CASSCF, CASPT2, and MP2 calculations. The options for density printing can be specified in ``moprint`` block (see below for example and :ref:`here <moprint>` for the keywords).
   | **Datatype:** bool
   | **Default:** false

Example
=======
The benzophenone molecule

.. figure:: benzophenone.png
    :width: 200px
    :align: center
    :alt: alternate text
    :figclass: align-center

    The benzophenone molecule with carbon atoms in grey, oxygen in red, and hydrogen in white.

Sample input: force
-------------------

.. code-block:: javascript

  { "bagel" : [

  {
    "title" : "molecule",
    "basis" : "cc-pvdz",
    "df_basis" : "cc-pvdz-jkfit",
    "angstrom" : false,
    "geometry" : [
    { "atom" : "C", "xyz" : [     -2.002493,     -2.027773,      0.004882 ] },
    { "atom" : "C", "xyz" : [     -2.506057,     -4.613700,      0.009896 ] },
    { "atom" : "C", "xyz" : [      0.536515,     -1.276360,      0.003515 ] },
    { "atom" : "C", "xyz" : [     -0.558724,     -6.375134,      0.013503 ] },
    { "atom" : "H", "xyz" : [     -4.396140,     -5.341490,      0.011057 ] },
    { "atom" : "C", "xyz" : [      2.478233,     -3.024614,      0.007049 ] },
    { "atom" : "H", "xyz" : [      0.959539,      0.714937,     -0.000292 ] },
    { "atom" : "C", "xyz" : [      1.936441,     -5.592475,      0.012127 ] },
    { "atom" : "H", "xyz" : [     -1.012481,     -8.367883,      0.017419 ] },
    { "atom" : "H", "xyz" : [      4.418042,     -2.380738,      0.005919 ] },
    { "atom" : "H", "xyz" : [      3.448750,     -6.968581,      0.014980 ] },
    { "atom" : "C", "xyz" : [     -6.758666,     -0.057378,      0.001157 ] },
    { "atom" : "C", "xyz" : [     -8.231109,     -2.241648,      0.000224 ] },
    { "atom" : "C", "xyz" : [     -8.022986,      2.269249,      0.001194 ] },
    { "atom" : "C", "xyz" : [    -10.853532,     -2.110536,     -0.000769 ] },
    { "atom" : "H", "xyz" : [     -7.410047,     -4.093049,      0.000224 ] },
    { "atom" : "C", "xyz" : [    -10.632155,      2.405932,      0.000369 ] },
    { "atom" : "H", "xyz" : [     -6.913797,      3.976253,      0.001805 ] },
    { "atom" : "C", "xyz" : [    -12.064741,      0.207004,     -0.000695 ] },
    { "atom" : "H", "xyz" : [    -11.941318,     -3.840822,     -0.001614 ] },
    { "atom" : "H", "xyz" : [    -11.548963,      4.232744,      0.000447 ] },
    { "atom" : "H", "xyz" : [    -14.107194,      0.302907,     -0.001460 ] },
    { "atom" : "C", "xyz" : [     -3.892311,      0.136360,      0.001267 ] },
    { "atom" : "O", "xyz" : [     -3.026383,      2.227189,     -0.001563 ] }
    ]
  },

  {
    "title" : "force",
     "method" : [ {
      "title" : "hf",
      "thresh" : 1.0e-12
    } ]
  }
 ]}


Using the same molecule block, a XMS-CASPT2 analytical gradient calculation can be performed.
In this particular example as is often the case, the active keyword is used to select the orbitals for the active space that includes 4 electrons and 3 orbitals.
Three sets of  :math:`\pi` and :math:`\pi^*` orbitals localized on the phenyl rings are included along with one non-bonding orbital (oxygen lone pair).
The CASSCF orbitals are state-averaged over two states.

.. code-block:: javascript

  {
    "title" : "casscf",
    "nstate" : 2,
    "nclosed" : 46,
    "nact" : 3,
    "active" : [37, 44, 49]
  },

  {
    "title" : "force",
    "target" : 0,
    "method" : [ {
      "title" : "caspt2",
      "smith" : {
        "method" : "caspt2",
        "ms" : "true",
        "xms" : "true",
        "sssr" : "true",
        "shift" : 0.2,
        "frozen" : true
      },
      "nstate" : 2,
      "nact" : 3,
      "nclosed" : 46
    } ]
  }

Sample input: NACME and DGRAD
-----------------------------

.. code-block:: javascript

  {
   "title" : "nacme",
   "target" : 0,
   "target2" : 1,
   "method" : [ {
     "title" : "caspt2",
     "smith" : {
       "method" : "caspt2",
       "ms" : "true",
       "xms" : "true",
       "sssr" : "true",
       "shift" : 0.2,
       "frozen" : true
     },
     "nstate" : 3,
     "nact" : 7,
     "nclosed" : 44
   } ]
  }

Using the keyword ``forces``, you can run multiple gradient or derivative coupling calculations without repeating the energy calculations. The example below evaluates the nuclear gradient of the energy of the ground state, the first excited state, and the interstate coupling vector (``nacmtype`` is ``interstate``) between these two states.

.. code-block:: javascript

  {
   "title" : "forces",
   "grads" : [
     { "title" : "force", "target" : 0 },
     { "title" : "force", "target" : 1 },
     { "title" : "nacme", "target" : 0, "target2" : 1, "nacmtype" : "interstate" }
   ],
   "method" : [ {
     "title" : "caspt2",
     "smith" : {
       "method" : "caspt2",
       "ms" : "true",
       "xms" : "true",
       "sssr" : "true",
       "shift" : 0.2,
       "frozen" : true
     },
     "nstate" : 3,
     "nact" : 7,
     "nclosed" : 44
   } ]
  }

Sample input: Printing dipole moments
-------------------------------------

You can use the keyword ``dipole`` to compute the dipole moments for all the states, the transition dipole moments and oscillator strengths of all possible pairs of electronic states.

.. code-block:: javascript

  {
   "title" : "forces",
   "dipole" : "true",
   "method" : [ {
     "title" : "caspt2",
     "smith" : {
       "method" : "caspt2",
       "ms" : "true",
       "xms" : "true",
       "sssr" : "true",
       "shift" : 0.2,
       "frozen" : true
     },
     "nstate" : 3,
     "nact" : 7,
     "nclosed" : 44
   } ]
  }

Alternatively, for CASPT2, you can use the keyword ``"ciderivative" : "false"`` to evaluate the dipole moments, transition dipole moments and oscillator strengths of selected states,
without computing the CI derivative term and gradient.
The example below evaluates the dipole moments of state 0, state 1, and transition dipole moment of the state 0 -- state 1 transition.

.. code-block:: javascript

  {
   "title" : "forces",
   "grads" : [
     { "title" : "force", "target" : 0, "ciderivative" : "false" },
     { "title" : "force", "target" : 1, "ciderivative" : "false" },
     { "title" : "nacme", "target" : 0, "target2" : 1, "nacmtype" : "interstate", "ciderivative" : "false" }
   ],
   "method" : [ {
     "title" : "caspt2",
     "smith" : {
       "method" : "caspt2",
       "ms" : "true",
       "xms" : "true",
       "sssr" : "true",
       "shift" : 0.2,
       "frozen" : true
     },
     "nstate" : 3,
     "nact" : 7,
     "nclosed" : 44
   } ]
  }


Sample input: Printing relaxed density
--------------------------------------

The relaxed density in MP2, SA-CASSCF, and CASPT2 gradient calculations can be printed in the Gaussian Cube format (``density.cub`` file) by setting the keyword ``density_print`` to true
(see the example below).

.. code-block:: javascript

  { "bagel" : [

  {
    "title" : "molecule",
    "basis" : "svp",
    "df_basis" : "svp-jkfit",
    "geometry" : [
      { "atom" : "Li", "xyz" : [ 0.000000, 0.000000, 6.000000] },
      { "atom" : "F",  "xyz" : [ 0.000000, 0.000000, 0.000000] }
    ]
  },

  {
    "title" : "force",
    "target" : 0,
    "density_print" : true,
    "moprint" : {
      "inc_size" : [ 0.20, 0.20, 0.20]
    },
    "method" : [ {
      "title" : "caspt2",
      "smith" : {
        "method" : "caspt2",
        "shift" : 0.2,
        "frozen" : true
      },
      "nstate" : 4,
      "nact" : 4,
      "nclosed" : 3
    } ]
  }

  ]}

Using the keyword ``forces``, you can print multiple relaxed densities without repeating the energy calculations.
The example below prints all the XMS-CASPT2 relaxed densities.


.. code-block:: javascript

  {
    "title" : "forces",
    "grads" : [
      { "title" : "force", "target" : 0, "density_print" : true, "moprint" : { "density_filename" : "density_0" } },
      { "title" : "force", "target" : 1, "density_print" : true, "moprint" : { "density_filename" : "density_1" } },
      { "title" : "force", "target" : 2, "density_print" : true, "moprint" : { "density_filename" : "density_2" } },
      { "title" : "force", "target" : 3, "density_print" : true, "moprint" : { "density_filename" : "density_3" } }
    ],
    "method" : [ {
      "title" : "caspt2",
      "smith" : {
        "method" : "caspt2",
        "ms" : "true",
        "xms" : "true",
        "sssr" : "true",
        "shift" : 0.2,
        "frozen" : true
      },
      "nstate" : 4,
      "nact" : 4,
      "nclosed" : 3
    } ]
  }


References
==========

BAGEL References
----------------
+-----------------------------------------------+---------------------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                              |
+===============================================+=================================================================================+
| SS-CASPT2 gradient                            | M\. K. MacLeod and T. Shiozaki, J. Chem. Phys. **142**, 051103 (2015).          |
+-----------------------------------------------+---------------------------------------------------------------------------------+
| (X)MS-CASPT2 gradient                         | B\. Vlaisavljevich and T. Shiozaki, J. Chem. Theory Comput. **12**, 3781 (2016).|
+-----------------------------------------------+---------------------------------------------------------------------------------+
| (X)MS-CASPT2 derivative coupling              | J\. W. Park and T. Shiozaki, J. Chem. Theory Comput. **13**, 2561 (2017).       |
+-----------------------------------------------+---------------------------------------------------------------------------------+
| Current (X)MS-CASPT2 gradient algorithm       | J\. W. Park and T. Shiozaki, J. Chem. Theory Comput. **13**, 3676 (2017).       |
+-----------------------------------------------+---------------------------------------------------------------------------------+

General References
------------------

+-----------------------------------------------+--------------------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                             |
+===============================================+================================================================================+
| General review of gradient methods            | P\. Pulay, WIREs Comput. Mol. Sci. **4**, 169-181 (2014).                      |
+-----------------------------------------------+--------------------------------------------------------------------------------+


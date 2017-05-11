.. _dimer_asd:

*********
Dimer ASD
*********


Description
===========
The active space decomposition algorithm for molecular dimers allows for efficient computation of the dimer's complete-active-space wave functions. The current algorithm works for dimer molecules whose fragments are not covalently linked. Full-CI and restricted-active-space-CI can be used to obtain the fragment state wave functions. ASD calculation starts with dimer molecule construction, see `Dimer Construction`_ section for more information.


Dimer Construction
==================
.. toctree::
   :maxdepth: 1

   dimer.rst


Keywords
========

Required Keywords
-----------------

.. topic:: ``method``

   | **Description:** Method to compute active subspaces.
   | **Datatype:** string
   | **Value:**
   |   ``cas`` or ``fci``: use full configuration interaction method
   |   ``ras``: use restricted active space configuration interaction method

.. topic:: ``fci``

   | **Description:** If ``fci`` is used, specify the implementations here.
   | **Recommendation:** See :ref:`fci` for details.

.. topic:: ``ras``

   | **Description:** If ``ras`` is used, specify the implementations here.
   | **Recommendation:** See :ref:`rasci` for details.

.. topic:: ``space``

   | **Description:** Specify important fragment states with the following keys:
   |  ``charge``, ``spin``, ``nstate``
   | **Recommendation:** See sample input for details.

Optional Keywords
-----------------

.. topic:: ``nstates``

   | **Description:** Number of target states.
   | **Datatype:** int
   | **Default** 10

.. topic:: ``charge``

   | **Description:** Dimer charge.
   | **Datatype:** int
   | **Default:** 0

.. topic:: ``nspin``

   | **Description:** Number of dimer total spin.
   | **Datatype:** int
   | **Default:** 0

.. topic:: ``nguess``

   | **Description:** Number of initial guess state vectors for Davidson diagonalization.
   | **Datatype:** int
   | **Default:** :math:`10\times nstates`

.. topic:: ``Davidson_subspace``

   | **Description:** Size of Davidson subspace.
   | **Datatype:** int
   | **Default:** 10

.. topic:: ``max_iter``

   | **Description:** Maximum number of iterations for Davidson diagonalization.
   | **Datatype:** int
   | **Default:** 50

.. topic:: ``dipoles``

   | **Description:** Whether to calculate dipole moment.
   | **Datatype:** bool
   | **Default:** false

.. topic:: ``thresh``

   | **Description:** Threshold for convergence in Davidson diagonalization.
   | **Datatype:** double
   | **Default:** :math:`1.0\times 10^{-7}`

.. topic:: ``print_thresh``

   | **Description:** Threshold for printing out important configurations.
   | **Datatype:** double
   | **Default:** 0.01

.. topic:: ``store matrix``

   | **Description:** Whether the Hamiltonian matrix is stored.
   | **Datatye:** bool
   | **Default:** false

.. topic:: ``print_info``

   | **Description:** Whether print out information (e.g. reduced density matrix and energy).
   | **Datatype:** bool
   | **Default:** false


Example
=======
Here is a sample calculation of a benzene dimer molecule.

.. figure:: ../figure/benzene_dimer.png
    :width: 100px
    :align: center
    :alt: alternate text
    :figclass: align-center

Sample input
------------

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : false,
     "cartesian" : false,
     "geometry" : [
       {"atom" :"C", "xyz" : [    0.00000000000,     0.00000000000,     2.64112304663] },
       {"atom" :"C", "xyz" : [    2.28770766388,     0.00000000000,     1.32067631141] },
       {"atom" :"C", "xyz" : [    2.28770047235,     0.00000000000,    -1.32071294538] },
       {"atom" :"C", "xyz" : [    0.00000000000,     0.00000000000,    -2.64114665444] },
       {"atom" :"C", "xyz" : [   -2.28770047235,     0.00000000000,    -1.32071294538] },
       {"atom" :"C", "xyz" : [   -2.28770766388,     0.00000000000,     1.32067631141] },
       {"atom" :"H", "xyz" : [    4.07221260176,     0.00000000000,     2.35164689765] },
       {"atom" :"H", "xyz" : [    4.07221517814,     0.00000000000,    -2.35163163881] },
       {"atom" :"H", "xyz" : [    0.00000000000,     0.00000000000,    -4.70191324441] },
       {"atom" :"H", "xyz" : [   -4.07221517814,     0.00000000000,    -2.35163163881] },
       {"atom" :"H", "xyz" : [   -4.07221260176,     0.00000000000,     2.35164689765] },
       {"atom" :"H", "xyz" : [    0.00000000000,     0.00000000000,     4.70197960246] }
     ]
   },

   {
     "title" : "hf"
   },

   {
     "title" : "dimerize",
     "angstrom" : true,
     "translate" : [0.0, 4.0, 0.0],
     "dimer_active" : [17, 20, 21, 22, 23, 24],
     "hf" : {
       "thresh" : 1.0e-12
     },
     "localization" : {
       "max_iter" : 50,
       "thresh" : 1.0e-8
     }
   },

   {
     "title" : "asd",
     "method" : "cas",
     "store_matrix" : false,
     "space" : [
       { "charge" : 0, "spin" : 0, "nstate" : 3},
       { "charge" : 0, "spin" : 2, "nstate" : 3},
       { "charge" : 1, "spin" : 1, "nstate" : 3},
       { "charge" :-1, "spin" : 1, "nstate" : 3}
     ],
     "fci" : {
       "thresh" : 1.0e-6,
       "algorithm" : "kh",
       "nguess" : 400
     },
     "nstates" : 5
   }

   ]}


Reference
=========
+-----------------------------------------------+--------------------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                             |
+===============================================+================================================================================+
| Active Space Decompotion Method               | S\. M. Parker, T. Seideman, M. A. Ratner, and T. Shiozaki,                     |
|                                               | J. Chem. Phys. **139**, 021108 (2013).                                         |
+-----------------------------------------------+--------------------------------------------------------------------------------+



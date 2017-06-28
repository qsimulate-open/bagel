.. _asd_dmrg:

********
ASD-DMRG
********


Description
===========
ASD-DMRG algorithm is able to deal with molecular aggregates with more than two sites.

Multisite construction
======================
.. toctree::
   :maxdepth: 1

   multisite.rst


Keywords
========

.. topic:: ``method``

   | **Description:** Method to compute active subspaces.
   | **Datatype:** string
   | **Value:**
   |   ``ras``: use restricted active space configuration interaction method

.. topic:: ``ras``

   | **Description:** Specify restricted active space method settings.
   | **Recommendation:** See :ref:`rasci` for details.

.. topic:: ``restricted``

   | **Description:** Specify occupation restriction in active space.
   | **Recommendation:** See sample input for details.

.. topic:: ``spaces``

   | **Description:** Specify important single site states with the following keys:
   |  ``charge``, ``spin``, ``nstate``
   | **Recommendation:** See sample input for details.

.. topic:: ``nstates``

   | **Description:** Number of target states.
   | **Datatype:** int
   | **Default** 1

.. topic:: ``ntrunc``

   | **Description:** Number of dmrg states to keep.
   | **Datatype:** int

.. topic:: ``thresh``

   | **Description:** Threshold for convergence in Davidson diagonalization.
   | **Datatype:** double

.. topic:: ``maxiter``

   | **Description:** Maximum number of iterations for Davidson diagonalization.
   | **Datatype:** int
   | **Default:** 50

.. topic:: ``perturb``

   | **Description:** Initial perturbation value.
   | **Datatype:** double
   | **Default:** 0.001
   | **Recommendation:** Use default.

.. topic:: ``perturb_thresh``

   | **Description:** Threshold for energy convergence when perturbation is applied.
   | **Datatype:** double
   | **Default:** 0.0001
   | **Recommendation:** Use default.

.. topic:: ``perturb_min``

   | **Description:** Minimum perturbation to be applied.
   | **Datatype:** double
   | **Default:** 0.00001
   | **Recommendation:** Use default.


Example
=======
Here is a sample calculation of Helium trimer aggregate with ASD-DMRG.


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
       {"atom" :"He", "xyz" : [    0.00000000000,     0.00000000000,     0.00000000000] }
     ]
   },

   {
     "title" : "hf",
     "saveref" : "A"
   },

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : false,
     "cartesian" : false,
     "geometry" : [
       {"atom" :"He", "xyz" : [    0.00000000000,     0.00000000000,     3.00000000000] }
     ]
   },

   {
     "title" : "hf",
     "saveref" : "B"
   },

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : false,
     "cartesian" : false,
     "geometry" : [
       {"atom" :"He", "xyz" : [    0.00000000000,     0.00000000000,     6.00000000000] }
     ]
   },

   {
     "title" : "hf",
     "saveref" : "C"
   },

   {
     "title" : "multisite",
     "refs" : ["A", "B", "C"],
     "active" : [ [1, 2, 3, 4, 5] ],
     "hf" : {
       "thresh" : 1.0e-12
     },
     "localization" : {
       "max_iter" : 50,
       "thresh" : 1.0e-12
     }
   },

   {
     "title" : "asd_dmrg",
     "nstate" : 2,
     "ntrunc" : 100,
     "method" : "ras",
     "thresh" : 1.0e-8,
     "perturb" : 0.001,
     "perturb_min" : 1.0e-7,
     "perturb_thresh" : 1.0e-7,
     "ras" : {
       "nguess" : 5,
       "maxiter" : 50,
       "thresh" : 1.0e-8
     },
     "maxiter" : 50,
     "spaces" : [ [ {"charge" : 0, "nspin" : 0, "nstate" : 1},
                    {"charge" : 0, "nspin" : 2, "nstate" : 1},
                    {"charge" : 1, "nspin" : 1, "nstate" : 1},
                    {"charge" : -1, "nspin" : 1, "nstate" : 1} ] ],
     "restricted" : [ { "orbitals" : [ 1, 0, 4], "max_holes" : 1, "max_particles" : 1 } ]
   }

   ]}


Reference
=========

+------------------------------------------------+--------------------------------------------------------------------------------+
|          Description of Reference              |                          Reference                                             |
+================================================+================================================================================+
| Active space decomposition with multiple sites:| S\. M. Parker and T. Shiozaki, J. Chem. Phys. **141**, 211102 (2014).          |
| density matrix renormalization group algorithm |                                                                                |
+------------------------------------------------+--------------------------------------------------------------------------------+


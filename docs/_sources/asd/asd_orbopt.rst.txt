.. _asd_orbopt:

************************
ASD orbital optimization
************************


Description
===========
Active space decomposition with orbital optimization is implemented. Set ``"title" = "asd_orbopt"`` to start the input section. The two fragments of the molecule can be covalently linked, but the bridging orbitals cannot be active. Orbital rotations between active subspaces are included to partition the total active space into subspaces. The two active subspaces should be well separated. The calculation starts with dimer construction which is slightly different from that in :ref:`dimer_asd`.


Dimer Construction
==================
.. toctree::
   :maxdepth: 1

   dimer.rst


Keywords
========

Required Keywords
-----------------

.. topic:: ``algorithm``

   | **Description:** 'bfgs' orbital optimization algorithm is implemented.
   | **Datatype:** string
   | **Value:** bfgs

.. topic:: ``maxiter``

   | **Description:** Maximum number of orbital optimization iterations.
   | **Datatype:** int
   | **Default:** 50

.. topic:: ``asd``

   | **Description:** ``asd`` block specifies asd options.
   | **Recommendation:** See :ref:`dimer_asd` for details.

.. topic:: ``gradient_thresh``

   | **Description:** Threshold for calculating gradient.
   | **Datatype:** double
   | **Default:** 0.0001
   | **Recommendation:** Use default.

.. topic:: ``rotation_thresh``

   | **Description:** Threshold for orbital rotation.
   | **Datatype:** double
   | **Default:** 0.0001
   | **Recommendation:** Use default.

.. topic:: ``energy_thresh``

   | **Description:** Convergence threshold for calculating energy.
   | **Datatype:** double
   | **Default:** :math:`1.0\times 10^{-6}`
   | **Recommendation:** Use default.

.. topic:: ``semi-canonicalize``

   | **Description:** Whether to semi-canonicalize localized orbitals.
   | **Datatype:** bool
   | **Default:** false


Reference
=========
+-----------------------------------------------+--------------------------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                                   |
+===============================================+======================================================================================+
|  ASD orbital optimization                     | I\. Kim, S. M. Parker, and T. Shiozaki, J. Chem. Theory Comput. **11**, 3636 (2015). |
+-----------------------------------------------+--------------------------------------------------------------------------------------+


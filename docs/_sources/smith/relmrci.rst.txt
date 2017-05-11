.. _relmrci:

***************************************************************
Relativistic Multireference configuration interaction (RelMRCI)
***************************************************************


Description
===========
Relativistic version of the fully internally contracted multireference configuration interaction.
This is implemented with the automatical code generator SMITH3.
By default, the Davidson corrected (+Q) energy is also printed.

Keywords
========
CASSCF keywords
---------------
See :ref:`casscf`.

SMITH keywords
--------------

The default values are recommended unless mentioned otherwise.

.. topic:: ``title``

   | **Description:** Method to use in SMITH-generated code.
   | **Datatype:** string
   | **Values:**
   |    ``caspt2``: Request (Rel)CASPT2 calculation.
   |    ``mrci``: Request (Rel)MRCI calculation.
   |    ``casa``: Request (Rel)CAS/A calculation.
   | **Default:** N/A.

.. topic:: ``thresh``

   | **Description:** Convergence threshold.
   | **Datatype:** double precision
   | **Default:** 1.0e-8 (gradient) 1.0e-6 (otherwise)

.. topic:: ``thresh_overlap``

   | **Description:** Overlap threshold.
   | **Datatype:** double precision
   | **Default:** 1.0-9

.. topic:: ``frozen``

   | **Description**: Freeze core orbitals. 
   | **Datatype**: bool
   | **Default**: true

.. topic:: ``ncore``

   | **Description:** Number of frozen core orbitals.
   | **Datatype:** int 
   | **Default:** If ``frozen`` is true, subvalence orbitals are frozen. If false, zero. 

.. topic:: ``nfrozenvirt``

   | **Description:** Number of frozen virtual orbitals.
   | **Datatype:** int
   | **Default:** 0

.. topic:: ``maxiter``

   | **Description:** Maximum number of iterations in MRCI calculations.
   | **Datatype:** int
   | **Default:** 50

.. topic:: ``maxtile``

   | **Description:** Maximum number of orbitals in a single data tile used in SMITH3.
   | **Datatype:** int
   | **Default:** 10


Example
=======

References
==========

Fully internally contracted MRCI itself is not particularly new.

BAGEL References
----------------

+---------------------------------------------------+-------------------------------------------------------------------------------------+
|          Description of Reference                 |                         Reference                                                   |
+===================================================+=====================================================================================+
| SMITH3                                            | M\. K. MacLeod and T. Shiozaki, J. Chem. Phys. **142**, 010507 (2015)               |
+---------------------------------------------------+-------------------------------------------------------------------------------------+
|  Relativistic extension                           | T\. Shiozaki and W. Mizukami, J. Chem. Theory Comput. **11**, 4733 (2015)           |
+---------------------------------------------------+-------------------------------------------------------------------------------------+

General References
------------------
+---------------------------------------------------+-------------------------------------------------------------------------------------+
|          Description of Reference                 |                         Reference                                                   |
+===================================================+=====================================================================================+
|  Fully internally contracted MRCI                 | H\.-J\. Werner and E.-A. Reinsch, J. Chem. Phys. **76**, 3144 (1982)                |
+---------------------------------------------------+-------------------------------------------------------------------------------------+
|  Modern implementation                            | M\. Saitow, Y. Kurashige, and T. Yanai, J. Chem. Phys. **139**, 044118 (2013)       |
+---------------------------------------------------+-------------------------------------------------------------------------------------+


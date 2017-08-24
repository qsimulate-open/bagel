.. _relcaspt2:

************************************************************************
Relativistic multireference second-order perturbation theory (RelCASPT2)
************************************************************************


Description
===========
The relativistic version of the fully internally contracted CASPT2. As in its non-relativistic analogue, single-state calculations,
multi-state calculations, and its extended variant are available. Of course, this is implemented with the automatic code generator SMITH3.

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

.. topic:: ``ms``

   | **Description:** Do multistate CASPT2.
   | **Datatype:** bool
   | **Values:**
   |    ``true``: Do MS-CASPT2.
   |    ``false``: Do SS-CASPT2.
   | **Default:** true.

.. topic:: ``xms``

   | **Description:** Do extended multistate CASPT2.
   | **Datatype:** bool
   | **Values:**
   |    ``true``: Do XMS-CASPT2.
   |    ``false``: Do MS-CASPT2.
   | **Default:** true.

.. topic:: ``sssr``

   | **Description:** Use SS-SR contraction scheme.
   | **Datatype:** bool
   | **Values:**
   |    ``true``: use SS-SR contraction scheme.
   |    ``false``: use MS-MR contraction scheme.
   | **Default:** true.

.. topic:: ``shift``

   | **Description:** Vertical shift.
   | **Datatype:** double precision
   | **Default:** 0.0
   | **Recommendation:** Use default. Increase the value if the convergence problem presents.

.. topic:: ``thresh``

   | **Description:** Convergence threshold.
   | **Datatype:** double precision
   | **Default:** For single point energy calculation, 1.0e-6. Tight convergence for gradient calculation, 1.0e-8.

.. topic:: ``thresh_overlap``

   | **Description:** Overlap cutoff threshold for internally contracted basis.
   | **Datatype:** Double precision
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

   | **Description:** Maximum number of iterations in CASPT2 calculations.
   | **Datatype:** int
   | **Default:** 50

.. topic:: ``maxtile``

   | **Description:** Maximum number of orbitals in a single data tile used in CASPT2.
   | **Datatype:** int
   | **Default:** 10

.. topic:: ``cimaxtile``

   | **Description:** Maximum number of Slater determinants in a single data tile used in CASPT2 gradient.
   | **Datatype:** int
   | **Default:** 100 (When number of determinants is over 10000), 10 (otherwise)
   | **Recommendation:** Use default. Increase further when the number of determinants is larger.


Example
=======
XMS-CASPT2 calculation based on the two-state CASSCF reference function, with vertical shift of 0.2 :math:`E_h`. "SS-SR" contraction scheme is used.
The active space of (6e,6o), which comprises three :math:`\pi` and three :math:`\pi^*` orbitals, is used.

Sample input
------------

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
    "title" : "caspt2",
    "smith" : {
      "method" : "caspt2",
      "ms" : true,
      "xms" : true,
      "sssr" : true,
      "shift" : 0.2
    }
    "nstate" : 2,
    "nact" : 6,
    "nclosed" : 18,
    "active" : [17, 20, 21, 22, 23, 30]
  }
  ]}

References
==========

BAGEL References
----------------
+---------------------------------------------------+----------------------------------------------------------------------------------------------+
|          Description of Reference                 |                          Reference                                                           |
+===================================================+==============================================================================================+
| SMITH3                                            | M\. K. MacLeod and T. Shiozaki, J. Chem. Phys. **142**, 010507 (2015).                       |
+---------------------------------------------------+----------------------------------------------------------------------------------------------+
| Relativistic CASPT2                               | T\. Shiozaki and W. Mizukami, J. Chem. Theory Comput. **11**, 4733 (2015).                   |
+---------------------------------------------------+----------------------------------------------------------------------------------------------+

General References
------------------
+---------------------------------------------------+-------------------------------------------------------------------------------------------------------+
|          Description of Reference                 |                          Reference                                                                    |
+===================================================+=======================================================================================================+
| CASPT2                                            | K\. Andersson, P.-Å. Malmqvist, and B. O. Roos, J. Chem. Phys. **96**, 1218 (1992).                   |
+---------------------------------------------------+-------------------------------------------------------------------------------------------------------+
| MS-CASPT2                                         | J\. Finley, P.-Å. Malmqvist, B. O. Roos, and L. Serrano-Andres, Chem. Phys. Lett. **288**, 299 (1998).|
+---------------------------------------------------+-------------------------------------------------------------------------------------------------------+
| Extended multiconfigurational perturbation theory | A\. A. Granovsky, J. Chem. Phys. **134**, 214113 (2011).                                              |
+---------------------------------------------------+-------------------------------------------------------------------------------------------------------+
| XMS-CASPT2                                        | T\. Shiozaki, W. Győrffy, P. Celani, and H.-J. Werner, J. Chem. Phys. **135**, 081106 (2011).         |
+---------------------------------------------------+-------------------------------------------------------------------------------------------------------+

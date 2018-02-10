.. _relcaspt2:

************************************************************************
Relativistic multireference second-order perturbation theory (RelCASPT2)
************************************************************************


Description
===========
The relativistic version of the fully internally contracted CASPT2. As in its non-relativistic analogue, single-state calculations,
multi-state calculations, and its extended variant are available. Of course, this is implemented with the automatic code generator SMITH3.

Title: ``relsmith``

Keywords
========
CASSCF keywords
---------------
See :ref:`zcasscf`.

SMITH keywords
--------------

The default values are recommended unless mentioned otherwise.

.. topic:: ``method``

   | **Description:** Do multistate CASPT2.
   | **Datatype:** string
   | **Values:**
   |    ``caspt2``: Standard CASPT2.
   |    ``casa``  : Use Dyall's Hamiltonian

.. topic:: ``ms``

   | **Description:** Do multistate CASPT2.
   | **Datatype:** bool
   | **Default:** true.

.. topic:: ``xms``

   | **Description:** Do extended multistate CASPT2.
   | **Datatype:** bool
   | **Default:** true.

.. topic:: ``sssr``

   | **Description:** Use SS-SR contraction scheme.
   | **Datatype:** bool
   | **Default:** true.

.. topic:: ``shift``

   | **Description:** Vertical shift.
   | **Datatype:** double precision
   | **Default:** 0.0

.. topic:: ``thresh``

   | **Description:** Convergence threshold.
   | **Datatype:** double precision
   | **Default:** For single point energy calculation, 1.0e-6. Tight convergence for gradient calculation, 1.0e-8.

.. topic:: ``thresh_overlap``

   | **Description:** Overlap cutoff threshold for internally contracted basis.
   | **Datatype:** Double precision
   | **Default:** 1.0e-9

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

.. topic:: ``block_diag_fock``

   | **Description:** Using a block-diagonal zeroth-order Hamiltonian
   | **Datatype:** bool
   | **Default:** false.

.. topic:: ``maxiter``

   | **Description:** Maximum number of iterations in CASPT2 calculations.
   | **Datatype:** int
   | **Default:** 50

.. topic:: ``maxtile``

   | **Description:** Maximum number of orbitals in a single data tile used in CASPT2.
   | **Datatype:** int
   | **Default:** 10

.. topic:: ``davidson_subspace``

   | **Description:**  Number of vectors retained in the limited-memory Davidson algorithm.
   | **Datatype:** int
   | **Default:** 10
   | **Recommendation:** Reduce if an insufficient amount of memory is available (do not reduce to a value lower than 3). 


Example
=======

Sample input
------------

.. code-block:: javascript

   { "bagel" : [
   
   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : true,
     "geometry" : [
       { "atom" : "O",  "xyz" : [ 0.000, 0.000, 1.210 ]},
       { "atom" : "O",  "xyz" : [ 0.000, 0.000, 0.000 ]}
     ]
   },
   
   {
     "title" : "hf",
     "charge" : "+2"
   },
   
   {
     "title" : "dhf",
     "charge" : "+2"
   },
   
   {
     "title"  : "zcasscf",
     "algorithm" : "second",
     "charge" : "0",
     "state" : [0, 0, 1],
     "thresh" : 1.0e-7,
     "nclosed"  : 7,
     "nact"   : 2
   },
   
   {
     "title" : "relsmith",
     "method" : "caspt2",
     "sssr" : "false",
     "ncore" : 2,
     "ms" : "true",
     "xms" : "true",
     "shift" : "0.0"
   }
   
   ]}

from which one obtains:

.. code-block:: javascript

    * Zeroth order energy : state  0       -0.4311850266
    * Zeroth order energy : state  1       -0.4310720127
    * Zeroth order energy : state  2       -0.4310720127

      ---- iteration ----

        0    -0.40810539     0.00207948      1.70
        1    -0.40848834     0.00016074      1.74
        2    -0.40849019     0.00000849      1.71
        3    -0.40849020     0.00000034      1.70

        0    -0.40811484     0.00207947      1.71
        1    -0.40849779     0.00016074      1.68
        2    -0.40849965     0.00000849      1.69
        3    -0.40849965     0.00000034      1.75

        0    -0.40811484     0.00207947      1.73
        1    -0.40849779     0.00016074      1.69
        2    -0.40849965     0.00000849      1.69
        3    -0.40849965     0.00000034      1.71

      -------------------


    * RelCASPT2 energy : state  0     -149.9873960122
    * RelCASPT2 energy : state  1     -149.9873936904
    * RelCASPT2 energy : state  2     -149.9873936904

    * MS-RelCASPT2 Heff
      (-149.9873960122,-0.0000000000)(0.0000000000,0.0000000000)(0.0000000000,-0.0000000000)
      (0.0000000000,-0.0000000000)(-149.9873936904,-0.0000000000)(0.0000000000,0.0000000000)
      (0.0000000000,0.0000000000)(0.0000000000,-0.0000000000)(-149.9873936904,-0.0000000000)


    * MS-RelCASPT2 rotation matrix
      (1.0000000000,0.0000000000)(0.0000000000,0.0000000000)(0.0000000000,0.0000000000)
      (0.0000000000,0.0000000000)(-0.7112091498,0.4043385965)(0.0936202837,-0.5673861886)
      (0.0000000000,0.0000000000)(-0.4948015439,-0.2930243620)(-0.7669166768,0.2848630658)

    * MS-RelCASPT2 energy : state  0     -149.9873960122
    * MS-RelCASPT2 energy : state  1     -149.9873936904
    * MS-RelCASPT2 energy : state  2     -149.9873936904


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
| XMCQDPT2                                          | A\. A. Granovsky, J. Chem. Phys. **134**, 214113 (2011).                                              |
+---------------------------------------------------+-------------------------------------------------------------------------------------------------------+
| XMS-CASPT2                                        | T\. Shiozaki, W. Győrffy, P. Celani, and H.-J. Werner, J. Chem. Phys. **135**, 081106 (2011).         |
+---------------------------------------------------+-------------------------------------------------------------------------------------------------------+

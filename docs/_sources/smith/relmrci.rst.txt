.. _relmrci:

***************************************************************
Relativistic Multireference configuration interaction (RelMRCI)
***************************************************************


Description
===========
Relativistic version of the fully internally contracted multireference configuration interaction.
This is implemented with the automatical code generator SMITH3.
By default, the Davidson corrected (+Q) energy is also printed.

Title: ``relsmith``, method: ``mrci``

Keywords
========
CASSCF keywords
---------------
See :ref:`zcasscf`.

SMITH keywords
--------------

The default values are recommended unless mentioned otherwise.

.. topic:: ``thresh``

   | **Description:** Convergence threshold.
   | **Datatype:** double precision
   | **Default:** 1.0e-8 (gradient) 1.0e-6 (otherwise)

.. topic:: ``thresh_overlap``

   | **Description:** Overlap threshold.
   | **Datatype:** double precision
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

.. topic:: ``maxiter``

   | **Description:** Maximum number of iterations in MRCI calculations.
   | **Datatype:** int
   | **Default:** 50

.. topic:: ``maxtile``

   | **Description:** Maximum number of orbitals in a single data tile used in SMITH3.
   | **Datatype:** int
   | **Default:** 10

.. topic:: ``davidson_subspace``

   | **Description:**  Number of vectors retained in the limited-memory Davidson algorithm.
   | **Datatype:** int
   | **Default:** 10
   | **Recommendation:** Reduce if an insufficient amount of memory is available (do not reduce to a value lower than 3). 


Example
=======

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
     "method" : "mrci",
     "ncore" : 2 
   }
   
   ]}

from which one obtains

.. code-block:: javascript


    * Zeroth order energy : state  0       -0.4311850170
    * Zeroth order energy : state  1       -0.4310720094
    * Zeroth order energy : state  2       -0.4310720094

      ---- iteration ----

        0   0  -149.90092252     0.00065697     67.36
        0   1  -149.90091061     0.00065698      0.41
        0   2  -149.90091061     0.00065698      0.41

        1   0  -149.92100585     0.00012080     69.68
        1   1  -149.92099523     0.00012080      0.37
        1   2  -149.92099523     0.00012080      0.39

        2   0  -149.92179156     0.00003649     68.67
        2   1  -149.92178091     0.00003649      0.36
        2   2  -149.92178091     0.00003649      0.37

        3   0  -149.92186506     0.00000959     70.71
        3   1  -149.92185440     0.00000959      0.36
        3   2  -149.92185440     0.00000959      0.37

        4   0  -149.92187033     0.00000231     70.11
        4   1  -149.92185967     0.00000230      0.39
        4   2  -149.92185967     0.00000230      0.39

        5   0  -149.92187064     0.00000056     68.53
        5   1  -149.92185998     0.00000055      0.00
        5   2  -149.92185998     0.00000055      0.00


      -------------------

       - MRCI energy evaluation                  420.87

        0   0  -149.95278455     0.00000000      0.00
        0   1  -149.95277356     0.00000000      0.00
        0   2  -149.95277356     0.00000000      0.00

       - MRCI+Q energy evaluation                  0.44

    * METHOD: RELSMITH                           422.34



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


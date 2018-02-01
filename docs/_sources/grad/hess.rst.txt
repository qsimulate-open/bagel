.. _hess:

****************************************
Molecular Hessian and frequency analysis
****************************************

Description
===========
The Hessian section can be used to compute the numerical Hessian by central gradient differences. The Hessian, mass weighted Hessian, and symmetrized mass weighted Hessian are printed in the output by default. The rotational and translational degrees of freedom have been projected out. Vibrational frequencies, infrared intensities, and the Cartesian eigenvectors of each normal mode are also computed. The masses are averaged over the natural occurrence of isotopes unless otherwise specified (see other keywords below).

Keywords
========
Required Keywords
-----------------
.. topic:: ``hessian``

   | **Description:** Requests that the Numerical Hessian be computed

.. topic:: ``method``

   | **Description:** The method array allows the user to specify one or more methods to be used in the Hessian calculation. See section on input structure for more information.

.. topic:: ``dx``

   | **Description:** The step size used in the displacements in the gradient difference calculations. The units are bohr.
   | **Datatype:** double precision
   | **Default:** 1.0e-3

Optional Keywords
-----------------

.. topic:: ``nproc``

   | **Description:** The Hessian code is embarrassingly parallelized so that the displacements in the central gradient difference calculations can be run at the same time. The nproc keyword allows the user to specify the number of MPI processes to be used for each gradient calculation.
   | **Datatype:** int
   | **Default:** 1

Other Keywords
--------------

.. topic:: ``mass``

   | **Description:** Atomic masses can be specified by the users in the atom line in the molecule block. 
   | **Datatype:** double
   | **Default:** averaged over the natural occurrence of isotopes 


Example
=======

Sample input
------------
A sample input for the benzene molecule.

.. code-block:: javascript

  { "bagel" : [

  {
    "title" : "molecule",
    "basis" : "svp",
    "df_basis" : "cc-pvdz-jkfit",
    "angstrom" : false,
    "geometry" : [
    { "atom" : "C", "xyz" : [     -0.072972,      2.554311,     -0.000005 ] },
    { "atom" : "C", "xyz" : [      2.551272,      2.554150,     -0.000005 ] },
    { "atom" : "C", "xyz" : [      3.863393,      4.826435,     -0.001030 ] },
    { "atom" : "C", "xyz" : [      2.551429,      7.099204,     -0.002139 ] },
    { "atom" : "C", "xyz" : [     -0.072483,      7.099365,     -0.002170 ] },
    { "atom" : "C", "xyz" : [     -1.384772,      4.826799,     -0.001078 ] },
    { "atom" : "H", "xyz" : [     -1.095603,      0.783105,      0.000799 ] },
    { "atom" : "H", "xyz" : [      3.573593,      0.782757,      0.000939 ] },
    { "atom" : "H", "xyz" : [      5.908616,      4.826483,     -0.000944 ] },
    { "atom" : "H", "xyz" : [      3.574399,      8.870211,     -0.002887 ] },
    { "atom" : "H", "xyz" : [     -1.095141,      8.870560,     -0.003114 ] },
    { "atom" : "H", "xyz" : [     -3.429999,      4.827102,     -0.001033 ] }
    ]
  },

  {
    "title" : "hf"
  },

  {
    "title" : "casscf",
    "nstate" : 2,
    "nclosed" : 18,
    "nact" : 6,
    "active" : [17, 20, 21, 22, 23, 30]
  },

  {
    "title" : "hessian",
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
       "nact" : 6,
       "nclosed" : 18
    } ]
  }

  ]}

If you are running a Hessian calculation using the embarassingly parallel implementation, it is recommended to only have the Hessian calculation in your input.
A molden file generated from a previous calculation can be read at the start of the calculation.

.. code-block:: javascript

  { "bagel" : [

  {
    "title" : "molecule",
    "basis" : "molden",
    "df_basis" : "cc-pvdz-jkfit",
    "molden_file" : "restart.molden"
  },

  {
    "title" : "hessian",
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
       "nact" : 6,
       "nclosed" : 18
    } ]
  }

  ]}

The atomic masses can be specified as shown below.

.. code-block:: javascript

    { "atom" : "C", "xyz" : [     -0.072972,      2.554311,     -0.000005 ],  "mass" : 12.0 },


References
==========

+----------------------------------------------------+----------------------------------------------------------------------------------------------------------+
|          Description of Reference                  |                          Reference                                                                       |
+====================================================+==========================================================================================================+
| General description of vibrational spectroscopy    | E\. B. Wilson, Jr., J. C. Decius, and P. C. Cross, *Molecular Vibrations* (McGraw-Hill, New York, 1955). |
+----------------------------------------------------+----------------------------------------------------------------------------------------------------------+


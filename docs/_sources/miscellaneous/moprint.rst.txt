.. _moprint:

****************************************
Printing orbital densities to cube files
****************************************

===========
Description
===========

This module prints molecular orbital densities in the Gaussian Cube format.  
It can be used to view the shape and extent of relativistic or gauge-including molecular orbitals, which cannot 
be viewed in Molden due to the use of complex basis functions or the four-component framework.  

A separate .cub file is generated for each printed orbital, plus one for the total electron density.  
The axis vectors are simply the Cartesian *x*, *y*, and *z* axes.

In addition, one can print the relaxed electron density in MP2, SA-CASSCF, and CASPT2 force calculations as well.
See :ref:`nuclear gradient and derivative coupling <force>` and the example below.

Command: ``moprint``

========
Keywords
========

.. topic:: ``paired``

   | **Description:** Determined whether we plot spatial MOs (true) or spin MOs.  
   | **Datatype:** bool
   | **Default:** True, unless we are printing 4-component orbitals generated with an external magnetic field.   
   | **Recommendation:** Use the default.

.. topic:: ``orbitals``

   | **Description:**  Indices of the molecular orbitals to be printed.
   | **Datatype:** Vector of integers
   | **Default:** Prints the active orbitals from CASSCF, and the frontier orbitals from Hartree--Fock

.. topic:: ``start_pos``

   | **Description:** Coordinates for one corner of the box within which densities are printed.
   | **Datatype:** Array of 3 doubles
   | **Default:** A position is chosen so that all atoms (except the dummy atoms) are at least :math:`4 a_0` from the edges of the box.
   | **Recommendation:** Use the default.

.. topic:: ``inc_size``

   | **Description:** Distances between adjacent gridpoints in each of the three dimensions.
   | **Datatype:** Array of 3 doubles
   | **Default:** :math:`0.25 a_0` in each direction.
   | **Recommendation:** Use the default.

.. topic:: ``angstrom``

   | **Description:** Unit of the "inc_size" parameter
   | **Datatype:** bool
   | **Default:** False (meaning Bohr; set to true for angstrom)

.. topic:: ``mo_filename``

   | **Description:** Name of the MO cube file
   | **Datatype:** string
   | **Default:** "mo"

.. topic:: ``density_filename``

   | **Description:** Name of the density cube file
   | **Datatype:** string
   | **Default:** "density"

=======
Example
=======

Sample input: Print MO
----------------------

Write molecular orbitals:

.. code-block:: javascript

   { "bagel" : [

   .... energy calculation....

   {
     "title" : "moprint",
     "inc_size" : [ 0.20, 0.20, 0.20 ],
     "orbitals" : [ 14, 15, 16, 17, 18, 19, 20, 21, 22 ]
   }

   ]}

Sample input: Print relaxed density
-----------------------------------

Write relaxed density to the file ``density_0.cub`` from the XMS-CASPT2 force calculation:

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
      "density_filename" : "density_0",
      "inc_size" : [ 0.20, 0.20, 0.20 ]
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

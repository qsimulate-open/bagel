.. _molden:

******
MOLDEN
******

===========
Description
===========
A molden file format can be generated as an output for orbitals, geometries, and vibrational frequencies. If a geometry optimization is performed, a molden file is generated that contains the geometry at each step in the optimization by default.

A molden file can also be read to restart calculations. See the molecule section for more information.

========
Keywords
========

.. topic:: ``file``

   | **Description:** Name of the file to be generated
   | **Datatype:** string

.. topic:: ``orbitals``

   | **Description:** Molecular orbitals are written to the molden file.
   | **Datatype:** bool
   | **Default:** false

.. topic:: ``vibration``

   | **Description:** The vibrational frequencies, infrared intensities, and the cartesian eigenvectors of the normal modes are written to the molden file.
   | **Datatype:** bool
   | **Default:** false

=======
Example
=======

Sample input
------------

Write molecular orbitals:

.. code-block:: javascript

   { "bagel" : [

   .... energy calculation....

   {
     "title" : "print",
     "file" : "orbitals.molden",
     "orbitals" : true
   }

   ]}

Write results from a vibrational frequency calculation:

.. code-block:: javascript

   { "bagel" : [

   .... energy and Hessian calculation ....

   {
     "title" : "print",
     "file" : "freq.molden",
     "vibration" : true
   }

   ]}

One can use a Molden file as an input as follows. Note that the Molden format does not specify the density fitting basis set,
and therefore, the molecule block does need to contain the df_basis keyword.

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "molden",
     "df_basis" : "svp-jkfit",
     "molden_file" : "hf.molden"
   }

   ]}

==========
References
==========

+----------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
|          Description of Reference                  |                          Reference                                                                        |
+====================================================+===========================================================================================================+
| Molden and the model file format                   |   www.cmbi.ru.ml/molden                                                                                   |
+----------------------------------------------------+-----------------------------------------------------------------------------------------------------------+


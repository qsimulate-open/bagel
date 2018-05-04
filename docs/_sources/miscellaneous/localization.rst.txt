.. _localization:

********************
Orbital localization
********************

===========
Description
===========
Localized molecular orbitals can be generated using the Pipek-Mezey (PM) or regional localized molecular orbital (RLMO) procedures.

========
Keywords
========

.. topic:: ``occupied``

   | **Description:** Localize the occupied orbitals
   | **Datatype:** bool
   | **Default:** true

.. topic:: ``active``

   | **Description:** Localize the active orbitals
   | **Datatype:** bool
   | **Default:** true

.. topic:: ``virtual``

   | **Description:** Localize the virtual orbitals
   | **Datatype:** bool
   | **Default:** false

.. topic:: ``algorithm``

   | **Description:** The localization scheme being used.
   | **Datatype:** string
   | **Values:**
   |    ``pm``: Uses Pipek-Mezey localization
   |    ``region`` : Orthogonalize based on regions
   | **Default:** pm
   | **Recommendation:** Defining regions is particularly useful when studying dimers or trimers. For standard cases, use default

.. topic:: ``type``

   | **Description:** The type of localzation used in the Pipek-Mezey localization
   | **Datatype:** string
   | **Values:**
   |    ``region`` : localize to a region (a group of atoms defined by the user)
   |    ``atomic``: localize to the atoms
   | **Default:** atomic
   | **Recommendation:** Defining regions is particularly useful when studying dimers or trimers. For standard cases, use default

.. topic:: ``region_size``

   | **Description:** Define the regions used if type is set to region.
   | **Datatype:** vector<int>
   | **Values:** Define vector. For example, 3 regions of 2 atoms each would be [2,2,2]

.. topic:: ``lowdin``

   | **Description:** Lowdin charges are used in the localization
   | **Datatype:** bool
   | **Values:**
   |    ``true``: Uses Lowdin charges
   |    ``false`` : Uses Mulliken charges
   | **Default:** true
   | **Recommendation:** : Use default

=======
Example
=======

Sample input
------------

.. code-block:: javascript

   { "bagel" : [

   .... energy calculation....

   {
     "title" : "localize",
     "virtual" : true,
     "algorithm" : "pm"
   }

   ]}

==========
References
==========

+----------------------------------------------------+----------------------------------------------------------------------------------------------+
|          Description of Reference                  |                          Reference                                                           |
+====================================================+==============================================================================================+
| Pipek-Mezey orbital localization                   | J\. Pipek and P. G. Mezey, J. Chem. Phys. **90**, 4916 (1989).                               |
+----------------------------------------------------+----------------------------------------------------------------------------------------------+
| Orthogonalize based on regions                     | P\. de Silva, M. Giebultowski, and J. Korchowiec, Phys. Chem. Chem. Phys. **14**, 546 (2012).|
+----------------------------------------------------+----------------------------------------------------------------------------------------------+



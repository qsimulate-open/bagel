.. _mp2:

*****************************************
Møller--Plesset perturbation theory (MP2)
*****************************************

===========
Description
===========
Second-order Møller--Plesset perturbation theory (MP2) calculations are performed with density fitting using
the keyword ``"title" : "mp2"``.

========
Keywords
========

.. topic:: ``frozen``

   | **Description**: to have frozen orbitals or not
   | **Datatype**: bool
   | **Default**: true

.. topic:: ``ncore``

   | **Description**: manually specify number of frozen orbitals, used when 'frozen' is turned on
   | **Datatype**: int

.. topic:: ``aux_basis``

   | **Description**: specify an auxiliary basis set for MP2
   | **Datatype**: string
   | **Default**: use the same density fitting basis as in :ref:`molecule`
   | **Recommendation**: use MP2-fit auxiliary basis (auxiliary basis ends with 'ri')

=======
Example
=======

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "cc-pvdz",
     "df_basis" : "cc-pvdz-jkfit",
     "angstrom" : "true",
     "geometry" : [
       { "atom" : "C", "xyz" : [ -1.20433891360,  0.54285096106, -0.04748199659] },
       { "atom" : "C", "xyz" : [ -1.20543291352, -0.83826393986,  0.12432899108] },
       { "atom" : "C", "xyz" : [ -0.00000600000, -1.52953889027,  0.20833398505] },
       { "atom" : "C", "xyz" : [  1.20544091352, -0.83825393987,  0.12432799108] },
       { "atom" : "C", "xyz" : [  1.20433091360,  0.54284396106, -0.04748099659] },
       { "atom" : "C", "xyz" : [  0.00000400000,  1.23314191154, -0.13372399041] },
       { "atom" : "H", "xyz" : [ -2.13410484690,  1.07591192282, -0.12500499103] },
       { "atom" : "H", "xyz" : [ -2.13651384673, -1.37179190159,  0.18742198655] },
       { "atom" : "H", "xyz" : [  0.00000000000, -2.59646181374,  0.33932597566] },
       { "atom" : "H", "xyz" : [  2.13651384673, -1.37179290159,  0.18742198655] },
       { "atom" : "H", "xyz" : [  2.13410684690,  1.07591292282, -0.12500599103] },
       { "atom" : "H", "xyz" : [ -0.00000000000,  2.29608983528, -0.28688797942] }
     ]
   },

   {
     "title" : "hf"
   },

   {
     "title" : "mp2",
     "aux_basis" : "cc-pvdz-ri",
     "frozen" : true
   }

   ]}

The SCF calculation should converge after 11 iterations with energy -230.72159477. The DF-MP2 output is as follows:

.. code-block:: javascript

  === DF-MP2 calculation ===

    * freezing 6 orbitals
  Number of auxiliary basis functions:      420

  Since a DF basis is specified, we compute 2- and 3-index integrals:
    o Being stored without compression. Storage requirement is 0.044 GB
       - 3-index ints prep                         0.00
       - 3-index ints                              0.06
       - 2-index ints                              0.00
       - computing inverse                         0.02
        elapsed time:        0.09 sec.


  Number of basis functions:      114
  Number of electrons      :       42

    * 3-index integral transformation done
    * ncache = 20
    * assembly done

      MP2 correlation energy:   -0.7813195576      0.19

      MP2 total energy:       -231.5029143295


    * METHOD: MP2                                  0.34

==========
References
==========

+-----------------------------------------------+-----------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                    |
+===============================================+=======================================================================+
| Original reference for MP2                    | C\. Møller and M. S. Plesset, Phys. Rev. **46**, 618 (1934).          |
+-----------------------------------------------+-----------------------------------------------------------------------+


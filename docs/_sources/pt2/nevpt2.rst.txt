.. _nevpt2:

*****************************************************
N-electron valence state perturbation theory (NEVPT2)
*****************************************************

===========
Description
===========
Calculations using the strongly contracted state-specific n-electron valence state perturbation theory (NEVPT2)
are done using the keyword ``"title" : "nevpt2"``.

========
Keywords
========

The default values are recommended unless mentioned otherwise.

.. topic:: ``frozen``

   | **Description:** to have frozen orbitals or not
   | **Datatype:** bool
   | **Default:** true.

.. topic:: ``ncore``

   | **Description**: manually specify number of frozen orbitals, used when 'frozen' is turned on with ``"frozen" : "true"``.
   | **Datatype**: int

.. topic:: ``aux_basis``

   | **Description**: specify an auxiliary basis set for MP2
   | **Datatype**: string
   | **Default**: use the same density fitting basis as in :ref:`molecule`
   | **Recommendation**: use MP2-fit auxiliary basis (auxiliary basis ends with 'ri')

.. topic:: ``istate``

   | **Description**: specific state used in the evaluation of the dynamical correlation
   | **Datatype**: int
   | **Default**: 0 (ground state)

=======
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
       { "atom" : "O",  "xyz" : [ 0.00000000000, 0.0,  0.000000000000]},
       { "atom" : "H",  "xyz" : [ 1.45860189536, 0.0,  0.504283963824]},
       { "atom" : "H",  "xyz" : [ 0.75860194558, 0.0, -0.504283963824]}
     ]
   },

   {
     "title" : "nevpt2",
     "nact" : 2,
     "nclosed" : 4,
     "frozen" : true,
     "thresh" : 1.0e-8,
     "thresh_scf" : 1.0e-8,
     "thresh_fci" : 1.0e-10
   }

   ]}


from which one obtains

.. code-block:: javascript

  === DF-NEVPT2 calculation ===

       - RDMs, hole RDMs, others                   0.00
       - Fock computation                          0.01
    * 3-index integral transformation done
       - K matrices                                0.00
       - A, B, C, and D matrices                   0.00
    * ncache = 20
    * assembly done

      NEVPT2 correlation energy:   -0.1703158878      0.02

          (+0)     -0.0977514346
          (+0)'    -0.0257993320
          (+1)     -0.0086109214
          (+1)'    -0.0000000000
          (+2)     -0.0007336183
          (-1)     -0.0341475759
          (-1)'     0.0000000000
          (-2)     -0.0032730057

      NEVPT2 total energy:        -76.0205249507



==========
References
==========

BAGEL references
================
+-----------------------------------------------+-----------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                    |
+===============================================+=======================================================================+
| Relativistic implementation of NEVPT2 in      | T\. Shiozaki and W. Mizukami, J. Chem. Theory Comput. **11**, 4733    |
| BAGEL                                         | (2015).                                                               |
+-----------------------------------------------+-----------------------------------------------------------------------+

General references
==================

+-----------------------------------------------+-----------------------------------------------------------------------+
|          Description of Reference             |                          Reference                                    |
+===============================================+=======================================================================+
| Original reference  for NEVPT2                | C\. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, and J.-P.    |
|                                               | Malrieu, J. Chem. Phys. **114**, 10252 (2001).                        |
+-----------------------------------------------+-----------------------------------------------------------------------+
| Spin-free formulation                         | C\. Angeli, R. Cimiraglia, and J.-P. Malrieu, J. Chem. Phys. **117**, |
|                                               | 9138 (2002).                                                          |
+-----------------------------------------------+-----------------------------------------------------------------------+


.. _molecule:

Required keywords
=================

.. topic:: ``geometry``

   | **Description**: Specify atoms and their Cartesian coordinates
   | **Datatype**: vector
   | **Values**:
   |    Vector of atoms provided in the following format ``{ "atom" : "atom symbol",  "xyz" : [x, y, z] }``
        Please see the end of the file for some examples.

.. topic:: ``basis``

   | **Description**: Define default basis set used for the system
   | **Datatype**: string
   | **Values**:
   |    Please refer to `Basis sets`_ and `Effective core potential (ECP) basis sets`_ for possible arguments.
        `User defined basis sets`_ are also possible.

.. topic:: ``df_basis``

   | **Description**: Basis sets used for density fitting
   | **Datatype**: string
   | **Values**:
   |     Please refer to `Density fitting basis sets`_ for possible arguments

.. note::
   The use of mixed basis sets and/or density fitting basis sets is possible by specifying a different
   basis set other than the default for each atom (see example for `Basis sets`_ below).

Optional keywords
=================

.. topic:: ``angstrom``

   | **Description**: Specify units for atomic coordinates (Angstrom or Bohr)
   | **Datatype**: bool
   |    ``true``: use Angstrom
   |    ``false``: use Bohr
   | **Default**: false (Bohr)

.. topic:: ``finite_nucleus``

   | **Description**: Represent the nucleus as a Gaussian charge distribution with default exponents
   | **Datatype**: bool
   | **Default**: false
   | **Note**:
   |    Within the ``geometry`` block, the ``exponent`` keyword can be used to specify a different exponent for the charge distribution of a particular atom,
        where a value of 0.0 indicates a point charge.  

.. topic:: ``molden_file``

   | **Description**: Filename of the molden file, which is required if ``"basis" : "molden"`` is specified.
   | **Datatype**: string
   | **Recommendation**: restarting from a molden file is not recommended, which is nevertheless sometimes useful.

.. topic:: ``cfmm``

   | **Description**: Turn on RHF-FMM; for more details, refer to :ref:`hf` section.
   | **Datatype**: bool
   | **Default**: false
   | **Recommendation**: Use for calculations on very large systems.

.. topic:: ``schwarz_thresh``

   | **Description**: Schwarz screening integral threshold, only used in RHF-FMM ``"cfmm" : "true"``.
                      For more details, refer to :ref:`hf` section.
   | **Datatype**: double
   | **Default**: :math:`1.0\times 10^{-12}`
   | **Recommendation**: Default, looser thresholds reduce accuracy but potentially increase speed.

.. topic:: ``dkh``

   | **Description**: Option to use the second-order Douglas--Kroll--Hess Hamiltonian (DKH2).
   | **Datatype**: bool
   | **Default**: false

.. topic:: ``magnetic_field``

   | **Description**: External magnetic field.  External magnetic fields are available only for the :ref:`hf`, :ref:`dhf`, and :ref:`zcasscf` modules.
   | **Datatype**: Array of three doubles (x, y, z)
   | **Default**: (0.0, 0.0, 0.0)

.. topic:: ``tesla``

   | **Description**: External magnetic field in units of Tesla.
   | **Datatype**: bool
   | **Default**: false (i.e., atomic units are used)

.. topic:: ``basis_type``

   | **Description**: Specify the type of atomic orbital basis functions,
        either standard Gaussian functions or gauge-including atomic orbitals (GIAOs).
   | **Datatype**: string
   | **Values**: ``gaussian`` / ``giao``, ``london``
   | **Default**: ``gaussian`` at zero magnetic field; ``giao`` when a field is applied

.. topic:: ``skip_self_interaction``

   | **Description**: Skip the electrostatic interactions between the dummy atoms.
   | **Datatype**: bool
   |    ``true``: skip the electrostatic interactions between the dummies.
   |    ``false``: explicitly calculate the electrostatic interactions between the dummies.
   | **Default**: true


Basis sets
==========

==================
Orbital basis sets 
==================

The following basis sets are available in BAGEL library. The basis set name can be used with the ``basis`` keyword.

.. hlist::
   :columns: 3

   * sto-3g
   * 3-21g
   * 6-31g
   * svp
   * tzvpp
   * qzvpp
   * cc-pvdz
   * cc-pvtz
   * cc-pvqz
   * cc-pv5z
   * cc-pv6z
   * cc-pcvdz
   * cc-pcvtz
   * cc-pcvqz
   * cc-pcv5z
   * cc-pcvdz-dk
   * cc-pcvtz-dk
   * aug-cc-pvdz
   * aug-cc-pvtz
   * aug-cc-pvqz
   * aug-cc-pv5z
   * aug-cc-pv6z
   * aug-cc-pcvdz
   * aug-cc-pcvtz
   * aug-cc-pcvqz
   * aug-cc-pcv5z
   * aug-cc-pcvdz-dk
   * aug-cc-pcvtz-dk
   * aug-cc-pcvqz-dk
   * aug-cc-pwcvdz
   * aug-cc-pwcvtz
   * aug-cc-pwcvqz
   * aug-cc-pwcv5z
   * d-aug-cc-pvdz
   * d-aug-cc-pvtz
   * d-aug-cc-pvqz
   * d-aug-cc-pv5z
   * ano-rcc

==========================
Density fitting basis sets
==========================

The following density fitting basis sets are available in BAGEL library. The basis set name can be used with the ``df_basis`` keyword.

.. hlist::
   :columns: 3

   * svp-jkfit
   * tzvpp-jkfit
   * qzvpp-jkfit
   * cc-pvdz-jkfit
   * cc-pvtz-jkfit
   * cc-pvqz-jkfit
   * cc-pv5z-jkfit

Examples
--------

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : false,
     "geometry" : [
         {"atom" : "H", "xyz" : [ -0.22767998367, -0.82511994081,  -2.66609980874]; },
         {"atom" : "O", "xyz" : [  0.18572998668, -0.14718998944,  -3.25788976629]; },
         {"atom" : "H", "xyz" : [  0.03000999785,  0.71438994875,  -2.79590979943]; }
     ]
   },

   {
     "title" : "hf"
   }

   ]}

Example with mixed basis sets and density fitting basis sets:

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : "false",
     "geometry" : [
       { "atom" : "F",  "xyz" : [ -0.000000,     -0.000000,      2.720616]},
       { "atom" : "H",  "xyz" : [ -0.000000,     -0.000000,      0.305956],
                        "basis" : "cc-pvqz", "df_basis" : "cc-pvqz-jkfit" }
     ]
   },

   {
     "title" : "hf"
   }

   ]}

Example with running a calculation from a molden file using the keyword ``"basis" : "molden"``
and providing a value for ``"molden_file"``:

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "molden",
     "df_basis" : "svp-jkfit",
     "cartesian" : true,
     "molden_file" : "hf_write_mol_cart.molden"
   }

   ]}

(refer to Molden in :ref:`misc` for more details)

Example with external magnetic field using Gauge-invariant atomic orbitals (GIAO):

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "svp",
     "df_basis" : "svp-jkfit",
     "angstrom" : "false",
     "basis_type" : "giao",
     "tesla" : "false",
     "magnetic_field" : [  0.2000,   0.3000,  -0.1500   ],
     "geometry" : [
       { "atom" : "F",  "xyz" : [ -1.200000,      2.500000,      2.720616]},
       { "atom" : "H",  "xyz" : [ -1.200000,      2.500000,      0.305956]}
     ]
   },

   {
     "title" : "hf"
   }

   ]}

====================
Auxiliary basis sets
====================

The following MP2-fit basis sets are available in BAGEL. The basis set name can be used with the ``aux_basis`` keyword
in the method block (refer to :ref:`mp2` for more details).

* cc-pvdz-ri
* cc-pvtz-ri
* cc-pvqz-ri
* cc-pv5z-ri

Example
-------

An example using cc-pvdz-ri in MP2 calculation.

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
     "title" : "mp2",
     "aux_basis" : "cc-pvdz-ri",
     "frozen" : true
   }

   ]}

=========================================
Effective core potential (ECP) basis sets
=========================================
The following auxiliary basis sets are available in BAGEL library. The basis set name can be used with the ``basis`` keyword.

.. hlist::
   :columns: 3

   * ecp10mdf
   * ecp28mdf
   * ecp46mdf
   * ecp60mdf
   * ecp78mdf
   * def2-SVP-ecp
   * def2-SVP-2c-ecp
   * lanl2dz-ecp

.. note::
   User-defined ECP basis sets need to contain the keyword "ecp" in the names.
   Refer to `User defined basis sets`_ for more details.

Example
-------

Example for CuH2 using cc-pvtz basis set for H and lanl2dz-ecp for the heavy atom Cu

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "lanl2dz-ecp",
     "df_basis" : "svp-jkfit",
     "angstrom" : "true",
     "geometry" : [
       { "atom" : "Cu",  "xyz" : [  0.000000,      0.000000,      0.000000]},
       { "atom" :  "H",  "xyz" : [  0.000000,      0.000000,     -1.560000],
                         "basis" : "cc-pvtz"},
       { "atom" :  "H",  "xyz" : [  0.000000,      0.000000,      1.560000],
                         "basis" : "cc-pvtz"}
     ]
   },

   {
     "title" : "hf",
     "charge" : "-1"
   }

   ]}

========================
User defined basis sets
========================

The basis set file is in the following format

.. code-block:: javascript

 {
  "H" : [
    {
      "angular" : "s",
      "prim" : [5.4471780, 0.8245470],
      "cont" : [[0.1562850, 0.9046910]]
    }, {
      "angular" : "s",
      "prim" : [0.1831920],
      "cont" : [[1.0000000]]
    }
  ],
  "He" : [
    {
      "angular" : "s",
      "prim" : [13.6267000, 1.9993500],
      "cont" : [[0.1752300, 0.8934830]]
    }, {
      "angular" : "s",
      "prim" : [0.3829930],
      "cont" : [[1.0000000]]
    }
  ]
 }

The file is essentially one large array, the elements of which are further arrays, each corresponding to the basis set for a given element.
The basis set for associated with each element is then made up of further arrays, each of which  contains information specifying the properties
of a single basis function.

  * ``angular`` defines the kind of orbital (s,p,d,f...) .
  * ``prim`` is a array containing the exponents of the primitive orbitals from which the basis function is composed.
  * ``cont`` is an array containing the coefficients associated with each of these primitive orbitals.

The user can specify their own basis set using the above format, or use one of the predefined basis sets listed in `Basis sets`_.

.. note::
   Not all of the the basis sets are defined for all atoms;  an error message of form "No such node(X)", where X is the element, typically means that the relevant element was not found in the basis set file. Refer to the EMSL Basis set exchange library for more basis sets (https://bse.pnl.gov/bse/portal).

To use a user specified basis the explicit path to the basis set file must be specified in the basis set block.

Example
-------

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "/path/to/my/basis",
     "df_basis" : "/path/to/my/basis",
     "angstrom" : false,
     "geometry" : [
         {"atom" : "H", "xyz" : [ -0.22767998367, -0.82511994081,  -2.66609980874]; },
         {"atom" : "O", "xyz" : [  0.18572998668, -0.14718998944,  -3.25788976629]; },
         {"atom" : "H", "xyz" : [  0.03000999785,  0.71438994875,  -2.79590979943]; }
     ]
   },

   {
     "title" : "hf"
   }

   ]}

Other features
==============

===========
Dummy atoms
===========
Artificial point charges can be included in the calculation.
They introduce a user specified charge into the system, but  have no associated basis functions.
Introduction of such a charge is accomplished by inclusion of an additional line in the geometry block for an atom of  element "Q".

Example
-------

A dihydrogen molecule with a nearby dummy charge of +0.2. Note that the charge specified in the "hf" block does not include the charge associated with the dummy atom.

.. code-block:: javascript

   { "bagel" : [

   {
     "title" : "molecule",
     "basis" : "tzvpp",
     "df_basis" : "tzvpp-jkfit",
     "angstrom" : "true",
     "geometry" : [
       { "atom" :  "Q",  "xyz" : [  0.000000,   0.000000,   2.0000], "charge" : "0.2"},
       { "atom" :  "H",  "xyz" : [  0.000000,   0.000000,   0.7414]},
       { "atom" :  "H",  "xyz" : [  0.000000,   0.000000,   0.0000]}
     ]
   },

   {
     "title" : "hf"
   }

   ]}


from which one obtains

.. code-block:: javascript


  === RHF iteration (tzvpp) ===

               o Fock build                                  0.01
      0         -1.12552716          0.00743295           0.01
               o DIIS                                        0.00
               o Diag                                        0.00
               o Post process                                0.00
               o Fock build                                  0.01
      1         -1.12987462          0.00139213           0.01
               o DIIS                                        0.00
               o Diag                                        0.00
               o Post process                                0.00
               o Fock build                                  0.01
      2         -1.13008781          0.00009095           0.01
               o DIIS                                        0.00
               o Diag                                        0.00
               o Post process                                0.00
               o Fock build                                  0.01
      3         -1.13008889          0.00000614           0.01
               o DIIS                                        0.00
               o Diag                                        0.00
               o Post process                                0.00
               o Fock build                                  0.01
      4         -1.13008889          0.00000054           0.01
               o DIIS                                        0.00
               o Diag                                        0.00
               o Post process                                0.00
               o Fock build                                  0.01
      5         -1.13008889          0.00000007           0.01
               o DIIS                                        0.00
               o Diag                                        0.00
               o Post process                                0.00
               o Fock build                                  0.01
      6         -1.13008889          0.00000000           0.01

    * SCF iteration converged.

    * Permanent dipole moment:
           (    0.000000,    -0.000000,    -0.427736) a.u.

Updates to Molecule Specification
=================================

The examples above each provide a complete set of information for the molecule to be constructed from scratch.  
Alternatively, if a molecule block has already been provided, a new molecule input can be provided with a subset of 
the required inputs to be altered.  This is often useful when optimized orbitals from one calculation are to be 
used to provide an initial guess for a more difficult calculation.  

===========================
Useful Parameters to Change
===========================

.. topic:: ``geometry``

   | **Description**: Specify atoms and their Cartesian coordinates
   | **Note**:
   |    If atom positions have changed, new molecular orbitals will be generated by projecting the old orbitals into the space spanned by the new basis functions.  
        The atom positions and basis set cannot be changed simultaneously.   
        If units of Angstrom are to be used, the ``angstrom`` keyword must be supplied within the same ``molecule`` block.  

.. topic:: ``basis``

   | **Description**: Change the basis set
   | **Note**:
   |    The new basis set will be applied to all atoms, unless exceptions are specified in a new ``geometry`` block.   
        The new basis set must be identical to or larger than that used in any previous calculation.  
        Molecular orbitals will be projected into the expanded space.  

.. topic:: ``df_basis``

   | **Description**: Change the auxiliary basis set used for density fitting
   | **Note**:
        The new fitting basis set will be applied to all atoms, unless exceptions are specified in a new ``geometry`` block.   

.. topic:: ``magnetic_field``

   | **Description**: Change the external magnetic field
   | **Note**:
   |    Because giao basis orbitals are field-dependent, a projection will be performed to update the molecular orbitals.  
        If units of Tesla are to be used, the ``tesla`` keyword must be supplied within the same ``molecule`` block. 

.. topic:: ``basis_type``

   | **Description**: For relativistic calculations, this parameter can be used to convert from gaussian-type to giao-type orbitals.   
   | **Note**:
   |    When using the results of standard zero-field calculations to provide guess orbitals for a calculation with finite magnetic field, 
        this parameter is used to indicate a change to the giao basis set and prepare for magnetic field.  
        The atom positions, basis set, and magnetic field must not be changed simultaneously with this parameter.  
        (This means that to use a nonzero magnetic field, a second ``molecule`` block must be supplied.)  


References
==========

+-----------------------------------------------+----------------------------------------------------------------------------------+
|          Description of Reference             |                               Reference                                          |
+===============================================+==================================================================================+
| General text on electronic structure theory   | A\. Szabo and N. S. Ostlund,                                                     |
|                                               | *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory* |
|                                               | (McGraw-Hill, New York, 1989).                                                   |
+-----------------------------------------------+----------------------------------------------------------------------------------+
| Gauge invariant atomic orbitals               | R\. Ditchfield, Mol. Phys. **27**, 789 (1974).                                   |
+-----------------------------------------------+----------------------------------------------------------------------------------+
| Finite nuclear charge distribution and        | L\. Visscher and K. G. Dyall, At. Data Nucl. Data Tables **67**, 207 (1997).     |
| default exponents                             |                                                                                  |
+-----------------------------------------------+----------------------------------------------------------------------------------+

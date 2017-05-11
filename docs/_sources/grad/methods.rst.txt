.. _methods:

******************************
Description of input structure
******************************

Description
===========
One key feature that is heavily used in the nuclear gradients and relevant functionalities in Bagel is to use method array within the force (or optimize, or hessian) section. Since this array allow for quite a bit of flexibility in the input, several examples are given below to demonstrate its use.

Examples
========

In a force calculation, the methods can be nested in the followoing way. This is particularly useful for performing a CASPT2 gradient calculation where a CASSCF calculation must also be performed. The options for CASPT2 are specified in the smith block, while the CASSCF parameters (in this case nstate, nact, and nclosed) are specified outside of the smith block.

.. code-block:: javascript

  {
    "title" : "force",
    "method" : [ {
      "title" : "caspt2",
      "smith" : {
        "method" : "caspt2",
        "shift" : 0.2,
        "frozen" : true
      },
      "nstate" : 3,
      "nact" : 7,
      "nclosed" : 44
    } ]
  }

The same calculation could also be performed by setting up the input specifying both CASSCF and CASPT2 within in the method array but by first requesting a CASSCF calculation and subsequently asking for a CASPT2 calculation.

.. code-block:: javascript

  {
    "title" : "force",
    "method" : [
      {
        "title" : "casscf",
        "nstate" : 3,
        "nact" : 7,
        "nclosed" : 44
      },

      {
        "title" : "smith",
        "method" : "caspt2",
        "shift" : 0.2,
        "frozen" : true
      }
    ]
  }



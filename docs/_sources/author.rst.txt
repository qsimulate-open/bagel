.. _author:

********************
Author contributions
********************

Version 1.1
===========
**Jae Woo Park** has substantially improved the CASPT2 energy/gradient code. **Nils E. Strand** has implemented the spin-free DKH2 nuclear gradient code (available in v1.1.1).

Version 1.0
===========

**Toru Shiozaki** (TS) has designed the density fitting utilities and implemented as a parallel program. He also implemented most of the molecular integrals and the non-relativistic HF, MP2, CASSCF, CASPT2, NEVPT2, and MRCI programs. He is also responsible for the code for HF, MP2, CASSCF, and SA-CASSCF nuclear gradients and the *Z*-vector equations for HF and CASSCF reference functions. The determinantal FCI code has been developed by **TS** and **Shane M. Parker** (SMP). The parallel FCI is due to **TS**. All of the code for restricted active spaces is due to **SMP**. All of the FMM capabilities and ECP integrals were implemented by **Hai-Anh Le** (HAL). Douglas--Kroll--Hess integrals were implemented by **Yiqun Wang** (YW).

The code for state-specific CASPT2 nuclear energy gradients was developed by **Matthew K. MacLeod** (MKM) and **TS**. The contributions by **MKM** included the extension of SMITH3 to enable computation of the source terms in the *Z*-CASSCF equations. This code was later extended by **Bess Vlaisavljevich** (BV) to XMS-CASPT2 nuclear energy gradients and by **Jae Woo Park** (JWP) to include the computation of derivative couplings. Note that **JWP** has significantly improved the CASPT2 gradient algorithm, which is currently used in the released version. **BV** implemented the numerical differentiation programs for nuclear Hessian elements and vibrational analyses.

The Dirac–Hartree–Fock program was developed by **TS** with contributions from **Matthew S. Kelley**. The relativistic FCI was written by **TS**. The initial relativistic CASSCF program was implemented by **Jefferson E. Bates**. The program has since been replaced by a more stable second-order algorithm by **TS**, who also implemented the relativistic CASPT2 and MRCI programs. **Ryan D. Reynolds** (RDR) extended the relativistic CASPT2 and MRCI code to enable multi-state variants. **RDR** also developed a module that extracts EPR Hamiltonian parameters from the relativistic computations. The pilot implementation of nuclear gradients within the relativistic framework was written by **TS**.

**JWP** implemented all of the geometry optimizers (for equilibrium geometries, transition states,  and conical intersections) that are currently used, replacing inefficient code due to **TS**. **JWP** is also responsible for the interfaces to dynamics programs. All of the code for the ASD algorithms was developed by **SMP** (except orbital optimization part written by **Inkoo Kim**). The GIAO code was written by **RDR**.

**Peter J. Cherry**, **JWP**, **BV**, **HAL**, **RDR**, **Yeonjun Jeong**, **Jheng-Wei Li**, and **YW** wrote the manual at a group retreat in Wisconsin (May 1-3, 2017) while **TS** "supervised".

.. figure:: _static/group.jpg
   :width: 100 %
   :align: center


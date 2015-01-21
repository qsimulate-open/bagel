<<<<<<< HEAD
SUBDIRS = util input parallel math molecule integral io df opt grad mp2 wfn scf pscf pmp2 casscf nevpt2 smith prop dimer ciutil fci zfci zcasscf ras asd rel ks basis alglib london asdscf
AM_CXXFLAGS=-I$(top_srcdir)

INTLIBS = opt/libbagel_opt.la ks/libbagel_ks.la rel/libbagel_rel.la zcasscf/libbagel_zcasscf.la zfci/libbagel_zfci.la fci/libbagel_fci.la ras/libbagel_ras.la scf/libbagel_scf.la mp2/libbagel_mp2.la pscf/libbagel_pscf.la pmp2/libbagel_pmp2.la casscf/libbagel_casscf.la nevpt2/libbagel_nevpt2.la grad/libbagel_grad.la wfn/libbagel_wfn.la df/libbagel_df.la smith/libbagel_smith.la prop/libbagel_prop.la asd/libbagel_asd.la dimer/libbagel_dimer.la ciutil/libbagel_ciutil.la io/libbagel_io.la molecule/libbagel_molecule.la integral/libbagel_integral.la alglib/libbagel_alglib.la london/libbagel_london.la asdscf/libbagel_asdscf.la
=======
SUBDIRS = util molecule integral df opt grad pt2 wfn scf multi smith prop ci asd mat1e basis alglib
AM_CXXFLAGS=-I$(top_srcdir)

INTLIBS = opt/libbagel_opt.la multi/zcasscf/libbagel_zcasscf.la ci/zfci/libbagel_zfci.la ci/fci/libbagel_fci.la ci/ras/libbagel_ras.la scf/libbagel_scf.la multi/casscf/libbagel_casscf.la pt2/libbagel_pt2.la grad/libbagel_grad.la wfn/libbagel_wfn.la df/libbagel_df.la smith/libbagel_smith.la prop/libbagel_prop.la asd/libbagel_asd.la asd/dmrg/libbagel_asd_dmrg.la asd/dimer/libbagel_dimer.la asd/multisite/libbagel_multisite.la ci/ciutil/libbagel_ciutil.la mat1e/libbagel_mat1e.la util/io/libbagel_io.la molecule/libbagel_molecule.la integral/libbagel_integral.la alglib/libbagel_alglib.la
>>>>>>> origin/master

lib_LTLIBRARIES = libbagel.la
libbagel_la_SOURCES = static.cc
libbagel_la_LIBADD = $(INTLIBS)

bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc
BAGEL_LDADD = libbagel.la $(INTLIBS)

check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = libbagel.la $(INTLIBS)

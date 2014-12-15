SUBDIRS = util molecule integral df opt grad pt2 wfn scf multi smith prop ci asd rel basis alglib london
AM_CXXFLAGS=-I$(top_srcdir)

INTLIBS = opt/libbagel_opt.la rel/libbagel_rel.la multi/zcasscf/libbagel_zcasscf.la ci/zfci/libbagel_zfci.la ci/fci/libbagel_fci.la ci/ras/libbagel_ras.la scf/libbagel_scf.la multi/casscf/libbagel_casscf.la pt2/libbagel_pt2.la grad/libbagel_grad.la wfn/libbagel_wfn.la df/libbagel_df.la smith/libbagel_smith.la prop/libbagel_prop.la asd/libbagel_asd.la asd/dmrg/libbagel_asd_dmrg.la asd/dimer/libbagel_dimer.la asd/multisite/libbagel_multisite.la ci/ciutil/libbagel_ciutil.la util/io/libbagel_io.la molecule/libbagel_molecule.la integral/libbagel_integral.la alglib/libbagel_alglib.la london/libbagel_london.la

lib_LTLIBRARIES = libbagel.la
libbagel_la_SOURCES = static.cc
libbagel_la_LIBADD = $(INTLIBS)

bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc
BAGEL_LDADD = libbagel.la $(INTLIBS)

check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = libbagel.la $(INTLIBS)

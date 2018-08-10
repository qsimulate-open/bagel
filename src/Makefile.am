SUBDIRS = util molecule integral df opt grad pt2 wfn scf multi smith prop ci response asd mat1e basis
AM_CXXFLAGS=-I$(top_srcdir)

INTLIBS = opt/libbagel_opt.la multi/zcasscf/libbagel_zcasscf.la ci/zfci/libbagel_zfci.la ci/fci/libbagel_fci.la ci/ras/libbagel_ras.la response/libbagel_response.la scf/libbagel_scf.la multi/casscf/libbagel_casscf.la pt2/libbagel_pt2.la grad/libbagel_grad.la wfn/libbagel_wfn.la df/libbagel_df.la smith/libbagel_smith.la prop/libbagel_prop.la asd/libbagel_asd.la asd/dmrg/libbagel_asd_dmrg.la asd/orbital/libbagel_asd_orbital.la asd/dimer/libbagel_dimer.la asd/multisite/libbagel_multisite.la ci/ciutil/libbagel_ciutil.la mat1e/libbagel_mat1e.la integral/libbagel_integral.la molecule/libbagel_molecule.la util/io/libbagel_io.la util/libbagel_util.la util/parallel/libbagel_parallel.la util/math/libbagel_math.la util/input/libbagel_input.la

lib_LTLIBRARIES = libbagel.la
libbagel_la_SOURCES = run_bagel.cc static.cc
libbagel_la_LIBADD = $(INTLIBS)
include_HEADERS = bagel.h

bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc
BAGEL_LDADD = libbagel.la $(INTLIBS)

check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = libbagel.la $(INTLIBS)

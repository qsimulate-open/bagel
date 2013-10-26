SUBDIRS = util input parallel math molecule integral io df opt grad mp2 wfn scf pscf pmp2 input casscf smith prop dimer fci zfci zcasscf ras meh rel ks opt basis alglib
AM_CXXFLAGS=-I$(top_srcdir)

INTLIBS = opt/libbagel_opt.la ks/libbagel_ks.la rel/libbagel_rel.la zfci/libbagel_zfci.la fci/libbagel_fci.la ras/libbagel_ras.la scf/libbagel_scf.la mp2/libbagel_mp2.la pscf/libbagel_pscf.la pmp2/libbagel_pmp2.la zcasscf/libbagel_zcasscf.la casscf/libbagel_casscf.la grad/libbagel_grad.la wfn/libbagel_wfn.la df/libbagel_df.la smith/libbagel_smith.la prop/libbagel_prop.la meh/libbagel_meh.la dimer/libbagel_dimer.la io/libbagel_io.la molecule/libbagel_molecule.la integral/libbagel_integral.la alglib/libbagel_alglib.la

lib_LTLIBRARIES = libbagel.la
libbagel_la_SOURCES = static.cc
libbagel_la_LIBADD = $(INTLIBS)

bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc
BAGEL_LDADD = libbagel.la $(INTLIBS)

check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = libbagel.la $(INTLIBS)

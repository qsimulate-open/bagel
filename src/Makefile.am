INTLIBS = parallel/libbagel_parallel.la scf/libbagel_scf.la grad/libbagel_grad.la rysint/libbagel_rysint.la osint/libbagel_osint.la pscf/libbagel_pscf.la pmp2/libbagel_pmp2.la util/libbagel_util.la casscf/libbagel_casscf.la wfn/libbagel_wfn.la df/libbagel_df.la smith/libbagel_smith.la mp2/libbagel_mp2.la prop/libbagel_prop.la fci/libbagel_fci.la dimer/libbagel_dimer.la io/libbagel_io.la rel/libbagel_rel.la
# opt/libbagel_opt.a

SUBDIRS = parallel io df opt grad mp2 wfn rysint osint scf pscf pmp2 util casscf smith prop dimer fci rel
bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc
BAGEL_LDADD = $(INTLIBS) $(FLIBS)
AM_CXXFLAGS=-I$(top_srcdir)


check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = $(INTLIBS) $(FLIBS)

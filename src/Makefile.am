INTLIBS = parallel/libparallel.la scf/libscf.la grad/libgrad.la rysint/librysint.la osint/libosint.la pscf/libpscf.la pmp2/libpmp2.la util/libutil.la casscf/libcasscf.la wfn/libwfn.la df/libdf.la smith/libsmith.la mp2/libmp2.la prop/libprop.la fci/libfci.la dimer/libdimer.la io/libio.la rel/librel.la
# opt/libopt.a

SUBDIRS = parallel io df opt grad mp2 wfn rysint osint scf pscf pmp2 util casscf smith prop dimer fci rel
bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc
BAGEL_LDADD = $(INTLIBS) $(FLIBS)
AM_CXXFLAGS=-I$(top_srcdir)


check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = $(INTLIBS) $(FLIBS)

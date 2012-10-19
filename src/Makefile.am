INTLIBS = scf/libscf.a grad/libgrad.a rysint/librysint.a osint/libosint.a slater/libslater.a pscf/libpscf.a pmp2/libpmp2.a util/libutil.a casscf/libcasscf.a wfn/libwfn.a df/libdf.a smith/libsmith.a mp2/libmp2.a prop/libprop.a fci/libfci.a dimer/libdimer.a io/libio.a rel/librel.a -L./slater -lslater -L./grad -lgrad -L./osint -losint
# opt/libopt.a

SUBDIRS = df opt grad mp2 wfn rysint osint scf slater pscf pmp2 util casscf smith prop dimer fci rel io
bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc
BAGEL_LDADD = $(INTLIBS) $(FLIBS)
AM_CXXFLAGS=-I$(top_srcdir)


check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = $(INTLIBS) $(FLIBS)

INTLIBS = scf/libscf.a grad/libgrad.a rysint/librysint.a osint/libosint.a slater/libslater.a pscf/libpscf.a pmp2/libpmp2.a util/libutil.a fci/libfci.a casscf/libcasscf.a wfn/libwfn.a df/libdf.a smith/libsmith.a mp2/libmp2.a prop/libprop.a molden/libmolden.a dimer/libdimer.a -L./slater -lslater 
# opt/libopt.a

SUBDIRS = df opt grad mp2 wfn rysint osint scf slater pscf pmp2 util fci casscf smith prop molden dimer
bin_PROGRAMS = PolyR12
PolyR12_SOURCES = main.cc
PolyR12_LDADD = $(INTLIBS) -lboost_regex $(FLIBS)
AM_CXXFLAGS=-I$(top_srcdir)


check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = $(INTLIBS) -lboost_unit_test_framework -lboost_regex $(FLIBS) 

SUBDIRS = df grad mp2 wfn rysint osint scf slater pscf pmp2 util fci casscf smith
bin_PROGRAMS = PolyR12
PolyR12_SOURCES = main.cc
PolyR12_LDADD = scf/libscf.a grad/libgrad.a rysint/librysint.a osint/libosint.a slater/libslater.a pscf/libpscf.a pmp2/libpmp2.a util/libutil.a fci/libfci.a casscf/libcasscf.a wfn/libwfn.a df/libdf.a smith/libsmith.a mp2/libmp2.a -L./slater -lslater -lboost_regex -lblas -llapack $(FLIBS)
AM_CXXFLAGS=-I$(top_srcdir)

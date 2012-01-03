SUBDIRS = wfn rysint osint scf slater pscf pmp2 util fci casscf df
bin_PROGRAMS = PolyR12
PolyR12_SOURCES = main.cc
PolyR12_LDADD = scf/libscf.a rysint/librysint.a osint/libosint.a slater/libslater.a pscf/libpscf.a pmp2/libpmp2.a util/libutil.a fci/libfci.a casscf/libcasscf.a wfn/libwfn.a df/libdf.a -L./slater -lslater -llapack $(FLIBS)
AM_CXXFLAGS=-I$(top_srcdir)

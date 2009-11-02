SUBDIRS = rysint osint scf slater pscf pmp2 util
bin_PROGRAMS = PolyR12
PolyR12_SOURCES = main.cc
PolyR12_LDADD = scf/libscf.a rysint/librysint.a osint/libosint.a slater/libslater.a pscf/libpscf.a pmp2/libpmp2.a util/libutil.a $(FLIBS) -lblas -llapack


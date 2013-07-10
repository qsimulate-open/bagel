INTLIBS = opt/libbagel_opt.la ks/libbagel_ks.la rel/libbagel_rel.la fci/libbagel_fci.la zfci/libbagel_zfci.la scf/libbagel_scf.la mp2/libbagel_mp2.la pscf/libbagel_pscf.la pmp2/libbagel_pmp2.la casscf/libbagel_casscf.la grad/libbagel_grad.la wfn/libbagel_wfn.la df/libbagel_df.la smith/libbagel_smith.la prop/libbagel_prop.la meh/libbagel_meh.la dimer/libbagel_dimer.la io/libbagel_io.la osint/libbagel_osint.la rysint/libbagel_rysint.la molecule/libbagel_molecule.la input/libbagel_input.la math/libbagel_math.la parallel/libbagel_parallel.la util/libbagel_util.la

SUBDIRS = util parallel math io df opt grad mp2 wfn rysint osint scf pscf pmp2 molecule input casscf smith prop dimer fci zfci meh rel ks opt
bin_PROGRAMS = BAGEL
BAGEL_SOURCES = main.cc static.cc
BAGEL_LDADD = $(INTLIBS)
AM_CXXFLAGS=-I$(top_srcdir)


check_PROGRAMS = TestSuite
TestSuite_SOURCES = test_main.cc
TestSuite_LDADD = $(INTLIBS)

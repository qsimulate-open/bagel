//
// Author: Toru Shiozaki
// Date  : April 2009
//
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <memory>
#if 0
#include <src/pscf/pgeometry.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/pscf.h>
#include <src/pscf/pscf_disk.h>
#include <src/pmp2/pmp2.h>
#endif
#include <src/scf/geometry.h>
#include <src/scf/scf.h>
#include <src/wfn/reference.h>
#include <src/fci/fci.h>
#include <src/casscf/superci.h>
#include <src/casscf/werner.h>
#include <src/global.h>
#include <src/stackmem.h>

#include <src/util/input.h>

using namespace std;

StackMem* stack;

// debug
extern void a();

int main(int argc, char** argv) {
a();
  // openmp is broken now due to the use of stack.
  // What we need is a proper thread model.
  #ifdef _OPENMP
  assert(false); // trap
  #endif

  try {
    print_header();

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw runtime_error("no input file provided");
    }
    const string input = argv[1];

    shared_ptr<InputData> idata(new InputData(input));

    const bool fci_card = idata->exist("fci"); 
    const bool casscf_card = idata->exist("casscf");

    shared_ptr<Geometry> geom(new Geometry(idata));
    list<pair<string, multimap<string, string> > > keys = idata->data();

    bool scf_done = false;
    bool casscf_done = false;
    shared_ptr<SCF_base> scf;
    shared_ptr<CASSCF> casscf;
    shared_ptr<FCI> fci;
    shared_ptr<Reference> ref;

    for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
      const string method = iter->first;

      if (method == "hf") {

        shared_ptr<SCF<0> > scf_(new SCF<0>(iter->second, geom)); scf = scf_;
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-hf") {

        if (!geom->df()) throw runtime_error("It seems that DF basis was not specified in Geometry");
        shared_ptr<SCF<1> > scf_(new SCF<1>(iter->second, geom)); scf = scf_;
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "casscf") {
        if (!ref) throw runtime_error("CASSCF needs a reference");

        string algorithm = read_input<string>(iter->second, "algorithm", ""); 
        if (algorithm == "superci" || algorithm == "") {
          shared_ptr<CASSCF> casscf_(new SuperCI(iter->second, geom, ref)); casscf = casscf_;
          casscf->compute();
          ref = casscf->conv_to_ref();
        } else if (algorithm == "werner" || algorithm == "knowles") {
          shared_ptr<CASSCF> werner(new WernerKnowles(iter->second, geom, ref));
          werner->compute();
          ref = werner->conv_to_ref();
        } else {
          throw runtime_error("unknown CASSCF algorithm specified.");
        }

      } else if (method == "fci") {
        if (!ref) throw runtime_error("FCI needs a reference");

        shared_ptr<FCI> fci_(new FCI(iter->second, geom, ref)); fci = fci_;
        fci->compute();

      }
    }
    print_footer();

    // end of the main file

#if 0
    const bool use_hy2 = (bool)count_string(input, "HY2");
    const int depth_basis = count_string(input, "Basis");
    const bool periodic = (bool)count_string(input, "Periodic");
    const bool domp2 = (bool)count_string(input, "MP2");
    typedef std::shared_ptr<Geometry> RefGeom;
    typedef std::shared_ptr<PGeometry> RefPGeom;
    typedef std::shared_ptr<PSCF_DISK> RefPSCF_DISK;

    RefPSCF_DISK pscf;

    if (depth_basis == 2) {
      if (!periodic) {
        throw runtime_error("Haven't supported projection in molecular cases.");
      } else {
        RefPGeom pgeom(new PGeometry(input, 0));
        RefPGeom pgeom2(new PGeometry(input, 1));
        RefPSCF_DISK pscf2(new PSCF_DISK(pgeom2));
        pscf2->compute();

        // DO PROJECTION!
        shared_ptr<PMatrix1e> density_new;
        {
          RefPGeom uniongeom(new PGeometry(*pgeom2));
          uniongeom->merge_obs_cabs();
          POverlap union_overlap_r(uniongeom);
          PMatrix1e union_overlap(union_overlap_r.ft());
          const int nold = pgeom2->nbasis();
          const int nnew = pgeom2->ncabs();
          PMatrix1e s_new_old(union_overlap.split(nold, nnew).first, make_pair(nold, nold+nnew));
          PMatrix1e s_new_new(union_overlap.split(nold, nnew).second, make_pair(nold, nold+nnew));
          PMatrix1e s_new_new_inv = *s_new_new.inverse();
          PMatrix1e density_old = pscf2->coeff()->form_density_rhf(false);
          PMatrix1e transform = s_new_old * s_new_new_inv;
          PMatrix1e density_new_k = transform % density_old * transform;

          density_new_k.set_geom(pgeom);
          shared_ptr<PMatrix1e> density_new_i(new PMatrix1e(density_new_k.bft()));
          density_new = density_new_i;
        }
        RefPSCF_DISK tmp(new PSCF_DISK(pgeom, density_new));
        pscf = tmp;
        pscf->compute();
      }
    } else if (depth_basis == 1) {
      if (!periodic) {
        RefGeom geom(new Geometry(input, 0));
        if (!fci_card && !casscf_card) {
          SCF scf(geom);
          scf.compute();
        } else if (casscf_card) {
          SuperCI casscf(geom);
          casscf.compute();
        } else {
          FCI fci(geom);
          fci.compute();
        }
      } else {
        RefPGeom pgeom(new PGeometry(input, 0));
        RefPSCF_DISK tmp(new PSCF_DISK(pgeom));
        pscf = tmp;
        pscf->compute();
      }
    } else {
      throw runtime_error("Haven't supported triple basis sets..");
    }

    if (periodic && domp2) {
      PMP2 pmp2(pscf->geom(), pscf->coeff(), pscf->eig(), pscf->ao_eri(), use_hy2);
      pmp2.compute();
    }
    print_footer();
#endif


  } catch (const std::exception &e) {
    cout << "  ERROR: EXCEPTION RAISED:" << e.what() << endl;
    throw;
  } catch (...) {
    throw;
  }

  return 0;
}


/*
This file computes the ERI over London orbitals.
--Currently, it can only handle s-orbitals, but I will probably expand it to include angular momentum.
--A basis set is defined in an external input file passed as an argument to the main function
--We also define a set of molecular orbitals, with coefficients given in
  "orbital1", "orbital2", and so on, then it computes the full ERI over those two molecular orbitals by expanding
  them out in the basis and adding together the needed terms.
--Complex Gaussian quadrature is used, so the compiler must be able to see the interpolation files in BAGEL/src/integral/comprys
*/

#include <iostream>
#include <iomanip>
#include <cassert>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <utility>
#include "polynomial.h"
#include "ericompute.h"

using namespace std;
using namespace ryan;


int main (int argc, char*argv[]) {

  cout << setprecision(14) << endl;
  if (argc==1) {
    cout << "Please specify an input file." << endl;
    return 1;
  }

  ifstream ifs (argv[1]);
  string content ((istreambuf_iterator<char>(ifs)),(istreambuf_iterator<char>()));
  stringstream ss (content);

  int nbasis_contracted = 0;
  int natom = 0;
  char calculation = '0';
  bool normalize_basis = 0;
  bool scale_input = 0;
  bool orthogonalize = 0;
  bool first = 1;

  vector<double> field = {};
  vector<double> positions = {};
  vector<int> angular = {};
  vector<double> exponents = {};
  vector<double> contraction_coefficients = {};
  vector<int> nprimitive = {};
  vector<double> nuclear_positions = {};
  vector<int> nuclear_charge = {};
  vector<complex<double>> orbital1 = {};
  vector<complex<double>> orbital2 = {};
  vector<complex<double>> orbital3 = {};
  vector<complex<double>> orbital4 = {};

  while(ss) {
    if (first) {
      string stuff;
      getline( ss, stuff, ':');
      ss >> calculation;
      getline( ss, stuff, ':');
      ss >> nbasis_contracted;
      getline( ss, stuff, ':');
      ss >> normalize_basis;
      getline( ss, stuff, ':');
      ss >> scale_input;
      getline( ss, stuff, ':');
      ss >> orthogonalize;

      field.resize(3);
      nprimitive.resize(nbasis_contracted);
      orbital1.resize(nbasis_contracted);
      orbital2.resize(nbasis_contracted);
      orbital3.resize(nbasis_contracted);
      orbital4.resize(nbasis_contracted);

      for (int k=0; k!=3; k++) {
        getline( ss, stuff, ':');
        ss >> field[k];
      }
      for (int i=0; i!=nbasis_contracted; i++) {
        getline( ss, stuff, ':');
        ss >> orbital1[i];
      }
      for (int i=0; i!=nbasis_contracted; i++) {
        getline( ss, stuff, ':');
        ss >> orbital2[i];
      }
      if (calculation == 'E') {
        for (int i=0; i!=nbasis_contracted; i++) {
          getline( ss, stuff, ':');
          ss >> orbital3[i];
        }
        for (int i=0; i!=nbasis_contracted; i++) {
          getline( ss, stuff, ':');
          ss >> orbital4[i];
        }
      }
      if (calculation == 'N') {
        nuclear_charge.resize(natom);
        nuclear_positions.resize(3*natom);
        for (int i=0; i!=natom; i++) {
          getline( ss, stuff, ':');
          ss >> nuclear_charge[i];
          for (int j=0; j!=3; j++) {
            getline( ss, stuff, ':');
            ss >> nuclear_positions[3*i+j];
          }
        }
      }
      for (int i=0; i!=nbasis_contracted; i++) {
        double coords[3];
        int angmom[3];
        double temp;
          for (int k=0; k!=3; k++) {
            getline( ss, stuff, ':');
            ss >> coords[k];
          }
          for (int k=0; k!=3; k++) {
            getline( ss, stuff, ':');
            ss >> angmom[k];
          }
        getline( ss, stuff, ':');
        ss >> nprimitive[i];
        for (int j=0; j!=nprimitive[i]; j++) {
          getline( ss, stuff, ':');
          ss >> temp;
          exponents.push_back(temp);
          for (int k=0; k!=3; k++) {
            positions.push_back(coords[k]);
            angular.push_back(angmom[k]);
          }
        }
        for (int j=0; j!=nprimitive[i]; j++) {
          getline( ss, stuff, ':');
          ss >> temp;
          contraction_coefficients.push_back(temp);
        }
      }
    first = 0;
    }
    ss.ignore( numeric_limits<streamsize>::max(), '\n' );
  }

  const int nbasis = exponents.size();
  if (calculation == 'E') cout << "Full ERI summation over four molecular orbitals specified below,";
  else if (calculation == 'N') cout << "Full NAI summation over two molecular orbitals and " << natom << " nuclei specified below,";
  else throw runtime_error ("Desired calculation not recognized.  Enter 'E' for ERI or 'N' for NAI.");
  if (scale_input) cout << "with coefficients scaled to one, ";
  else cout << "with coefficients taken as-is, ";
  cout << endl << "using a basis of " << nbasis;
  if (orthogonalize) cout << " orthogonal combinations of ";
  else cout << " atom-centered (nonorthogonal) ";
  if (normalize_basis) cout << "normalized";
  else cout << "unnormalized";
  cout << " London orbitals." << endl;

  //if (nbasis < 4) throw runtime_error ("Need at least four basis orbitals, but some can be unused if you set their coefficients to zero in the MO.");
#if 1
  // to print out input parameters
  cout << "  nbasis = " << nbasis << endl;
  cout << "  normalize_basis = " << normalize_basis << endl;
  cout << "  scale_input = " << scale_input << endl;
  cout << "  orthogonalize = " << orthogonalize << endl;
  cout << "  field = " << field[0] << ", " << field[1] << ", " << field[2] << endl;

  cout << "  molecular orbital 1 coefficients = ";
  for (int j=0; j!=nbasis_contracted; j++) {
    cout << orbital1[j];
    if (j!=nbasis_contracted-1) cout << ", ";
    else cout << endl;
  }
  cout << "  molecular orbital 2 coefficients = ";
  for (int j=0; j!=nbasis_contracted; j++) {
    cout << orbital2[j];
    if (j!=nbasis_contracted-1) cout << ", ";
    else cout << endl;
  }
  if (calculation == 'E') {
  cout << "  molecular orbital 3 coefficients = ";
    for (int j=0; j!=nbasis_contracted; j++) {
      cout << orbital3[j];
      if (j!=nbasis_contracted-1) cout << ", ";
      else cout << endl;
    }
    cout << "  molecular orbital 4 coefficients = ";
    for (int j=0; j!=nbasis_contracted; j++) {
      cout << orbital4[j];
      if (j!=nbasis_contracted-1) cout << ", ";
      else cout << endl;
    }
  }
  for (int i = 0; i!=nbasis; i++) {
    cout << "  Atomic orbital " << i << endl;
    cout << "    exponent = " << exponents[i] << endl;
    cout << "    positions = " << positions[3*i+0] << ", " << positions[3*i+1] << ", " << positions[3*i+2] << endl;
    cout << "    angular = " << angular[3*i+0] << ", " << angular[3*i+1] << ", " << angular[3*i+2] << endl;
  }
#endif
cout << endl;

  if (calculation == 'E') {
    complex<double> FULL_ERI = compute_eri (nbasis_contracted, normalize_basis, scale_input, orthogonalize, field,
                               positions, angular, exponents, contraction_coefficients, nprimitive, orbital1, orbital2, orbital3, orbital4);
    cout << "Final result = " << FULL_ERI << endl;
  } else if (calculation == 'N') {
    cout << "NAI not coded yet." << endl;
  }

  return 0;

}

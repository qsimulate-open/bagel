//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: petite.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __SRC_MOLECULE_PETITE_H
#define __SRC_MOLECULE_PETITE_H

#include <tuple>
#include <src/molecule/atom.h>
#include <src/util/serialization.h>

namespace bagel {

class Petite {
  protected:
    int natom_;
    int nshell_;
    int nirrep_;
    int nsymop_;
    std::string sym_;

    std::vector<std::vector<double>> symop_;

    std::vector<std::vector<int>> sym_atommap_;
    std::vector<std::vector<int>> sym_shellmap_;
    std::vector<int> p1_;
    std::vector<int> lambda_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & natom_ & nshell_ & nirrep_ & nsymop_ & sym_ & symop_
         & sym_atommap_ & sym_shellmap_ & p1_ & lambda_;
    }

  public:
    Petite() { }
    Petite(const std::vector<std::shared_ptr<const Atom>>&, const std::string);
    ~Petite();

    std::vector<double> symop(const int i) const { return symop_[i]; };

    std::vector<int> sym_shellmap(const int i) const { return sym_shellmap_[i]; };
    std::vector<int> sym_atommap(const int i) const { return sym_atommap_[i]; };

    int nirrep() const { return nirrep_; };
    int nsymop() const { return nsymop_; };

    bool in_p1(int i) const { return (nirrep_ == 1 || p1_[i]); };
    bool in_p2(int ij) const { return (nirrep_ == 1 || lambda_[ij]); };
    int  in_p4(int ij, int kl, int i, int j, int k, int l) const {
      if (nirrep_ == 1) return 1;
      const int shift = sizeof(int) * 4;
      const int ijkl = ij < kl ? (ij << shift) + kl : (kl << shift) + ij;
      int nijkl = 1;
      for (int iop = 1; iop != nsymop_; ++iop) {
        const int gi = sym_shellmap_[i][iop];
        const int gj = sym_shellmap_[j][iop];
        const int gk = sym_shellmap_[k][iop];
        const int gl = sym_shellmap_[l][iop];
        const int gij = std::min(gi, gj) * nshell_ + std::max(gi, gj);
        const int gkl = std::min(gk, gl) * nshell_ + std::max(gk, gl);
        const int gijkl =  gij < gkl ? (gij << shift) + gkl : (gkl << shift) + gij;

        if (gijkl < ijkl) return 0;
        else if (gijkl == ijkl) ++nijkl;
      }
      return nsymop_ / nijkl;
    };

};


class Symmetry {
  protected:
    // vector of 9 doubles (xyz * xyz)
    std::vector<std::vector<double>> symop_;
  public:
    Symmetry() {};
    ~Symmetry() {};

    virtual int nirrep() const = 0;
    std::vector<std::vector<double>> symop() const { return symop_; };

};


class SymCs : public Symmetry {
  protected:

  public:
    SymCs() : Symmetry() {
      std::vector<double> E, sigmaxy;
      double Ea[9] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
      E.insert(E.begin(), Ea, Ea + 9);
      symop_.push_back(E);

      double sxya[9] = { 1.0,  0.0,  0.0,
                         0.0,  1.0,  0.0,
                         0.0,  0.0, -1.0};
      sigmaxy.insert(sigmaxy.begin(), sxya, sxya + 9);
      symop_.push_back(sigmaxy);

    };

    ~SymCs() {};
    int nirrep() const { return 2; };
};


class SymCi : public Symmetry {
  protected:

  public:
    SymCi() : Symmetry() {
      std::vector<double> E, I;
      double Ea[9] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
      E.insert(E.begin(), Ea, Ea + 9);
      symop_.push_back(E);

      double Ia[9] = {-1.0,  0.0,  0.0,
                       0.0, -1.0,  0.0,
                       0.0,  0.0, -1.0};
      I.insert(I.begin(), Ia, Ia + 9);
      symop_.push_back(I);

    };

    ~SymCi() {};
    int nirrep() const { return 2; };
};


class SymC2 : public Symmetry {
  protected:

  public:
    SymC2() : Symmetry() {
      std::vector<double> E, C2;
      double Ea[9] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
      E.insert(E.begin(), Ea, Ea + 9);
      symop_.push_back(E);

      double c2a[9] = {-1.0,  0.0,  0.0,
                        0.0, -1.0,  0.0,
                        0.0,  0.0,  1.0};
      C2.insert(C2.begin(), c2a, c2a + 9);
      symop_.push_back(C2);

    };

    ~SymC2() {};
    int nirrep() const { return 2; };
};


class SymC2h : public Symmetry {
  protected:

  public:
    SymC2h() : Symmetry() {
      std::vector<double> E, C2, I, sigmaxy;
      double Ea[9] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
      E.insert(E.begin(), Ea, Ea + 9);
      symop_.push_back(E);

      double C2a[9] = {-1.0,  0.0,  0.0,
                        0.0, -1.0,  0.0,
                        0.0,  0.0,  1.0};
      C2.insert(C2.begin(), C2a, C2a + 9);
      symop_.push_back(C2);

      double Ia[9] = {-1.0,  0.0,  0.0,
                       0.0, -1.0,  0.0,
                       0.0,  0.0, -1.0};
      I.insert(I.begin(), Ia, Ia + 9);
      symop_.push_back(I);

      double sza[9] = { 1.0,  0.0,  0.0,
                        0.0,  1.0,  0.0,
                        0.0,  0.0, -1.0};
      sigmaxy.insert(sigmaxy.begin(), sza, sza + 9);
      symop_.push_back(sigmaxy);
    };

    ~SymC2h() {};
    int nirrep() const { return 4; };
};


class SymD2 : public Symmetry {
  protected:

  public:
    SymD2() : Symmetry() {
      std::vector<double> E, C2x, C2y, C2z;
      double Ea[9] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
      E.insert(E.begin(), Ea, Ea + 9);
      symop_.push_back(E);

      double C2xa[9] = { 1.0,  0.0,  0.0,
                         0.0, -1.0,  0.0,
                         0.0,  0.0, -1.0};
      C2x.insert(C2x.begin(), C2xa, C2xa + 9);
      symop_.push_back(C2x);

      double C2ya[9] = {-1.0,  0.0,  0.0,
                         0.0,  1.0,  0.0,
                         0.0,  0.0, -1.0};
      C2y.insert(C2y.begin(), C2ya, C2ya + 9);
      symop_.push_back(C2y);

      double C2za[9] = {-1.0,  0.0,  0.0,
                         0.0, -1.0,  0.0,
                         0.0,  0.0,  1.0};
      C2z.insert(C2z.begin(), C2za, C2za + 9);
      symop_.push_back(C2z);

    };

    ~SymD2() {};
    int nirrep() const { return 4; };
};


class SymC2v : public Symmetry {
  protected:

  public:
    SymC2v() : Symmetry() {
      std::vector<double> E, C2, sigmaxz, sigmayz;
      double Ea[9] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
      E.insert(E.begin(), Ea, Ea + 9);
      symop_.push_back(E);

      double C2a[9] = {-1.0,  0.0,  0.0,
                        0.0, -1.0,  0.0,
                        0.0,  0.0,  1.0};
      C2.insert(C2.begin(), C2a, C2a + 9);
      symop_.push_back(C2);

      double sxa[9] = { 1.0,  0.0,  0.0,
                        0.0, -1.0,  0.0,
                        0.0,  0.0,  1.0};
      sigmaxz.insert(sigmaxz.begin(), sxa, sxa + 9);
      symop_.push_back(sigmaxz);

      double sya[9] = {-1.0,  0.0,  0.0,
                        0.0,  1.0,  0.0,
                        0.0,  0.0,  1.0};
      sigmayz.insert(sigmayz.begin(), sya, sya + 9);
      symop_.push_back(sigmayz);
    };

    ~SymC2v() {};
    int nirrep() const { return 4; };
};

class SymD2h : public Symmetry {
  protected:

  public:
    SymD2h() : Symmetry() {
      std::vector<double> E, C2z, C2y, C2x, I, sigmaxy, sigmayz, sigmazx;
      double Ea[9] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
      E.insert(E.begin(), Ea, Ea + 9);
      symop_.push_back(E);

      double C2za[9] = {-1.0,  0.0,  0.0,
                         0.0, -1.0,  0.0,
                         0.0,  0.0,  1.0};
      C2z.insert(C2z.begin(), C2za, C2za + 9);
      symop_.push_back(C2z);

      double C2ya[9] = {-1.0,  0.0,  0.0,
                         0.0,  1.0,  0.0,
                         0.0,  0.0, -1.0};
      C2y.insert(C2y.begin(), C2ya, C2ya + 9);
      symop_.push_back(C2y);

      double C2xa[9] = { 1.0,  0.0,  0.0,
                         0.0, -1.0,  0.0,
                         0.0,  0.0, -1.0};
      C2x.insert(C2x.begin(), C2xa, C2xa + 9);
      symop_.push_back(C2x);

      double Ia[9]  = {-1.0,  0.0,  0.0,
                        0.0, -1.0,  0.0,
                        0.0,  0.0, -1.0};
      I.insert(I.begin(), Ia, Ia + 9);
      symop_.push_back(I);

      double sxya[9] = { 1.0,  0.0,  0.0,
                         0.0,  1.0,  0.0,
                         0.0,  0.0, -1.0};
      sigmaxy.insert(sigmaxy.begin(), sxya, sxya + 9);
      symop_.push_back(sigmaxy);

      double syza[9] = {-1.0,  0.0,  0.0,
                         0.0,  1.0,  0.0,
                         0.0,  0.0,  1.0};
      sigmayz.insert(sigmayz.begin(), syza, syza + 9);
      symop_.push_back(sigmayz);

      double szxa[9] = { 1.0,  0.0,  0.0,
                         0.0, -1.0,  0.0,
                         0.0,  0.0,  1.0};
      sigmazx.insert(sigmazx.begin(), szxa, szxa + 9);
      symop_.push_back(sigmazx);

    };

    ~SymD2h() {};
    int nirrep() const { return 8; };
};

}

#endif



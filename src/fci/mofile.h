//
// author : Toru Shiozaki
//

#ifndef __NEWINT_FCI_MOFILE_H
#define __NEWINT_FCI_MOFILE_H

#include <fstream>
#include <string>
#include <memory>
#include <src/util/filename.h>
#include <src/scf/coeff.h>
#include <src/scf/geometry.h>

class MOFile {
  protected:
    const std::shared_ptr<Geometry> geom_;
    std::shared_ptr<std::fstream> file_;
    const std::shared_ptr<Coeff> coeff_; 
    size_t sizeij_;
    long filesize_;
    std::string filename_;
    std::vector<std::shared_ptr<Shell> > basis_;
    std::vector<int> offset_;

    std::vector<double> mo1e_;
    std::vector<double> mo2e_;

  public:
    MOFile(const std::shared_ptr<Geometry>, const std::shared_ptr<Coeff>);
    ~MOFile();

    void create_Jiiii(const int, const int);
    double mo1e(const size_t i) const { return mo1e_[i]; };
    double mo2e(const size_t i, const size_t j) const { return mo2e_[i+j*sizeij_]; };

};

#endif

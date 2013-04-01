//
//  dynapro_solver.h
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#ifndef lattice_binding_dynapro_solver_h
#define lattice_binding_dynapro_solver_h

#include <iostream>
#include <vector>
#include "solver.h"

namespace lattice {
  class DynaProSolver : public Solver {
  private:
    std::vector<double> Z_;
    std::vector<bool> Z_cached_;
    std::vector<std::vector<double>> dZdK_;
    std::vector<std::vector<bool>> dZdK_cached_;
    std::vector<std::vector<double>> c_;
    
    double Z(long n);
    double dZdK(long n, unsigned short g);
    
  public:
    DynaProSolver(const std::vector<Particle>& particles,
                  unsigned long N) : Solver(particles, N) {
      Z_ = std::vector<double>(N+1);
      Z_cached_ = std::vector<bool>(N+1);
      dZdK_ = std::vector<std::vector<double>>(N+1);
      dZdK_cached_ = std::vector<std::vector<bool>>(N+1);
      for (int i = 0; i <= N; i++) {
        Z_cached_[i] = false;
        dZdK_[i] = std::vector<double>(particles.size());
        dZdK_cached_[i] = std::vector<bool>(particles.size());
        for (int j = 0; j < particles.size(); j++) {
          dZdK_cached_[i][j] = false;
        }
      }
      c_ = std::vector<std::vector<double>>(N);
      for (int n = 0; n < N; n++) {
        c_[n] = std::vector<double>(particles.size());
      }
    }
    
    virtual double Z() override {
      double pf = Z(N_);
      std::cout << "Partition function Z = " << pf << std::endl;
      return pf;
    }
    
    virtual double c(unsigned long n, unsigned short g) override {      
      return dZdK(n, g) * K(n, g) / Z();
    }
  };
}

#endif

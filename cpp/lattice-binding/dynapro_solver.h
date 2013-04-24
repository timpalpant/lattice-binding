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
    std::vector<std::vector<double>> c_;
    
    double Z(long n);
    
  public:
    DynaProSolver(const std::vector<Particle>& particles,
                  unsigned long N) : Solver(particles, N) {
      Z_ = std::vector<double>(N+1);
      c_ = std::vector<std::vector<double>>(N);
      for (int n = 0; n < N; n++) {
        c_[n] = std::vector<double>(particles.size());
      }
    }
    
    virtual void solve() override;
    
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

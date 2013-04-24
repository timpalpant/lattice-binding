//
//  transfer_matrix_solver.h
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#ifndef lattice_binding_transfer_matrix_solver_h
#define lattice_binding_transfer_matrix_solver_h

#include <vector>
#include "Eigen/Sparse"
#include "solver.h"
#include "particle.h"

namespace lattice {
  class TransferMatrixSolver : public Solver {
  private:
    double Z_;
    std::vector<std::vector<double>> c_;
    
    void Q(Eigen::SparseMatrix<double>& Q, unsigned long n);
    void dQdK(Eigen::SparseMatrix<double>& dQdK, unsigned long n, unsigned short g);
    
    unsigned int num_states() const {
      unsigned int count = max_V_ + 3;
      for (const Particle& p : particles_) {
        count += p.size() + p.max_interaction_distance();
      }
      return count;
    }
    
  public:
    TransferMatrixSolver(const std::vector<Particle>& particles,
                         std::size_t N) : Solver(particles, N) { }
    
    virtual void solve() override;
    
    virtual double Z() override {
      if (Z_ == 0) {
        solve();
      }
      
      return Z_;
    }
    
    virtual double c(unsigned long n, unsigned short g) override {
      if (Z_ == 0) {
        solve();
      }
      
      return c_[n-1][g-1];
    }
  };
}

#endif

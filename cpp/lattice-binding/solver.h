//
//  lattice_solver.h
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#ifndef lattice_binding_lattice_solver_h
#define lattice_binding_lattice_solver_h

#include "particle.h"

namespace lattice {
  class Solver {
  protected:
    std::vector<Particle> particles_;
    std::size_t N_;
    unsigned int max_m_, max_V_;
    
    Solver(const std::vector<Particle>& particles,
           std::size_t N) : particles_(particles), N_(N) {
      max_m_ = 0;
      max_V_ = 0;
      for (const Particle& p: particles_) {
        if (p.max_interaction_distance() > max_V_) {
          max_V_ = p.max_interaction_distance();
        }
        if (p.size() > max_m_) {
          max_m_ = p.size();
        }
      }
    }
    
    double K(unsigned long n, unsigned short g) const {
      return particles_[g-1].K(n-1);
    }
    
    double K(unsigned long n, unsigned short g, unsigned short h) const {
      return particles_[g-1].K(n-1, h-1);
    }
    
    double K(unsigned long n, unsigned short g, unsigned short h1, unsigned short h2) const {
      double product = 1;
      for (unsigned short h = h1+1; h <= m(g)-h2; h++) {
        product *= K(N_-m(g)+h, g, h);
      }
      return product;
    }
    
    double dK(unsigned long n, unsigned short g, unsigned short h) const {
      return particles_[g-1].dK(n-1, h-1);
    }
    
    double c0(unsigned short g) const {
      return particles_[g-1].c0();
    }
    
    unsigned int m(unsigned short g) const {
      return particles_[g-1].size();
    }
    
    unsigned int V(unsigned short g) const {
      return particles_[g-1].max_interaction_distance();
    }
    
    double Unwrap(unsigned int h, unsigned short g) {
      if (h == 0) return 1;
      return particles_[g-1].unwrap();
    }
    
    double w(unsigned int l, unsigned short g1, unsigned short g2) {
      return 1;
    }
    
    virtual ~Solver() {}
    
  public:
    virtual double Z() = 0;
    virtual double c(unsigned long n, unsigned short g) = 0;
    
    const std::vector<Particle>& particles() const { return particles_; }
    std::size_t N() const { return N_; }
  };
}

#endif

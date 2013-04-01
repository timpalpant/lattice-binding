//
//  particle.h
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#ifndef lattice_binding_particle_h
#define lattice_binding_particle_h

#include <cmath>

namespace lattice {
  class Particle {
  private:
    std::string name_;
    unsigned int size_, max_interaction_distance_;
    double chemical_potential_, potential_, beta_, unwrap_;
    
  public:
    Particle(unsigned int size, unsigned int max_interaction_distance,
             double chemical_potential, double potential, double unwrap) {
      size_ = size;
      max_interaction_distance_ = max_interaction_distance;
      chemical_potential_ = chemical_potential;
      potential_ = potential;
      unwrap_ = unwrap;
    }
    
    std::string name() const { return name_; }
    unsigned int size() const { return size_; }
    unsigned int max_interaction_distance() const { return max_interaction_distance_; }
    double chemical_potential() const { return chemical_potential_; }
    double potential(unsigned long n) const { return potential_; }
    double K(unsigned long n) const { return exp(-beta() * potential(n)); }
    double K(unsigned long n, unsigned short h) const { return pow(K(n), 1.0/size()); }
    double dK(unsigned long n, unsigned short h) const {
      if (unwrap() == 0) return 1;
      return pow(K(n,h), 1.0/size()-1) / size();
    }
    double c0() const { return exp(beta() * chemical_potential()); }
    double unwrap() const { return unwrap_; }
    double beta() const { return beta_; }
    void set_beta(double beta) { beta_ = beta; }
  };
}

#endif

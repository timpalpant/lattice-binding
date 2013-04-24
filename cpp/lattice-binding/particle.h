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
#include <boost/property_tree/ptree.hpp>
#include "config_error.h"

namespace lattice {
  class particle_error : public config_error {
  public:
    explicit particle_error(const std::string& what_arg)
      : config_error(what_arg) { }
    explicit particle_error(const char* what_arg)
      : config_error(what_arg) { }
  };
  
  class Particle {
  private:
    std::string name_;
    unsigned int size_, max_interaction_distance_;
    double chemical_potential_, beta_;
    std::vector<double> potential_;
    
  public:
    Particle(const std::string& name) : name_(name) { }
    void configure(const boost::property_tree::ptree& pt) throw (particle_error);
    
    std::string name() const { return name_; }
    unsigned int size() const { return size_; }
    unsigned int max_interaction_distance() const { return max_interaction_distance_; }
    double chemical_potential() const { return chemical_potential_; }
    double potential(std::size_t n) const { return potential_[n]; }
    double K(std::size_t n) const { return exp(-beta() * potential(n)); }
    double K(std::size_t n, unsigned short h) const { return pow(K(n), 1.0/size()); }
    double dK(std::size_t n, unsigned short h) const {
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

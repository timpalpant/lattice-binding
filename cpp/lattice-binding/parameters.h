//
//  parameters.h
//  lattice-binding
//
//  Created by Timothy Palpant on 3/31/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#ifndef lattice_binding_parameters_h
#define lattice_binding_parameters_h

#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>
#include "particle.h"

namespace lattice {
  enum class SOLVER {DYNAPRO, TRANSFER_MATRIX};
  
  class config_error : public std::runtime_error {
  public:
    explicit config_error(const std::string& what_arg)
      : std::runtime_error(what_arg) { }
    explicit config_error(const char* what_arg)
      : std::runtime_error(what_arg) { }
  };
  
  class Parameters {
  private:
    double temperature_;
    std::size_t lattice_size_;
    std::vector<Particle> particles_;
    SOLVER solver_;
    boost::filesystem::path output_;
    
    Parameters() { }
    
  public:
    static Parameters load(const boost::filesystem::path& p) throw (config_error);
    static Parameters for_argv(int argc, const char* argv[]) throw (config_error);
    
    void update(const boost::property_tree::ptree& pt) throw (config_error);
    
    const std::vector<Particle>& particles() const { return particles_; }
    double temperature() const { return temperature_; }
    std::size_t lattice_size() const { return lattice_size_; }
    SOLVER solver() const { return solver_; }
    boost::filesystem::path output() const { return output_; }
  };
}

#endif

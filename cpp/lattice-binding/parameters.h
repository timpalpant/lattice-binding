//
//  parameters.h
//  lattice-binding
//
//  Created by Timothy Palpant on 3/31/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#ifndef lattice_binding_parameters_h
#define lattice_binding_parameters_h

#include <boost/filesystem.hpp>
#include "config_error.h"
#include "particle.h"

namespace lattice {
  enum class SolverType {kDynaPro, kTransferMatrix};
  enum class BoundaryCondition {kFixed, kPeriodic};
  
  class Parameters {
  private:
    double beta_;
    std::size_t lattice_size_;
    BoundaryCondition bc_;
    std::vector<Particle> particles_;
    SolverType solver_;
    boost::filesystem::path output_;
    
  public:
    static Parameters load(const boost::filesystem::path& p) throw (config_error);
    static Parameters for_argv(int argc, const char* argv[]) throw (config_error);
    
    const std::vector<Particle>& particles() const { return particles_; }
    void add_particle(const Particle& p) { particles_.push_back(p); }
    double beta() const { return beta_; }
    void set_beta(const double beta) { beta_ = beta; }
    std::size_t lattice_size() const { return lattice_size_; }
    void set_lattice_size(const std::size_t N) { lattice_size_ = N; }
    BoundaryCondition bc() const { return bc_; }
    void set_boundary_condition(const BoundaryCondition bc) { bc_ = bc; }
    SolverType solver() const { return solver_; }
    void set_solver(SolverType s) { solver_ = s; }
    boost::filesystem::path output() const { return output_; }
    void set_output(const boost::filesystem::path& p) { output_ = p; }
  };
}

#endif

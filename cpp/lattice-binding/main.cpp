//
//  main.cpp
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "parameters.h"
#include "transfer_matrix_solver.h"
#include "dynapro_solver.h"

lattice::Solver* initSolver(const lattice::Parameters& params) throw (lattice::config_error) {
  switch (params.solver()) {
    case lattice::SOLVER::DYNAPRO:
      return new lattice::DynaProSolver(params.particles(),
                                        params.lattice_size());
    case lattice::SOLVER::TRANSFER_MATRIX:
      return new lattice::TransferMatrixSolver(params.particles(),
                                               params.lattice_size());
    default:
      throw lattice::config_error("Unknown solver");
  }
}

int main(int argc, const char * argv[]) {
  if (argc < 3) {
    std::cerr << "USAGE: lattice-binding [--include ARK] [--cfg KEY=VALUE]" << std::endl;
    return 2;
  }
  
  lattice::Parameters params = lattice::Parameters::for_argv(argc, argv);
  lattice::Solver* solver = initSolver(params);
  
  std::ofstream of(params.output().string());
  of << "# Position";
  for (const lattice::Particle& p : solver->particles()) {
    of << '\t' << p.name();
  }
  of << std::endl;
  for (std::size_t n = 1; n <= solver->N(); n++) {
    of << n;
    for (std::size_t g = 1; g <= solver->particles().size(); g++) {
      of << '\t' << solver->c(n, g);
    }
    of << std::endl;
  }
}


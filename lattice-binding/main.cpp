//
//  main.cpp
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#include <iostream>
#include <vector>
#include <fstream>
#include "ark.h"
#include "transfer_matrix_solver.h"
#include "dynapro_solver.h"

const double BOLTZMANN = 0.0019872041;

unsigned long lattice_size(config::Ark& config) {
  return 8;
}

double temperature(config::Ark& config) {
  return 1.0 / BOLTZMANN;
}

std::vector<lattice::Particle> init_particles(config::Ark& config) {
  std::vector<lattice::Particle> particles;
  lattice::Particle p(3, 0, 0.0, 0.0, 1.0);
  p.set_beta(1.0 / (BOLTZMANN * temperature(config)));
  particles.push_back(p);
  return particles;
}

int main(int argc, const char * argv[]) {
  if (argc < 2) {
    std::cout << "USAGE: lattice-binding [--include ARK] [--cfg KEY=VALUE]" << std::endl;
    return 2;
  }
  
  //config::Ark ark = config::Ark::fromArgv(argv);
  config::Ark ark;
  unsigned long N = lattice_size(ark);
  std::vector<lattice::Particle> particles = init_particles(ark);
  lattice::Solver* solver = new lattice::TransferMatrixSolver(particles, N);
  
  std::string output("distribution.txt");
  std::ofstream of(output);
  for (int n = 1; n <= N; n++) {
    of << n;
    for (int g = 1; g <= particles.size(); g++) {
      of << '\t' << solver->c(n, g);
    }
    of << std::endl;
  }
}


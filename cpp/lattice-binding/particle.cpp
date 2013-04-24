//
//  particle.cpp
//  kmc
//
//  Created by Timothy Palpant on 3/30/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#include "particle.h"

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>

namespace lattice {
  std::vector<double> load_potential(const boost::filesystem::path& p) {
    std::cout << "Loading potential from: " << p << std::endl;
    std::vector<double> potential;
    std::ifstream f(p.c_str());
    std::string line;
    while (std::getline(f, line)) {
      potential.push_back(boost::lexical_cast<double>(line));
    }
    f.close();
    return potential;
  }
    
  void Particle::configure(const boost::property_tree::ptree& pt) throw (particle_error) {
    size_ = pt.get<std::size_t>("size");
    std::cout << "Particle has size = " << size_ << std::endl;
    for (std::size_t i = 0; i < size_; i++) {
      substates_.push_back(state_->substate(i));
    }

    boost::optional<boost::filesystem::path> p;
    p = pt.get_optional<boost::filesystem::path>("potential");
    if (p) {
      potential_ = load_potential(p.get());
    }
  }
  
}

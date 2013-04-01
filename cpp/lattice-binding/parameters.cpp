//
//  parameters.cpp
//  lattice-binding
//
//  Created by Timothy Palpant on 3/31/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#include "parameters.h"

namespace lattice {
  const double BOLTZMANN = 0.0019872041;
  
  Parameters Parameters::load(const boost::filesystem::path& p) throw (config_error) {
    
  }
  
  Parameters Parameters::for_argv(int argc, const char* argv[]) throw (config_error) {
    
  }
  
  void Parameters::update(const boost::property_tree::ptree& pt) throw (config_error) {
    
  }
}
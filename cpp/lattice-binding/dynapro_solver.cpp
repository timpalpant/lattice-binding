//
//  dynapro_solver.cpp
//  lattice-binding
//
//  Created by Timothy Palpant on 3/17/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#include "dynapro_solver.h"

namespace lattice {
  
  double DynaProSolver::Z(long n) {
    if (n < 1) {
      return 0;
    } else if (n < max_m_) {
      return 1;
    } else if (!Z_cached_[n]) {
      double Zn = Z(n-1);
      size_t f = particles_.size();
      for (unsigned short g = 1; g <= f; g++) {
        Zn += K(n-m(g)+1, g) * c0(g) * Z(n-m(g)-max_V_);
      }
      
      for (unsigned short g = 1; g <= f; g++) {
        for (unsigned short gPrime = 1; gPrime <= f; gPrime++) {
          for (unsigned int l = 0; l <= max_V_; l++) {
            Zn += w(l, gPrime, g) * K(n-m(g)+1, g) * c0(g) * (Z(n-m(g)-l) - Z(n-m(g)-l-1));
          }
        }
      }
      
      Z_[n] = Zn;
      Z_cached_[n] = true;
    }
    
    return Z_[n];
  }
  
  double DynaProSolver::dZdK(long n, unsigned short g) {
    if (n < m(g)) {
      return 0;
    } else if (!dZdK_cached_[n][g]) {
      double dZdKng = dZdK(n-1, g);
      size_t f = particles_.size();
      for (unsigned short g = 1; g <= f; g++) {
        dZdKng += K(n-m(g)+1, g) * c0(g) * Z(n-m(g)-max_V_);
      }
      
      for (unsigned short g = 1; g <= f; g++) {
        for (unsigned short gPrime = 1; gPrime <= f; gPrime++) {
          for (unsigned int l = 0; l <= max_V_; l++) {
            dZdKng += w(l, gPrime, g) * K(n-m(g)+1, g) * c0(g) * (Z(n-m(g)-l) - Z(n-m(g)-l-1));
          }
        }
      }
      
      dZdK_[n][g] = dZdKng;
      dZdK_cached_[n][g] = true;
    }
    
    return dZdK_[n][g];
  }
  
}
//
//  transfer_matrix_solver.cpp
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#include <iostream>
#include "transfer_matrix_solver.h"
#include "Eigen/StdVector"

namespace lattice {

  void TransferMatrixSolver::solve() {
    unsigned int R = num_states();
    std::size_t f = particles_.size();
    std::cout << "Transfer matrix for " << f << " particles has " << R << " states" << std::endl;
    Eigen::VectorXd J = Eigen::VectorXd::Constant(R, 1);
    Eigen::VectorXd A = Eigen::VectorXd::Constant(R, 1);
    std::vector<std::vector<Eigen::VectorXd>> dAdK(N_);
    for (unsigned long n = 1; n <= N_; n++) {
      dAdK[n-1] = std::vector<Eigen::VectorXd>(f);
    }
    
    std::cout << "Calculating partition function and its derivatives" << std::endl;
    Eigen::SparseMatrix<double> Qn(R, R);
    Eigen::SparseMatrix<double> dQdKn(R, R);
    for (unsigned long n = 1; n <= N_; n++) {
      Q(Qn, n);
      // std::cout << "Q(" << n << ") has ";
      // std::cout << Qn.nonZeros() << "/" << Qn.size() << " nonzero elements" << std::endl;
      // std::cout << Qn << std::endl;
      for (unsigned short g = 1; g <= f; g++) {
        dQdK(dQdKn, n, g);
        // std::cout << "dQdK(" << n << "," << g << ") has ";
        // std::cout << dQdKn.nonZeros() << "/" << dQdKn.size() << " nonzero elements" << std::endl;
        dAdK[n-1][g-1] = A.transpose() * dQdKn;
        for (unsigned long i = 1; i < n; i++) {
          dAdK[i-1][g-1] = dAdK[i-1][g-1].transpose() * Qn;
        }
      }
      
      A = A.transpose() * Qn;
    }
    Z_ = A.dot(J);
    std::cout << "Partition function Z = " << Z_ << std::endl;
    
    std::cout << "Calculating particle distributions" << std::endl;
    c_ = std::vector<std::vector<double>>(N_);
    for (unsigned long n = 1; n <= N_; n++) {
      c_[n-1] = std::vector<double>(f);
      for (unsigned short g = 1; g <= f; g++) {
        double dZdK = dAdK[n-1][g-1].dot(J);
        c_[n-1][g-1] = dZdK * K(n,g) / Z_;
      }
    }
  }
  
  // Helper method to compute sums over member functions
  template <typename T, typename return_type>
  return_type sum(TransferMatrixSolver* solver,
                  return_type (TransferMatrixSolver::*function)(T) const,
                  T start, T stop) {
    return_type total = 0;
    for (T k = start; k <= stop; k++) {
      total += (solver->*function)(k);
    }
    return total;
  }
  
  void TransferMatrixSolver::Q(Eigen::SparseMatrix<double>& Q, unsigned long n) {
    int i, j;
    size_t f = particles_.size();
    Q.setZero();
    
    // 1. First unit of g-type protein followed by second unit
    for (unsigned short g = 1; g <= f; g++) {
      i = 1;
      //sum(this, &TransferMatrixSolver::m, 1, g-1);
      for (int k = 1; k <= g-1; k++) {
        i += m(k);
      }
      j = i + 1;
      
      if (n == 1) {
        Q.insert(i-1, j-1) = K(n, g, 1) * c0(g);
      } else if (1 < n && n <= N_ - (m(g)-1)) {
        Q.insert(i-1, j-1) = K(n, g, 1);
      }
    }
    
    // 2. Bound unit of g-type protein, not at the ends of the protein,
    // followed by another bound unit (if m > 2)
    for (unsigned short g = 1; g <= f; g++) {
      for (unsigned int h = 1; h <= m(g)-2; h++) {
        i = h + 1;
        for (unsigned short k = 1; k <= g-1; k++) {
          i += m(k);
        }
        j = i + 1;
        
        if (h < n && n <= N_ - (m(g)-h-1)) {
          Q.insert(i-1, j-1) = K(n,g,h+1);
        }
      }
    }
    
    // 3. Bound unit (m(g)-h) of g-type protein followed by a right
    // free DNA end (h protein units are unbound)
    for (unsigned short g = 1; g <= f; g++) {
      for (unsigned int h = 0; h < m(g); h++) {
        i = -h;
        for (unsigned short k = 1; k <= g; k++) {
          i += m(k);
        }
        j = 2;
        for (unsigned short k = 1; k <= f; k++) {
          j += m(k);
        }
        
        if (m(g)-h <= n && n <= N_-h) {
          if (n == 1) {
            Q.insert(i-1, j-1) = K(n, g, m(g)-h) * c0(g) * Unwrap(h, g);
          } else if (n > 1) {
            Q.insert(i-1, j-1) = K(n, g, m(g)-h) * Unwrap(h, g);
          }
        }
      }
    }
    
    // 4. Bound unit (m_g1 - h_1) of g1-type protein followed by unit 1
    // of g2-protein (no gap between proteins; h1 units of g1-protein
    // are unbound
    for (unsigned short g1 = 1; g1 <= f; g1++) {
      for (unsigned short g2 = 1; g2 <= f; g2++) {
        j = 1;
        for (unsigned short k = 1; k <= g2-1; k++) {
          j += m(k);
        }
        
        for (unsigned int h1 = 0; h1 < m(g1); h1++) {
          i = -h1;
          for (unsigned short k = 1; k <= g1; k++) {
            i += m(k);
          }
          
          if (m(g1)-h1 <= n && n <= N_-m(g2)) {
            if (n == 1) {
              Q.insert(i-1, j-1) = K(n, g1, m(g1)-h1) * w(0, g1, g2) * Unwrap(h1, g1) * c0(g1) * c0(g2);
            } else if (n > 1) {
              Q.insert(i-1, j-1) = K(n, g1, m(g1)-h1) * w(0, g1, g2) * Unwrap(h1, g1) * c0(g1);
            }
          }
        }
      }
    }
    
    // 5. Last unit of g1-type protein followed by unit h2 + 1 of g2-protein
    // (no gap between proteins; h2 units of g2-protein are unbound)
    for (unsigned short g1 = 1; g1 <= f; g1++) {
      i = 0;
      for (unsigned short k = 1; k <= g1; k++) {
        i += m(k);
      }
      
      for (unsigned short g2 = 1; g2 <= f; g2++) {
        for (unsigned int h2 = 1; h2 < m(g2); h2++) {
          j = h2 + 1;
          for (unsigned short k = 1; k <= g2-1; k++) {
            j += m(k);
          }
          
          if (m(g1) <= n && n <= N_ - (m(g2)-h2)) {
            Q.insert(i-1, j-1) = K(n, g1, m(g1)) * w(0, g1, g2) * Unwrap(h2, g2) * c0(g2);
          }
        }
      }
    }
    
    // 6. Left free DNA end continues
    i = 1;
    for (unsigned short k = 1; k <= f; k++) {
      i += m(k);
    }
    j = i;
    Q.insert(i-1, j-1) = 1;
    
    // 7. Right free DNA end continues
    i = 2;
    for (unsigned short k = 1; k <= f; k++) {
      i += m(k);
    }
    j = i;
    if (n > 1) {
      Q.insert(i-1, j-1) = 1;
    }
    
    // 8. Left free DNA end followed by bound unit h+1 of g-type protein
    // (h unbound protein units)
    i = 1;
    for (unsigned short k = 1; k <= f; k++) {
      i += m(k);
    }
    for (unsigned short g = 1; g <= f; g++) {
      for (unsigned int h = 0; h < m(g); h++) {
        j = h + 1;
        for (unsigned short k = 1; k <= g-1; k++) {
          j += m(k);
        }
        
        if (n <= N_ - (m(g)-h)) {
          Q.insert(i-1, j-1) = w(0, 0, g) * Unwrap(h, g) * c0(g);
        }
      }
    }
    
    // 9. Bound unit (m(g) - h) of g-type protein followed by a
    // non-interacting gap longer V(g) (h unbound protein units)
    for (unsigned short g = 1; g <= f; g++) {
      for (unsigned int h = 0; h < m(g); h++) {
        i = -h;
        for (unsigned short k = 1; k <= g; k++) {
          i += m(k);
        }
        
        for (unsigned int l = V(g)+1; l <= max_V_+1; l++) {
          j = l + 2;
          for (unsigned short g = 1; g <= f; g++) {
            j += m(g) + V(g);
          }
          
          if (m(g)-h <= n && n < N_) {
            if (n == 1) {
              Q.insert(i-1, j-1) = w(V(g)+l, g, 0) * K(n, g, m(g)-h) * Unwrap(h, g) * c0(g);
            } else if (n > 1) {
              Q.insert(i-1, j-1) = w(V(g)+l, g, 0) * K(n, g, m(g)-h) * Unwrap(h, g);
            }
          }
        }
      }
    }
    
    // 10. Large noninteracting gap continues, units before max(V(g))
    i = 2;
    for (unsigned short g = 1; g <= f; g++) {
      i += m(g) + V(g);
    }
    for (unsigned int k = 1; k <= max_V_; k++) {
      j = i + k + 1;
      Q.insert(i-1, j-1) = 1;
    }
    
    // 11. Large noninteracting gap continues, units after max(V(g))
    i = max_V_ + 3;
    for (unsigned short g = 1; g <= f; g++) {
      i += m(g) + V(g);
    }
    j = i;
    if (1 < n && n < N_) {
      Q.insert(i-1, j-1) = 1;
    }
    
    // 12. Large noninteracting gap followed by bound unit h+1
    // of g-type protein (h unbound units)
    i = max_V_ + 3;
    for (unsigned short g = 1; g <= f; g++) {
      i += m(g) + V(g);
    }
    for (unsigned short g = 1; g <= f; g++) {
      for (unsigned int h = 0; h < m(g); h++) {
        j = h + 1;
        for (unsigned int k = 1; k <= g-1; k++) {
          j += m(k);
        }
        
        if (1 < n && n <= N_ - (m(g)-h)) {
          Q.insert(i-1, j-1) = Unwrap(h, g) * c0(g);
        }
      }
    }
    
    // 13. Unit (m(g1)-h) of g1-type protein followed by l-unit gap
    // followed by g2-protein
    for (unsigned short g1 = 1; g1 <= f; g1++) {
      for (unsigned int h = 0; h < m(g1); h++) {
        i = -h;
        for (unsigned short k = 1; k <= g1; k++) {
          i += m(k);
        }
        
        for (unsigned short g2 = 1; g2 <= f; g2++) {
          for (unsigned int l = 1; l <= V(g2); l++) {
            j = l + 2;
            for (unsigned short k = 1; k <= f; k++) {
              j += m(k);
            }
            for (unsigned short k = 1; k <= g2-1; k++) {
              j += V(k);
            }
            
            if (m(g1) <= n && n <= N_) {
              if (n == 1) {
                Q.insert(i-1, j-1) = w(l, g1, g2) * K(n, g1, m(g1)-h) * Unwrap(h, g1) * c0(g1);
              } else if (n > 1) {
                Q.insert(i-1, j-1) = w(l, g1, g2) * K(n, g1, m(g1)-h) * Unwrap(h, g1);
              }
            }
          }
        }
      }
    }
    
    // 14. g1-l-g2 gap continues (l free units before g2-type protein)
    for (unsigned short g2 = 1; g2 <= f; g2++) {
      for (unsigned int l = 2; l <= V(g2); l++) {
        i = l + 2;
        for (unsigned short g = 1; g <= f; g++) {
          i += m(g);
        }
        for (unsigned short g = 1; g <= g2-1; g++) {
          i += V(g);
        }
        j = i - 1;
        
        if (1 <= n && n <= N_) {
          Q.insert(i-1, j-1) = 1;
        }
      }
    }
    
    // 15. g1-l-g2 gap followed by bound unit h+1 of g2-type protein
    for (unsigned short g2 = 1; g2 <= f; g2++) {
      for (unsigned int h = 0; h < m(g2); h++) {
        j = h + 1;
        for (unsigned short k = 1; k <= g2-1; k++) {
          j += m(k);
        }
        i = 2;
        for (unsigned short g = 1; g <= f; g++) {
          i += m(g);
        }
        for (unsigned short g = 1; g <= g2-1; g++) {
          i += V(g);
        }
        
        if (1 < n && n < N_-(m(g2)-h)) {
          Q.insert(i-1, j-1) = Unwrap(h, g2) * c0(g2);
        }
      }
    }
  }
  
  void TransferMatrixSolver::dQdK(Eigen::SparseMatrix<double>& dQdK,
                                  unsigned long n, unsigned short g) {
    size_t f = particles_.size();
    int i, j;
    dQdK.setZero();
    
    // 1. First unit of g-type protein followed by second unit
    i = 1;
    //sum(this, &TransferMatrixSolver::m, 1, g-1);
    for (int k = 1; k <= g-1; k++) {
      i += m(k);
    }
    j = i + 1;
    
    if (n == 1) {
      dQdK.insert(i-1, j-1) = dK(n, g, 1) * c0(g);
    } else if (1 < n && n <= N_ - (m(g)-1)) {
      dQdK.insert(i-1, j-1) = dK(n, g, 1);
    }
    
    // 2. Bound unit of g-type protein, not at the ends of the protein,
    // followed by another bound unit (if m > 2)
    for (unsigned int h = 1; h <= m(g)-2; h++) {
      i = h + 1;
      for (unsigned short k = 1; k <= g-1; k++) {
        i += m(k);
      }
      j = i + 1;
      
      if (h < n && n <= N_ - (m(g)-h-1)) {
        dQdK.insert(i-1, j-1) = dK(n,g,h+1);
      }
    }
    
    // 3. Bound unit (m(g)-h) of g-type protein followed by a right
    // free DNA end (h protein units are unbound)
    for (unsigned int h = 0; h < m(g); h++) {
      i = -h;
      for (unsigned short k = 1; k <= g; k++) {
        i += m(k);
      }
      j = 2;
      for (unsigned short k = 1; k <= f; k++) {
        j += m(k);
      }
      
      if (m(g)-h <= n && n <= N_-h) {
        if (n == 1) {
          dQdK.insert(i-1, j-1) = dK(n, g, m(g)-h) * c0(g) * Unwrap(h, g);
        } else if (n > 1) {
          dQdK.insert(i-1, j-1) = dK(n, g, m(g)-h) * Unwrap(h, g);
        }
      }
    }
    
    // 4. Bound unit (m_g1 - h_1) of g1-type protein followed by unit 1
    // of g2-protein (no gap between proteins; h1 units of g1-protein
    // are unbound
    unsigned short g1 = g;
    for (unsigned short g2 = 1; g2 <= f; g2++) {
      j = 1;
      for (unsigned short k = 1; k <= g2-1; k++) {
        j += m(k);
      }
      
      for (unsigned int h1 = 0; h1 < m(g1); h1++) {
        i = -h1;
        for (unsigned short k = 1; k <= g1; k++) {
          i += m(k);
        }
        
        if (m(g1)-h1 <= n && n <= N_-m(g2)) {
          if (n == 1) {
            dQdK.insert(i-1, j-1) = dK(n, g1, m(g1)-h1) * w(0, g1, g2) * Unwrap(h1, g1) * c0(g1) * c0(g2);
          } else if (n > 1) {
            dQdK.insert(i-1, j-1) = dK(n, g1, m(g1)-h1) * w(0, g1, g2) * Unwrap(h1, g1) * c0(g1);
          }
        }
      }
    }
    
    // 5. Last unit of g1-type protein followed by unit h2 + 1 of g2-protein
    // (no gap between proteins; h2 units of g2-protein are unbound)
    i = 0;
    for (unsigned short k = 1; k <= g1; k++) {
      i += m(k);
    }
    
    for (unsigned short g2 = 1; g2 <= f; g2++) {
      for (unsigned int h2 = 1; h2 < m(g2); h2++) {
        j = h2 + 1;
        for (unsigned short k = 1; k <= g2-1; k++) {
          j += m(k);
        }
        
        if (m(g1) <= n && n <= N_ - (m(g2)-h2)) {
          dQdK.insert(i-1, j-1) = dK(n, g1, m(g1)) * w(0, g1, g2) * Unwrap(h2, g2) * c0(g2);
        }
      }
    }
    
    // 6. Left free DNA end continues
    
    // 7. Right free DNA end continues
    
    // 8. Left free DNA end followed by bound unit h+1 of g-type protein
    // (h unbound protein units)
    
    // 9. Bound unit (m(g) - h) of g-type protein followed by a
    // non-interacting gap longer V(g) (h unbound protein units)
    for (unsigned int h = 0; h < m(g); h++) {
      i = -h;
      for (unsigned short k = 1; k <= g; k++) {
        i += m(k);
      }
      
      for (unsigned int l = V(g)+1; l <= max_V_+1; l++) {
        j = l + 2;
        for (unsigned short g = 1; g <= f; g++) {
          j += m(g) + V(g);
        }
        
        if (m(g)-h <= n && n < N_) {
          if (n == 1) {
            dQdK.insert(i-1, j-1) = w(V(g)+l, g, 0) * dK(n, g, m(g)-h) * Unwrap(h, g) * c0(g);
          } else if (n > 1) {
            dQdK.insert(i-1, j-1) = w(V(g)+l, g, 0) * dK(n, g, m(g)-h) * Unwrap(h, g);
          }
        }
      }
    }
    
    // 10. Large noninteracting gap continues, units before max(V(g))
    
    // 11. Large noninteracting gap continues, units after max(V(g))
    
    // 12. Large noninteracting gap followed by bound unit h+1
    // of g-type protein (h unbound units)
    
    // 13. Unit (m(g1)-h) of g1-type protein followed by l-unit gap
    // followed by g2-protein
    for (unsigned short g1 = 1; g1 <= f; g1++) {
      for (unsigned int h = 0; h < m(g1); h++) {
        i = -h;
        for (unsigned short k = 1; k <= g1; k++) {
          i += m(k);
        }
        
        for (unsigned short g2 = 1; g2 <= f; g2++) {
          for (unsigned int l = 1; l <= V(g2); l++) {
            j = l + 2;
            for (unsigned short k = 1; k <= f; k++) {
              j += m(k);
            }
            for (unsigned short k = 1; k <= g2-1; k++) {
              j += V(k);
            }
            
            if (m(g1) <= n && n <= N_) {
              if (n == 1) {
                dQdK.insert(i-1, j-1) = w(l, g1, g2) * dK(n, g1, m(g1)-h) * Unwrap(h, g1) * c0(g1);
              } else if (n > 1) {
                dQdK.insert(i-1, j-1) = w(l, g1, g2) * dK(n, g1, m(g1)-h) * Unwrap(h, g1);
              }
            }
          }
        }
      }
    }
    
    // 14. g1-l-g2 gap continues (l free units before g2-type protein)
    
    // 15. g1-l-g2 gap followed by bound unit h+1 of g2-type protein
  }
  
//  void TransferMatrixSolver::Q(Eigen::SparseMatrix<double>& Q, unsigned long n) {
//    unsigned int R = num_states();
//    size_t f = particles_.size();
//    Q.setZero();
//    
//    // A free DNA unit followed by a free unit
//    Q.insert(R-1, R-1) = 1;
//    
//    // A free unit followed by g type protein
//    for (unsigned short g = 1; g <= f; g++) {
//      int j = 1 - m(g);
//      for (unsigned short k = 1; k <= g; k++) {
//        j += m(k);
//      }
//      
//      if (1 <= n && n <= N_ - (m(g)-1)) {
//        Q.insert(R-1, j-1) = 1;
//      }
//    }
//    
//    // g type protein followed by a free unit
//    for (unsigned short g = 1; g <= f; g++) {
//      int i = 0;
//      for (unsigned short k = 1; k <= g; k++) {
//        i += m(k);
//      }
//      int j = R;
//      
//      if (n >= m(g) && n <= N_) {
//        Q.insert(i-1, j-1) = 1;
//      }
//    }
//    
//    // g1 type protein followed by g2 type protein
//    for (unsigned short g1 = 1; g1 <= f; g1++) {
//      int i = 0;
//      for (unsigned short k = 1; k <= g1; k++) {
//        i += m(k);
//      }
//      
//      for (unsigned short g2 = 1; g2 <= f; g2++) {
//        int j = 1 - m(g2);
//        for (unsigned short k = 1; k <= g2; k++) {
//          j += m(k);
//        }
//        
//        if (m(g1) <= n && n <= N_-(m(g2)-1)) {
//          Q.insert(i-1, j-1) = w(0,g1,g2);
//        }
//      }
//    }
//    
//    // 1st unit inside g-type protein binding site
//    for (unsigned short g = 1; g <= f; g++) {
//      int i = 1 - m(g);
//      for (unsigned short k = 1; k <= g; k++) {
//        i += m(k);
//      }
//      int j = i + 1; // ?
//      
//      if (1 <= n && n <= N_-(m(g)-1)) {
//        Q.insert(i-1, j-1) = K(n,g) * c0(g);
//      }
//    }
//    
//    // hth unit inside g-type protein binding site (h > 1)
//    for (unsigned short g = 1; g <= f; g++) {
//      for (unsigned int h = 2; h <= m(g)-1; h++) {
//        int i = h - m(g);
//        for (unsigned short k = 1; k <= g; k++) {
//          i += m(k);
//        }
//        int j = i+1;
//        
//        if (h <= n && n <= N_-(m(g)-h)) {
//          Q.insert(i-1, j-1) = 1;
//        }
//      }
//    }
//  }
//  
//  void TransferMatrixSolver::dQdK(Eigen::SparseMatrix<double>& dQdK,
//                                  unsigned long n, unsigned short g) {
//    dQdK.setZero();
//    
//    // 1st unit inside g-type protein binding site
//    int i = 1 - m(g);
//    for (unsigned short k = 1; k <= g; k++) {
//      i += m(k);
//    }
//    int j = i + 1;
//    
//    if (1 <= n && n <= N_-(m(g)-1)) {
//      dQdK.insert(i-1, j-1) = c0(g);
//    }
//  }
}
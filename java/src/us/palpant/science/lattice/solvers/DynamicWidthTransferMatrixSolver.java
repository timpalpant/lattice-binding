package us.palpant.science.lattice.solvers;

import us.palpant.science.lattice.Particle;

public class DynamicWidthTransferMatrixSolver extends TransferMatrixSolver {

  public DynamicWidthTransferMatrixSolver(Particle[] particles, int N) {
    super(particles, N);
  }

  @Override
  public TransferMatrix Q(int n) {
    TransferMatrix Q = new TransferMatrix(nStates());
    int i, j;
    
    // 1. First unit of g-type protein followed by the second unit
    for (int g = 1; g <= f; g++) {
      i = 1;
      for (int k = 1; k <= g-1; k++) {
        i += m(k);
      }
      j = i+1;
      
      if (n == 1) {
        Q.setEntry(i, j, K(n,g,1)*c0(g));
      } else if (n <= N - (m(g)-1)) {
        Q.setEntry(i, j, K(n,g,1));
      }
    }
    
    // 2. Bound unit of g-type protein, not at the ends of the protein,
    // followed by another bound unit (if m > 2)
    for (int g = 1; g <= f; g++) {
      for (int h = 2; h <= m(g)-2; h++) {
        i = h + 1;
        for (int k = 1; k <= g-1; k++) {
          i += m(k);
        }
        j = i+1;
        
        if (h < n && n <= N - (m(g)-h-1)) {
          Q.setEntry(i, j, K(n,g,h+1));
        }
      }
    }
    
    // 3. Bound unit (m_g - h) of g-type protein followed by a right
    // free DNA end (h base pairs are unbound)
    j = 2;
    for (int k = 1; k <= f; k++) {
      j += m(k);
    }
    for (int g = 1; g <= f; g++) {
      for (int h = 0; h < m(g); h++) {
        i = -h;
        for (int k = 1; k <= g; k++) {
          i += m(k);
        }
        
        if (m(g)-h <= n && n <= N-h) {
          if (n == 1) {
            Q.setEntry(i, j, K(n,g,m(g)-h)*c0(g)*Unwrap(h,g));
          } else {
            Q.setEntry(i, j, K(n,g,m(g)-h)*Unwrap(h,g));
          }
        }
      }
    }
    
    // 4. Bound unit (m_{g1} - h_1) of g_1-type protein followed by
    // unit 1 of g_2-protein (no gap between proteins; h_1 units of
    // g_1-protein are unbound)
    for (int g1 = 1; g1 <= f; g1++) {
      for (int g2 = 1; g2 <= f; g2++) {
        j = 1;
        for (int k = 1; k <= g2-1; k++) {
          j += m(k);
        }
        
        for (int h1 = 0; h1 < m(g1); h1++) {
          i = -h1;
          for (int k = 1; k <= g1; k++) {
            i += m(k);
          }
          
          if (m(g1)-h1 <= n && n <= N-m(g2)) {
            if (n == 1) {
              Q.setEntry(i, j, K(n,g1,m(g1)-h1)*w(0,g1,g2)*Unwrap(h1,g1)*c0(g1)*c0(g2));
            } else {
              Q.setEntry(i, j, K(n,g1,m(g1)-h1)*w(0,g1,g2)*Unwrap(h1,g1)*c0(g1));
            }
          }
        }
      }
    }
    
    // 5. Last unit of g_1-type protein followed by unit h_2+1 of 
    // g_2-protein (no gap between proteins; h_2 units of g_2-protein 
    // are unbound)
    for (int g1 = 1; g1 <= f; g1++) {
      i = 0;
      for (int k = 1; k <= g1; k++) {
        i += m(k);
      }
      for (int g2 = 1; g2 <= f; g2++) {
        for (int h2 = 1; h2 < m(g2); h2++) {
          j = h2 + 1;
          for (int k = 1; k <= g2-1; k++) {
            j += m(k);
          }
          
          if (m(g1) <= n && n <= N-(m(g2)-h2)) {
            Q.setEntry(i, j, K(n,g1,m(g1))*w(0,g1,g2)*Unwrap(h2,g2)*c0(g2));
          }
        }
      }
    }
    
    // 6. Left free DNA end continues
    i = 1;
    for (int g = 1; g <= f; g++) {
      i += m(g);
    }
    j = i;
    Q.setEntry(i, j, 1);
    
    // 7. Right free DNA end continues
    i = 2;
    for (int g = 1; g <= f; g++) {
      i += m(g);
    }
    j = i;
    Q.setEntry(i, j, 1);
    
    // 8. Left free DNA end followed by bound unit h+1 of g-type
    // protein (h unbound protein units)
    i = 1;
    for (int k = 1; k <= f; k++) {
      i += m(k);
    }
    for (int g = 1; g <= f; g++) {
      for (int h = 0; h < m(g); h++) {
        j = h + 1;
        for (int k = 1; k <= g-1; k++) {
          j += m(k);
        }
        
        if (n <= N-(m(g)-h)) {
          Q.setEntry(i, j, w(0,0,g)*Unwrap(h,g)*c0(g));
        }
      }
    }
    
    // 9. Bound unit (m_g - h) of g-type protein followed by a non-
    // interacting gap longer V_g (h unbound protein units)
    j = 2;
    for (int k = 1; k <= f; k++) {
      j += m(k) + V(k);
    }
    for (int g = 1; g <= f; g++) {
      for (int h = 0; h < m(g); h++) {
        i = -h; // - 1;
        for (int k = 1; k <= g; k++) {
          i += m(k);
        }
        
        if (m(g)-h <= n && n < N) {
          if (n == 1) {
            Q.setEntry(i, j, w(V(g)+1,g,0)*K(n,g,m(g)-h)*Unwrap(h,g)*c0(g));
          } else {
            Q.setEntry(i, j, w(V(g)+1,g,0)*K(n,g,m(g)-h)*Unwrap(h,g));
          }
        }
      }
    }
    
    // 10. Large noninteracting gap continues, units before max(V_g)
    i = 2;
    for (int g = 1; g <= f; g++) {
      i += m(g) + V(g);
    }
    for (int k = 1; k <= maxV(); k++) {
      j = i + k;
      Q.setEntry(j, j+1, 1);
    }
    
    // 11. Large noninteracting gap continues, units after max(V_g)
    i = maxV() + 3;
    for (int g = 1; g <= f; g++) {
      i += m(g) + V(g);
    }
    j = i;
    if (1 < n && n < N) {
      Q.setEntry(i, j, 1);
    }
    
    // 12. Large noninteracting gap followed by bound unit h+1
    // of g-type protein (h unbound protein units)
    i = maxV() + 3;
    for (int g = 1; g <= f; g++) {
      i += m(g) + V(g);
    }
    for (int g = 1; g <= f; g++) {
      for (int h = 0; h < m(g); h++) {
        j = h + 1;
        for (int k = 1; k <= g-1; k++) {
          j += m(k);
        }
        
        if (1 < n && n <= N - (m(g)-h)) {
          Q.setEntry(i, j, Unwrap(h,g)*c0(g));
        }
      }
    }
    
    // 13. Unit (m_{g1} - h) of g_1-type protein followed by l-unit
    // gap followed by g_2-protein
    for (int g1 = 1; g1 <= f; g1++) {
      for (int g2 = 1; g2 <= f; g2++) {
        for (int h = 0; h < m(g1); h++) {
          i = -h;
          for (int k = 1; k <= g1; k++) {
            i += m(k);
          }
          for (int l = 1; l <= V(g2); l++) {
            j = l + 2;
            for (int k = 1; k <= f; k++) {
              j += m(k);
            }
            for (int k = 1; k <= g2-1; k++) {
              j += V(k);
            }
            
            if (m(g1) <= n && n <= N) {
              if (n == 1) {
                Q.setEntry(i, j, w(l,g1,g2)*K(n,g1,m(g1)-h)*Unwrap(h,g1)*c0(g1));
              } else {
                Q.setEntry(i, j, w(l,g1,g2)*K(n,g1,m(g1)-h)*Unwrap(h,g1));
              }
            }
          }
        }
      }
    }
    
    // 14. g_1-l-g_2 gap continues (l free units before g_2-type protein)
    for (int g1 = 1; g1 <= f; g1++) {
      for (int g2 = 1; g2 <= f; g2++) {
        for (int l = 2; l <= V(g2); l++) {
          i = l + 2;
          for (int k = 1; k <= f; k++) {
            i += m(k);
          }
          for (int k = 1; k <= g2-1; k++) {
            i += V(k);
          }
          j = i-1;
          
          // Is this condition necessary?
          if (1 <= n && n <= N) {
            Q.setEntry(i, j, 1);
          }
        }
      }
    }
    
    // 15. g_1-l-g_2 gap followed by bound unit h+1 of g_2-type protein
    for (int g1 = 1; g1 <= f; g1++) {
      for (int g2 = 1; g2 <= f; g2++) {
        i = 2;
        for (int k = 1; k <= f; k++) {
          i += m(k);
        }
        for (int k = 1; k <= g2-1; k++) {
          i += V(k);
        }
        for (int h = 0; h < m(g2); h++) {
          j = h + 1;
          for (int k = 1; k <= g2-1; k++) {
            j += m(k);
          }
          if (1 < n && n < N-m(g2)-h) {
            Q.setEntry(i, j, Unwrap(h,g2)*c0(g2));
          }
        }
      }
    }
    
    return Q;
  }
  
  @Override
  public TransferMatrix dQdK(int n, int g) {
    TransferMatrix dQdK = new TransferMatrix(nStates());
    int i, j;
    

    
    return dQdK;
  }
  
  @Override
  public int nStates() {
    int count = maxV() + 3;
    for (int g = 1; g <= f; g++) {
      count += m(g) + V(g);
    }
    
    return count;
  }

}

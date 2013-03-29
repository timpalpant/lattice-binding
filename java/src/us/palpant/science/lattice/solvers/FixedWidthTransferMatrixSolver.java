package us.palpant.science.lattice.solvers;

import us.palpant.science.lattice.Particle;

public class FixedWidthTransferMatrixSolver extends TransferMatrixSolver {

  public FixedWidthTransferMatrixSolver(Particle[] particles, int N) {
    super(particles, N);
  }

  @Override
  public TransferMatrix Q(int n) {
    int R = nStates();
    TransferMatrix Q = new TransferMatrix(R);
    
    // A free DNA unit followed by a free unit
    Q.setEntry(R, R, 1);
    
    // A free unit followed by g type protein
    for (int g = 1; g <= f; g++) {
      int j = 1 - m(g);
      for (int k = 1; k <= g; k++) {
        j += m(k);
      }
      
      if (1 <= n && n <= N - (m(g)-1)) {
        Q.setEntry(R, j, 1);
      }
    }
    
    // g type protein followed by a free unit
    for (int g = 1; g <= f; g++) {
      int i = 0;
      for (int k = 1; k <= g; k++) {
        i += m(k);
      }
      int j = R;
      
      if (n >= m(g) && n <= N) {
        Q.setEntry(i, j, 1);
      }
    }
    
    // g1 type protein followed by g2 type protein
    for (int g1 = 1; g1 <= f; g1++) {
      int i = 0;
      for (int k = 1; k <= g1; k++) {
        i += m(k);
      }
      
      for (int g2 = 1; g2 <= f; g2++) {
        int j = 1 - m(g2);
        for (int k = 1; k <= g2; k++) {
          j += m(k);
        }
        
        if (m(g1) <= n && n <= N-(m(g2)-1)) {
          Q.setEntry(i, j, w(0,g1,g2));
        }
      }
    }
    
    // 1st unit inside g-type protein binding site
    for (int g = 1; g <= f; g++) {
      int i = 1 - m(g);
      for (int k = 1; k <= g; k++) {
        i += m(k);
      }
      int j = i + 1; // ?
      
      if (1 <= n && n <= N-(m(g)-1)) {
        Q.setEntry(i, j, K(n,g)*c0(g));
      }
    }
    
    // hth unit inside g-type protein binding site (h > 1)
    for (int g = 1; g <= f; g++) {
      for (int h = 2; h <= m(g)-1; h++) {
        int i = h - m(g);
        for (int k = 1; k <= g; k++) {
          i += m(k);
        }
        int j = i+1;
        
        if (h <= n && n <= N-(m(g)-h)) {
          Q.setEntry(i, j, 1);
        }
      }
    }
    
    return Q;
  }

  @Override
  public TransferMatrix dQdK(int n, int g) {
    TransferMatrix dQdK = new TransferMatrix(nStates());
    
    // 1st unit inside g-type protein binding site
    int i = 1 - m(g);
    for (int k = 1; k <= g; k++) {
      i += m(k);
    }
    int j = i + 1;
    
    if (1 <= n && n <= N-(m(g)-1)) {
      dQdK.setEntry(i, j, c0(g));
    }
    
    return dQdK;
  }

  @Override
  public int nStates() {
    int count = 1;
    for (int g = 1; g <= f; g++) {
      count += m(g);
    }
    
    return count;
  }
  
}

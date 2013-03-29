package us.palpant.science.lattice.solvers;

import us.palpant.science.lattice.Particle;

/**
 * Solve lattice binding with a recursive reslations /
 * dynamic programming algorithm
 * 
 * The notation in this class is from:
 * 
 * Teif and Rippe (2011). Calculating transcription factor binding
 * maps for chromatin. Briefings in Bioinformatics 13:2, 187-201.
 * 
 * @author palpant
 *
 */
public class DynaPro implements LatticeSolver {
  
  private final Particle[] particles;
  private final int N, f, V;
  private Double[] Z;
  
  public DynaPro(Particle[] particles, int N) {
    this.particles = particles;
    this.N = N;
    f = particles.length;
    V = getMaxInteractionDistance();
    Z = new Double[N];
    Z[0] = 1.0;
  }

  @Override
  public double c(int n, int g) {
    // TODO Auto-generated method stub
    return 0;
  }

  @Override
  public double Z() {
    return Z(N-1);
  }
  
  /**
   * Solve for the partition function
   * See equation (12) of Teif and Rippe (2011).
   */
  private double Z(int n) {    
    if (Z[n] == null) {
      double Zn = Z(n-1);
      for (int g = 1; g <= f; g++) {
        for (int h1 = 0; h1 <= m(g)-1; h1++) {
          for (int h2 = 0; h2 <= m(g)-h1-1; h2++) {
            Zn += c0(g) * Z(n-m(g)+h1+h2-V) * K(n,g,h1,h2);
          }
        }
      }
      for (int l = 0; l <= V; l++) {
        for (int gPrime = 1; gPrime <= f; gPrime++) {
          for (int g = 1; g <= f; g++) {
            for (int h1 = 0; h1 <= m(g); h1++) {
              for (int h2 = 0; h2 <= m(g)-h1-1; h2++) {
                Zn += w(l,gPrime,g) * c0(g) * (Z(n-m(g)+h1+h2-l) - Z(N-m(g)+h1+h2-l-1)) * K(n,g,h1,h2);
              }
            }
          }
        }
      }
      Z[n] = Zn;
    }
    
    return Z[n];
  }
  
  private int getMaxInteractionDistance() {
    int max = 0;
    for (Particle p : particles) {
      if (p.getMaxInteractionDistance() > max) {
        max = p.getMaxInteractionDistance();
      }
    }
    return max;
  }
  
  private double K(int n, int g, int h1, int h2) {
    double product = 1.0;
    for (int h = h1+1; h <= m(g)-h2; h++) {
      product *= K(n-m(g)+h, g, h);
    }
    return product;
  }
  
  private double K(int n, int g, int h) {
    return particles[g-1].getK(n, h);
  }

  private int m(int g) {
    return particles[g-1].getSize();
  }
  
  private double c0(int g) {
    return particles[g-1].getC0();
  }
  
  private double w(int l, int gPrime, int g) {
    return 1;
  }
  
}

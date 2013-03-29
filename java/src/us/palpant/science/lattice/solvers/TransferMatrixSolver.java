package us.palpant.science.lattice.solvers;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.log4j.Logger;

import us.palpant.science.lattice.Particle;

/**
 * Construct the transfer matrices Q_n for particles on a lattice
 * of length N that can adsorb, desorb, and unwrap
 * 
 * The notations in this class are from:
 * Teif, Ettig, and Rippe. (2010). A Lattice Model for Transcription Factor
 * Access to Nucleosomal DNA. Biophysical J. 99:2597-2607.
 * 
 * @author palpant
 *
 */
public abstract class TransferMatrixSolver implements LatticeSolver {

  private static final Logger log = Logger.getLogger(TransferMatrixSolver.class);
  
  protected final Particle[] particles;
  protected final int f, N;
  protected double Z;
  protected double[][] c;
  
  public TransferMatrixSolver(Particle[] particles, int N) {
    f = particles.length;
    log.info("Initializing transfer matrix for "+f+" particle types on a "+N+" bp lattice");
    this.particles = particles;
    log.info("Transfer matrix will have "+nStates()+" states");
    this.N = N;
  }
  
  @Override
  public double c(int n, int g) {
    if (c == null) {
      solve();
    }
    
    return c[n-1][g-1];
  }
  
  @Override
  public double Z() {
    if (Z == 0) {
      solve();
    }
    
    return Z;
  }
  
  /**
   * Solve for the partition function and protein distributions
   * 
   * @return the protein distributions c[n][g]
   */
  public double[][] solve() {
    RealVector J = new ArrayRealVector(nStates(), 1);
    RealVector A = J.copy();
    RealVector[][] dAdK = new RealVector[N][f];
    
    log.info("Calculating partition function and its derivatives");
    for (int n = 1; n <= N; n++) {
      log.debug("n = "+n);
      TransferMatrix Qn = Q(n);
      log.debug(Qn);
      for (int g = 1; g <= f; g++) {
        dAdK[n-1][g-1] = dQdK(n,g).preMultiply(A);
        for (int i = 1; i < n; i++) {
          dAdK[i-1][g-1] = Qn.preMultiply(dAdK[i-1][g-1]);
        }
      }
      A = Qn.preMultiply(A);
      log.debug(A);
    }
    Z = A.dotProduct(J);
    log.info("Partition function Z = "+Z);
    
    log.info("Calculating particle distributions");
    c = new double[N][f];
    for (int n = 1; n <= N; n++) {
      for (int g = 1; g <= f; g++) {
        double dZdK = dAdK[n-1][g-1].dotProduct(J);
        c[n-1][g-1] = dZdK * K(n,g) / Z;
      }
    }
    
    return c;
  }
  
  /**
   * Get the transfer matrix Q_n for a lattice of latticeSize base pairs
   * @param n the transfer matrix for the nth base pair
   * @param latticeSize the total length of the lattice (for BCs)
   * @return the transfer matrix Q_n
   */
  public abstract TransferMatrix Q(int n);
  
  public abstract TransferMatrix dQdK(int n, int g);
  
  public abstract int nStates();
  
  protected int maxV() {
    int max = 0;
    for (Particle p : particles) {
      if (p.getMaxInteractionDistance() > max) {
        max = p.getMaxInteractionDistance();
      }
    }
    return max;
  }
  
  protected double K(int n, int g) {
    return particles[g-1].getK(n-1);
  }
  
  protected double K(int n, int g, int h) {
    return particles[g-1].getK(n-1, h-1);
  }
    
  protected double dK(int n, int g, int h) {
    return particles[g-1].getDK(n-1, h-1);
  }
  
  protected double c0(int g) {
    return particles[g-1].getC0();
  }
  
  protected int m(int g) {
    return particles[g-1].getSize();
  }
  
  protected int V(int g) {
    return particles[g-1].getMaxInteractionDistance();
  }
  
  protected double Unwrap(int h, int g) {
    return particles[g-1].getUnwrap(h);
  }
  
  protected double w(int l, int g1, int g2) {
    return 1;
  }
  
  protected int countNonZero(RealMatrix m) {
    int count = 0;
    for (int r = 0; r < m.getRowDimension(); r++) {
      for (int c = 0; c < m.getColumnDimension(); c++) {
        if (m.getEntry(r, c) != 0) {
          count++;
        }
      }
    }
    return count;
  }

}

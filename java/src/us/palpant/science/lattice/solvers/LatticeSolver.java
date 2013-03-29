package us.palpant.science.lattice.solvers;

/**
 * Solve lattice binding models
 * 
 * @author palpant
 *
 */
public interface LatticeSolver {

  /**
   * Get the concentration of protein g at n
   * 
   * @param n the base pair in the lattice in [1,N]
   * @param g the protein in [1,f]
   * @return the probability of protein g at base n
   */
  double c(int n, int g);
  
  /**
   * @return the partition function
   */
  double Z();
  
}

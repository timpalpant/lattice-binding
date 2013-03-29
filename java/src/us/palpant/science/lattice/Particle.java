package us.palpant.science.lattice;

import java.nio.file.Path;
import java.nio.file.Paths;

import us.palpant.Ark;
import us.palpant.ArkException;

public class Particle {
  
  private final int size;
  private final int maxInteractionDistance;
  private final double chemicalPotential;
  private final double potential;
  private final double unwrap;
  private final Path output;
  private double beta = 1;
  
  public Particle(int size, int maxInteractionDistance,
                  double chemicalPotential, double potential, double unwrap,
                  Path output) {
    this.size = size;
    this.maxInteractionDistance = maxInteractionDistance;
    this.chemicalPotential = chemicalPotential;
    this.potential = potential;
    this.unwrap = unwrap;
    this.output = output;
  }
  
  public static Particle parse(Ark config) throws ArkException {
    int size = Integer.parseInt((String) config.get("size"));
    double chemicalPotential = Double.parseDouble((String)config.get("chemical_potential"));
    double potential = Double.parseDouble((String)config.get("potential"));
    double unwrap = 0;
    if (config.has("unwrap")) {
      unwrap = Double.parseDouble((String)config.get("unwrap"));
    }
    Path output = null;
    if (config.has("output")) {
      output = Paths.get((String)config.get("output"));
    }
    return new Particle(size, 0, chemicalPotential, potential, unwrap, output);
  }

  /**
   * @return the size
   */
  public int getSize() {
    return size;
  }

  /**
   * @return the maximum interaction distance
   */
  public int getMaxInteractionDistance() {
    return maxInteractionDistance;
  }

  /**
   * @return the chemicalPotential
   */
  public double getChemicalPotential() {
    return chemicalPotential;
  }
  
  /**
   * Get the potential for base n in the lattice
   * @param n
   * @return
   */
  public double getPotential(int n) {
    return potential;
  }

  public double getK(int n) {
    return Math.exp(-beta*getPotential(n));
  }
  
  public double getK(int n, int h) {
    return Math.pow(getK(n), 1.0/getSize());
  }
  
  public double getDK(int n, int h) {
    return Math.pow(getK(n), 1.0/getSize()-1) / getSize();
  }
  
  public double getC0() {
    return Math.exp(beta*getChemicalPotential());
  }
  
  public double getUnwrap(int h) {
    return unwrap;
  }
  
  public final double getBeta() {
    return beta;
  }

  public final void setBeta(double beta) {
    this.beta = beta;
  }

  public Path getOutput() {
    return output;
  }

}

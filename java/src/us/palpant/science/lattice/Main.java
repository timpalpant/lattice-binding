package us.palpant.science.lattice;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;

import org.apache.log4j.Logger;

import us.palpant.Ark;
import us.palpant.ArkException;
import us.palpant.science.lattice.solvers.DynaPro;
import us.palpant.science.lattice.solvers.FixedWidthTransferMatrixSolver;
import us.palpant.science.lattice.solvers.LatticeSolver;

/**
 * Solve for equilibrium binding distributions in a lattice model
 * of factors that can adsorb, desorb, and unwrap.
 * 
 * @author palpant
 *
 */
public class Main {

  private static final Logger log = Logger.getLogger(Main.class);
  
  /**
   * Boltzmann's constant, in kcal/mol
   */
  private static final double BOLTZMANN = 0.0019872041;
  
  private final Ark config;
  
  public Main(Ark config) { 
    this.config = config;
  }
  
  public LatticeSolver getSolver(Particle[] particles, int N) throws ArkException {
    if (!config.has("solver")) {
      throw new ArkException("solver = {transfer_matrix|dynapro} is not specified!");
    }
    String solverName = (String) config.get("solver");
    switch(solverName) {
    case "transfer_matrix":
      log.info("Initializing transfer matrix solver");
      return new FixedWidthTransferMatrixSolver(particles, N);
    case "dynapro":
      log.info("Initializing dynamic programming solver");
      return new DynaPro(particles, N);
    default:
      throw new ArkException("Unknown solver: "+solverName);
    }
  }
  
  public int getLatticeSize() throws ArkException {
    return Integer.parseInt((String)config.get("lattice.length"));
  }
  
  public double getTemperature() throws ArkException {
    double temperature = 1.0 / BOLTZMANN;
    if (config.has("temperature")) {
      temperature = Double.parseDouble((String)config.get("temperature"));
    } else {
      log.warn("No temperature detected in config. Using beta = 1");
    }
    return temperature;
  }
  
  public Particle[] getParticles() throws ArkException {
    log.info("Initializing particles");
    Ark pConfigs = (Ark) config.get("particles");
    List<Particle> particles = new ArrayList<>();
    double beta = 1.0 / (BOLTZMANN * getTemperature());
    log.debug("beta = "+beta);
    for (Entry<String,Object> entry : pConfigs) {
      Ark pConfig = (Ark) entry.getValue();
      Particle p = Particle.parse(pConfig);
      log.info("Initialized particle "+entry.getKey());
      p.setBeta(beta);
      particles.add(p);
    }
    return particles.toArray(new Particle[particles.size()]);
  }
  
  public void run() throws IOException {
    Particle[] particles = getParticles();
    int N = getLatticeSize();
    LatticeSolver solver = getSolver(particles, N);
    
    log.info("Calculating particle distributions and writing to disk");
    for (int g = 1; g <= particles.length; g++) {
      if (particles[g-1].getOutput() != null) {
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(particles[g-1].getOutput(), 
                                                                          Charset.defaultCharset()))) {
          for (int n = 1; n <= N; n++) {
            writer.println(n+"\t"+solver.c(n, g));
          }
        }
      }
    }
  }

  public static void main(String[] args) throws IOException {
    if (args.length < 2) {
      System.err.println("USAGE: us.palpant.science.lattice.Main [--include config.ark] [--cfg PARAM=VALUE]");
      System.exit(2);
    }
    
    log.info("Loading configuration");
    Ark config = Ark.fromArgv(args);
    log.debug("Initializing application");
    Main app = new Main(config);
    app.run();
  }

}

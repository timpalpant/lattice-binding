package us.palpant.science.lattice.solvers;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.log4j.Logger;

public class TransferMatrix extends BlockRealMatrix {

  private static final long serialVersionUID = 1L;
  private static final Logger log = Logger.getLogger(TransferMatrix.class);

  public TransferMatrix(int nStates) {
    super(nStates, nStates);
  }
  
  @Override
  public void setEntry(int i, int j, double value) {
    if (getEntry(i-1,j-1) != 0) {
      log.warn("Setting entry ("+i+","+j+") that has already been set");
    }
    
    if (log.isDebugEnabled()) {
      log.debug("Setting ("+i+","+j+") = "+value);
    }
    super.setEntry(i-1, j-1, value);
  }
  
}

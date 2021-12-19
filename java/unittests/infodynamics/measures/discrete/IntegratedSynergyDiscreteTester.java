package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;
import infodynamics.utils.MatrixUtils;
import junit.framework.TestCase;
import java.util.Arrays;

/**
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class IntegratedSynergyDiscreteTester extends TestCase {

  protected RandomGenerator rg = new RandomGenerator();

  public void testDoubleXOR() throws Exception {

    int D = 4;
    int T = 20000;

    int[][] data_even = rg.generateRandomInts(T, D, 2);
    int[][] data_odd  = rg.generateRandomInts(T, D, 2);
    int[] even_sums = MatrixUtils.sumRows(data_even);
    int[] odd_sums  = MatrixUtils.sumRows(data_odd);
    for (int t=0; t < T; t++) {
      if ((even_sums[t] % 2) != 0) {
        int idx = rg.generateRandomInts(1,D)[0];
        data_even[t][idx] = 1 - data_even[t][idx];
      }
      if ((odd_sums[t] % 2) != 1) {
        int idx = rg.generateRandomInts(1,D)[0];
        data_odd[t][idx] = 1 - data_odd[t][idx];
      }
    }

    IntegratedSynergyCalculatorDiscrete synCalc = new IntegratedSynergyCalculatorDiscrete(2, D);
    synCalc.initialise();
    synCalc.startAddObservations();
    synCalc.addObservations(data_even);
    synCalc.addObservations(data_odd);
    synCalc.finaliseAddObservations();
    assertEquals(1.0, synCalc.computeAverageLocalOfObservations(), 0.01);

  }

}


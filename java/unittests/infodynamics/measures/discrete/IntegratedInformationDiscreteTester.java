package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;
import infodynamics.utils.MatrixUtils;
import junit.framework.TestCase;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class IntegratedInformationDiscreteTester extends TestCase {

  protected RandomGenerator rg = new RandomGenerator();

  public void testAddObservations() throws Exception {

    int[] t0 = { 0, 0, 0, 0};
    int[] t1 = { 1, 1, 1, 0};
    int[] t2 = { 1, 1, 1, 0};
    int[] t3 = { 0, 0, 1, 0};
    int[] t4 = { 1, 1, 1, 0};
    int[] t5 = { 0, 0, 1, 0};
    int[][] states0 = {t0, t1, t2, t3, t4, t5};
    IntegratedInformationCalculatorDiscrete eicd;
    int base = 2, dims = 4;
    boolean exceptionCaught = false;
    IntegratedInformationCalculatorDiscrete iicd;
    iicd = new IntegratedInformationCalculatorDiscrete(base, dims);
    try {
      iicd.initialise();
      iicd.setObservations(states0);
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);
    int[][] states1 = iicd.getData();
    assertTrue(Arrays.equals(states0, states1));

    // Test that adding multiple sets of observations is handled correctly
    IntegratedSynergyCalculatorDiscrete iscd;
    iscd = new IntegratedSynergyCalculatorDiscrete(2, 2);
    int[][] data1 = new int[][] {{1,1}, {0,0}, {0,0}};
    int[][] data2 = new int[][] {{0,1}, {1,1}};
    int[][] data3 = new int[][] {{1,0}, {1,1}};
    //
    iscd.initialise();
    iscd.startAddObservations();
    iscd.addObservations(data1);
    iscd.finaliseAddObservations();
    assertEquals(0.0, iscd.computeAverageLocalOfObservations(), 0.0001);
    //
    iscd.initialise();
    iscd.startAddObservations();
    iscd.addObservations(data1);
    iscd.addObservations(data2);
    iscd.addObservations(data3);
    iscd.finaliseAddObservations();
    assertEquals(1.0, iscd.computeAverageLocalOfObservations(), 0.0001);

  }

  public void testGenerativeModel() {
    // Example adapted from Wikipedia article on IIT.
    // https://en.wikipedia.org/wiki/Integrated_information_theory

    // Generate model.
    int duration = 1000000;

    int[] var0 = rg.generateRandomInts(duration, 2);
    var0[0] = 0;
    int[] var1 = new int[duration];
    var1[0] = 0;
    System.arraycopy(var0, 0, var1, 1, duration - 1);
    int[][] states = {var0, var1};
    int[][] data = MatrixUtils.transpose(states);

    int base = 2, dims = 2;
    boolean exceptionCaught = false;
    IntegratedInformationCalculatorDiscrete iicd;
    iicd = new IntegratedInformationCalculatorDiscrete(base, dims);
    try {
      iicd.initialise();
      iicd.setObservations(data);
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);

    // The integrated information for this model should be very close
    // to 1, but there is a margin of error that needs to be considered.
    double tol = 0.001;
    double res = Double.NaN, res2 = Double.NaN;
    try {
      res = iicd.computeAverageLocalOfObservations();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);
    assertEquals(res, 1, tol);

    // Now test that adding a dummy variable does not change the result
    int[] var2 = new int[duration];
    for (int i = 0; i < duration; i++) {
      var2[i] = 0;
    }
    int[][] states2 = {var0, var1, var2};
    data = MatrixUtils.transpose(states2);
    try {
      iicd = new IntegratedInformationCalculatorDiscrete(base, dims + 1);
      iicd.initialise();
      iicd.setObservations(data);
      res2 = iicd.computeAverageLocalOfObservations();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);
    assertEquals(res, res2, tol); 
  }

  public void testGetterMethods() {

    int base = 2;
    int dims = 2;
    int duration = 1000000;
    boolean exceptionCaught = false;

    IntegratedInformationCalculatorDiscrete iicd =
      new IntegratedInformationCalculatorDiscrete(base, dims);

    try {
      iicd.getSystemInformation();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertTrue(exceptionCaught);
    exceptionCaught = false;

    try {
      iicd.getMinimumInformationPartition();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertTrue(exceptionCaught);
    exceptionCaught = false;

    try {
      iicd.getMinimumInformationPartitionSize();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertTrue(exceptionCaught);
    exceptionCaught = false;

    int[][] data = rg.generateRandomInts(duration, dims, base);
    double res = Double.NaN, tol = 0.001;

    try {
      iicd.initialise();
      iicd.setObservations(data);
      res = iicd.computeAverageLocalOfObservations();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);
    assertEquals(res, 0, tol);

    try {
      iicd.getSystemInformation();
      iicd.getMinimumInformationPartition();
      iicd.getMinimumInformationPartitionSize();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);

  }

  // Utility function to compute integrated measures on data lists, for testing
  public double computeMeasure(IntegratedMeasureCalculatorDiscrete calc, List<int[][]> data)
      throws Exception {
    calc.initialise();
    calc.startAddObservations();
    for (int[][] d : data) { calc.addObservations(d); }
    calc.finaliseAddObservations();
    return calc.computeAverageLocalOfObservations();
  }

  public void testSingleAtoms() throws Exception {
    double tol = 0.00001;

    // Define data lists with a few single-atoms systems
    List<int[][]> rtr = new ArrayList<>(List.of(new int[][]{{0,0},{0,0}},
                                                new int[][]{{1,1},{1,1}}));

    List<int[][]> xtx = new ArrayList<>(List.of(new int[][]{{0,0},{0,0}},
                                                new int[][]{{0,0},{0,1}},
                                                new int[][]{{0,1},{0,0}},
                                                new int[][]{{0,1},{0,1}},
                                                new int[][]{{1,0},{1,0}},
                                                new int[][]{{1,0},{1,1}},
                                                new int[][]{{1,1},{1,0}},
                                                new int[][]{{1,1},{1,1}}));

    List<int[][]> xty = new ArrayList<>(List.of(new int[][]{{0,0},{0,0}},
                                                new int[][]{{0,0},{1,0}},
                                                new int[][]{{0,1},{0,0}},
                                                new int[][]{{0,1},{1,0}},
                                                new int[][]{{1,0},{0,1}},
                                                new int[][]{{1,0},{1,1}},
                                                new int[][]{{1,1},{0,1}},
                                                new int[][]{{1,1},{1,1}}));

    List<int[][]> sts = new ArrayList<>(List.of(new int[][]{{0,0},{0,0}},
                                                new int[][]{{0,0},{1,1}},
                                                new int[][]{{1,1},{0,0}},
                                                new int[][]{{1,1},{1,1}},
                                                new int[][]{{0,1},{0,1}},
                                                new int[][]{{0,1},{1,0}},
                                                new int[][]{{1,0},{0,1}},
                                                new int[][]{{1,0},{1,0}}));

    // Run single-atom systems with integrated synergy
    IntegratedSynergyCalculatorDiscrete psiCalc = new IntegratedSynergyCalculatorDiscrete(2,2);
    assertEquals(computeMeasure(psiCalc, rtr), 0, tol);
    assertEquals(computeMeasure(psiCalc, xtx), 0, tol);
    assertEquals(computeMeasure(psiCalc, xty), 0, tol);
    assertEquals(computeMeasure(psiCalc, sts), 1, tol);

    // Run single-atom systems with WMS integrated info
    IntegratedInformationCalculatorDiscrete phiCalc = new IntegratedInformationCalculatorDiscrete(2,2);
    assertEquals(computeMeasure(phiCalc, rtr), -1, tol);
    assertEquals(computeMeasure(phiCalc, xtx), 0, tol);
    assertEquals(computeMeasure(phiCalc, xty), 1, tol);
    assertEquals(computeMeasure(phiCalc, sts), 1, tol);

    // Run single-atom systems with integrated interaction
    IntegratedInteractionCalculatorDiscrete phiTCalc = new IntegratedInteractionCalculatorDiscrete(2,2);
    assertEquals(computeMeasure(phiTCalc, rtr), 0, tol);
    assertEquals(computeMeasure(phiTCalc, xtx), 0, tol);
    assertEquals(computeMeasure(phiTCalc, xty), 1, tol);
    assertEquals(computeMeasure(phiTCalc, sts), 1, tol);

    // Run single-atom systems with decoding integrated info
    DecoderIntegrationCalculatorDiscrete phiSCalc = new DecoderIntegrationCalculatorDiscrete(2,2);
    assertEquals(computeMeasure(phiSCalc, rtr), 0, tol);
    assertEquals(computeMeasure(phiSCalc, xtx), 0, tol);
    assertEquals(computeMeasure(phiSCalc, xty), 1, tol);
    assertEquals(computeMeasure(phiSCalc, sts), 1, tol);

    // Run single-atom systems with causal density
    CausalDensityCalculatorDiscrete cdCalc = new CausalDensityCalculatorDiscrete(2,2);
    assertEquals(computeMeasure(cdCalc, rtr), 0, tol);
    assertEquals(computeMeasure(cdCalc, xtx), 0, tol);
    assertEquals(computeMeasure(cdCalc, xty), 0.5, tol);
    assertEquals(computeMeasure(cdCalc, sts), 0, tol);

  }

  public void testPartitionScanMethods() throws Exception {
    // Test that partition scans and MIP getter methods work properly. To do
    // this, we create a 4-variable system where 3 of them are
    // parity-preserving random and one of them is fully random.

    // 0. Generate data. Variable 0 is random, the others are PPR
    int D = 4;
    int T = 20000;

    int[][] data_even = rg.generateRandomInts(T, D, 2);
    int[][] data_odd  = rg.generateRandomInts(T, D, 2);
    int[] even_sums = MatrixUtils.sumRows(MatrixUtils.selectColumns(data_even, 1, D-1));
    int[] odd_sums  = MatrixUtils.sumRows(MatrixUtils.selectColumns(data_odd, 1, D-1));
    for (int t=0; t < T; t++) {
      data_even[t][0] = rg.generateRandomInts(1,2)[0];
      if ((even_sums[t] % 2) != 0) {
        int idx = rg.generateRandomInts(1,D-1)[0] + 1;
        data_even[t][idx] = 1 - data_even[t][idx];
      }
      data_odd[t][0] = rg.generateRandomInts(1,2)[0];
      if ((odd_sums[t] % 2) != 1) {
        int idx = rg.generateRandomInts(1,D-1)[0] + 1;
        data_odd[t][idx] = 1 - data_odd[t][idx];
      }
    }

    IntegratedInformationCalculatorDiscrete iiCalc = new IntegratedInformationCalculatorDiscrete(2, D);
    ArrayList<ArrayList<Integer>> true_mip = new ArrayList<>(List.of(new ArrayList<>(List.of(0)), new ArrayList<>(List.of(1,2,3))));

    // 1. Atomic partition
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "ATOMIC");
    iiCalc.initialise();
    iiCalc.startAddObservations();
    iiCalc.addObservations(data_even);
    iiCalc.addObservations(data_odd);
    iiCalc.finaliseAddObservations();
    assertEquals(1.0, iiCalc.computeAverageLocalOfObservations(), 0.01);

    // 2. Even bipartitions
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "EVEN_BIPARTITIONS");
    iiCalc.initialise();
    iiCalc.startAddObservations();
    iiCalc.addObservations(data_even);
    iiCalc.addObservations(data_odd);
    iiCalc.finaliseAddObservations();
    assertEquals(1.0, iiCalc.computeAverageLocalOfObservations(), 0.01);

    // 3. All bipartitions
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "BIPARTITION");
    iiCalc.initialise();
    iiCalc.startAddObservations();
    iiCalc.addObservations(data_even);
    iiCalc.addObservations(data_odd);
    iiCalc.finaliseAddObservations();
    assertEquals(0.0, iiCalc.computeAverageLocalOfObservations(), 0.01);
    assertTrue(true_mip.equals(iiCalc.getMinimumInformationPartition()));

    // 4. All partitions
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "ALL");
    iiCalc.initialise();
    iiCalc.startAddObservations();
    iiCalc.addObservations(data_even);
    iiCalc.addObservations(data_odd);
    iiCalc.finaliseAddObservations();
    assertEquals(0.0, iiCalc.computeAverageLocalOfObservations(), 0.01);
    assertTrue(true_mip.equals(iiCalc.getMinimumInformationPartition()));

  }

}

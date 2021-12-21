/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.measures.continuous.gaussian;

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
public class IntegratedInformationGaussianTester extends TestCase {

  protected RandomGenerator rg = new RandomGenerator();

  public void testGetterMethods() {

    int duration = 1000000;
    boolean exceptionCaught = false;

    IntegratedInformationCalculatorGaussian iicg =
      new IntegratedInformationCalculatorGaussian();

    try {
      iicg.getSystemInformation();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertTrue(exceptionCaught);
    exceptionCaught = false;

    try {
      iicg.getMinimumInformationPartition();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertTrue(exceptionCaught);
    exceptionCaught = false;

    try {
      iicg.getMinimumInformationPartitionSize();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertTrue(exceptionCaught);
    exceptionCaught = false;

    int dims = 3;
    double[][] data = rg.generateNormalData(duration, dims, 0, 1);
    double res = Double.NaN, tol = 0.001;

    try {
      iicg.initialise(dims);
      iicg.setObservations(data);
      res = iicg.computeAverageLocalOfObservations();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);
    assertEquals(res, 0, tol);

    try {
      iicg.getSystemInformation();
      iicg.getMinimumInformationPartition();
      iicg.getMinimumInformationPartitionSize();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);

  }

  // Utility function to simulate AR processes
  public double[][] simulateAR(double[][] A, int T) throws Exception {
    double[][] X = rg.generateNormalData(T, A[0].length, 0, 1);
    for (int t=1; t < T; t++) {
      X[t] = MatrixUtils.add(X[t], MatrixUtils.matrixProduct(A, X[t-1]));
    }
    return X;
  }

  public void testAR() throws Exception {
    double tol = 0.001;
    int T = 20000;

    // Define data lists with a few single-atoms systems
    double[][] indep = simulateAR(new double[][]{{0,0},{0,0}}, T);
    double[][] xtx   = simulateAR(new double[][]{{0.8,0},{0,0}}, T);
    double[][] xty   = simulateAR(new double[][]{{0,0.8},{0,0}}, T);
    double[][] all   = simulateAR(new double[][]{{0.4,0.4},{0.4,0.4}}, T);

    // Run AR systems with integrated synergy
    IntegratedSynergyCalculatorGaussian psiCalc = new IntegratedSynergyCalculatorGaussian();
    assertEquals(psiCalc.compute(indep), 0, tol);
    assertEquals(psiCalc.compute(xtx),   0, tol);
    assertEquals(psiCalc.compute(xty),   0, tol);
    assertTrue(psiCalc.compute(all) > tol);

    // Run AR systems with WMS integrated info
    IntegratedInformationCalculatorGaussian phiCalc = new IntegratedInformationCalculatorGaussian();
    assertEquals(phiCalc.compute(indep), 0, tol);
    assertEquals(phiCalc.compute(xtx),   0, tol);
    assertTrue(phiCalc.compute(xty) > tol);
    assertTrue(phiCalc.compute(all) > tol);

    // Run AR systems with integrated interaction
    IntegratedInteractionCalculatorGaussian phiTCalc = new IntegratedInteractionCalculatorGaussian();
    assertEquals(phiTCalc.compute(indep), 0, tol);
    assertEquals(phiTCalc.compute(xtx),   0, tol);
    assertTrue(phiTCalc.compute(xty) > tol);
    assertTrue(phiTCalc.compute(all) > tol);

    // Run AR systems with decoding integrated info
    DecoderIntegrationCalculatorGaussian phiSCalc = new DecoderIntegrationCalculatorGaussian();
    assertEquals(phiSCalc.compute(indep), 0, tol);
    assertEquals(phiSCalc.compute(xtx),   0, tol);
    assertTrue(phiSCalc.compute(xty) > tol);
    assertTrue(phiSCalc.compute(all) > tol);

    // Run AR systems with causal density
    CausalDensityCalculatorGaussian cdCalc = new CausalDensityCalculatorGaussian();
    assertEquals(cdCalc.compute(indep), 0, tol);
    assertEquals(cdCalc.compute(xtx),   0, tol);
    assertTrue(cdCalc.compute(xty) > tol);
    assertTrue(cdCalc.compute(all) > tol);

  }

  public void testPartitionScanMethods() throws Exception {
    // Test that partition scans and MIP getter methods work properly. To do
    // this, we create a 4-variable system where 3 of them are
    // strongly coupled and one of them is fully random.

    // 0. Generate data. Variable 0 is random, the others are coupled
    int D = 4;
    int T = 20000;

    double[][] data = simulateAR(new double[][]{{0,0,0,0},{0,0.3,0.3,0.3},{0,0.3,0.3,0.3},{0,0.3,0.3,0.3}}, T);

    DecoderIntegrationCalculatorGaussian iiCalc = new DecoderIntegrationCalculatorGaussian();
    ArrayList<ArrayList<Integer>> true_mip = new ArrayList<>(List.of(new ArrayList<>(List.of(0)), new ArrayList<>(List.of(1,2,3))));
    double tol = 0.001;

    // 1. Atomic partition
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "ATOMIC");
    assertTrue(iiCalc.compute(data) > tol);

    // 2. Even bipartitions
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "EVEN_BIPARTITIONS");
    assertTrue(iiCalc.compute(data) > tol);

    // 3. All bipartitions
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "BIPARTITION");
    assertEquals(0.0, iiCalc.compute(data), tol);
    assertTrue(true_mip.equals(iiCalc.getMinimumInformationPartition()));

    // 4. All partitions
    iiCalc.setProperty("PARTITION_SCAN_METHOD", "ALL");
    assertEquals(0.0, iiCalc.compute(data), tol);
    assertTrue(true_mip.equals(iiCalc.getMinimumInformationPartition()));

  }

  public void testNonnegative() throws Exception {

    int D = 6;
    int T = 50;
    double tol = 0.001;
    double[][] data = rg.generateNormalData(T, D, 0, 1);

    IntegratedSynergyCalculatorGaussian psiCalc = new IntegratedSynergyCalculatorGaussian();
    psiCalc.setProperty("PARTITION_SCAN_METHOD", "ATOMIC");
    assertTrue(psiCalc.compute(data) > tol);

    IntegratedInteractionCalculatorGaussian phiTCalc = new IntegratedInteractionCalculatorGaussian();
    phiTCalc.setProperty("PARTITION_SCAN_METHOD", "ATOMIC");
    assertTrue(phiTCalc.compute(data) > tol);

    DecoderIntegrationCalculatorGaussian phiSCalc = new DecoderIntegrationCalculatorGaussian();
    phiSCalc.setProperty("PARTITION_SCAN_METHOD", "ATOMIC");
    assertTrue(phiSCalc.compute(data) > tol);

    CausalDensityCalculatorGaussian cdCalc = new CausalDensityCalculatorGaussian();
    cdCalc.setProperty("PARTITION_SCAN_METHOD", "ATOMIC");
    assertTrue(cdCalc.compute(data) > tol);
  }

  public void testTau() throws Exception {

    // Generate data with 0 Phi at tau = 1, but positive Phi at tau = 2
    int D = 2;
    int T = 20000;
    double tol = 0.001;
    ArrayList<double[][]> data = new ArrayList<>();
    for (int t = 0; t < T; t++) {
      double[][] X = rg.generateNormalData(3, D, 0, 1);
      X[2][0] += X[0][0] + X[0][1];
      X[2][1] += X[0][0] + X[0][1];
      data.add(X);
    }

    IntegratedSynergyCalculatorGaussian psiCalc = new IntegratedSynergyCalculatorGaussian();
    psiCalc.setProperty("TAU", "1");
    psiCalc.initialise(D);
    psiCalc.startAddObservations();
    for (double[][] X : data) { psiCalc.addObservations(X); }
    psiCalc.finaliseAddObservations();
    assertEquals(0.0, psiCalc.computeAverageLocalOfObservations(), tol);

    psiCalc.setProperty("TAU", "2");
    psiCalc.initialise(D);
    psiCalc.startAddObservations();
    for (double[][] X : data) { psiCalc.addObservations(X); }
    psiCalc.finaliseAddObservations();
    assertTrue(psiCalc.computeAverageLocalOfObservations() > tol);
  }

}

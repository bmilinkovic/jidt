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

package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.gaussian.IntegratedInformationCalculatorGaussian;
import infodynamics.measures.continuous.gaussian.IntegratedInteractionCalculatorGaussian;
import infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian;

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
public class IntegratedInformationKraskovTester extends TestCase {

  protected RandomGenerator rg = new RandomGenerator();

  public void testGetterMethods() {

    int duration = 100;
    boolean exceptionCaught = false;

    IntegratedInformationCalculatorKraskov iicg =
      new IntegratedInformationCalculatorKraskov();

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
    double res = Double.NaN;

    try {
      iicg.initialise(dims);
      iicg.setObservations(data);
      res = iicg.computeAverageLocalOfObservations();
    } catch (Throwable e) {
      exceptionCaught = true;
    }
    assertFalse(exceptionCaught);

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
    double tol = 0.01;
    int T = 20000;

    // Define data lists with a few single-atoms systems
    double[][] indep = simulateAR(new double[][]{{0,0},{0,0}}, T);
    double[][] xtx   = simulateAR(new double[][]{{0.8,0},{0,0}}, T);
    double[][] xty   = simulateAR(new double[][]{{0,0.8},{0,0}}, T);
    double[][] all   = simulateAR(new double[][]{{0.4,0.4},{0.4,0.4}}, T);

    // Run AR systems with WMS integrated info
    IntegratedInformationCalculatorKraskov phiCalc = new IntegratedInformationCalculatorKraskov();
    assertEquals(phiCalc.compute(indep), 0, tol);
    assertEquals(phiCalc.compute(xtx),   0, tol);
    assertTrue(phiCalc.compute(xty) > tol);
    assertTrue(phiCalc.compute(all) > tol);

    // Run AR systems with integrated interaction
    IntegratedInteractionCalculatorKraskov phiTCalc = new IntegratedInteractionCalculatorKraskov();
    assertEquals(phiTCalc.compute(indep), 0, tol);
    assertEquals(phiTCalc.compute(xtx),   0, tol);
    assertTrue(phiTCalc.compute(xty) > tol);
    assertTrue(phiTCalc.compute(all) > tol);

  }

  public void testMatchGaussian() throws Exception {
    // Test that the Kraskov calculators match the Gaussian ones for Gaussian data

    IntegratedInformationCalculatorGaussian phig  = new IntegratedInformationCalculatorGaussian();
    IntegratedInteractionCalculatorGaussian phiTg = new IntegratedInteractionCalculatorGaussian();

    IntegratedInformationCalculatorKraskov phik  = new IntegratedInformationCalculatorKraskov();
    IntegratedInteractionCalculatorKraskov phiTk = new IntegratedInteractionCalculatorKraskov();

    MutualInfoCalculatorMultiVariateGaussian micg = new MutualInfoCalculatorMultiVariateGaussian();

    int T = 20000;
    double[][] data = simulateAR(new double[][]{{0.1, 0.2}, {0.3, 0.4}}, T);

    micg.initialise(2,2);
    micg.setObservations(MatrixUtils.selectRows(data, 0, T-1), MatrixUtils.selectRows(data, 1, T-1));
    double tdmi = micg.computeAverageLocalOfObservations();

    double tol = 0.01;
    assertEquals(phig.compute(data), phik.compute(data), tol);
    assertEquals(phiTg.compute(data), phiTk.compute(data), tol);
    assertEquals(phig.getSystemInformation(), phik.getSystemInformation(), tol);

  }

  public void testTau() throws Exception {

    // Generate data with 0 Phi at tau = 1, but positive Phi at tau = 2
    int D = 2;
    int T = 20000;
    double tol = 0.01;
    ArrayList<double[][]> data = new ArrayList<>();
    for (int t = 0; t < T; t++) {
      double[][] X = rg.generateNormalData(3, D, 0, 1);
      X[2][0] += X[0][0] + X[0][1];
      X[2][1] += X[0][0] + X[0][1];
      data.add(X);
    }

    IntegratedInformationCalculatorKraskov phiCalc = new IntegratedInformationCalculatorKraskov();
    phiCalc.setProperty("TAU", "1");
    phiCalc.initialise(D);
    phiCalc.startAddObservations();
    for (double[][] X : data) { phiCalc.addObservations(X); }
    phiCalc.finaliseAddObservations();
    assertEquals(0.0, phiCalc.computeAverageLocalOfObservations(), tol);

    phiCalc.setProperty("TAU", "2");
    phiCalc.initialise(D);
    phiCalc.startAddObservations();
    for (double[][] X : data) { phiCalc.addObservations(X); }
    phiCalc.finaliseAddObservations();
    assertTrue(phiCalc.computeAverageLocalOfObservations() > tol);
  }

  public void testDeterministicNoise() throws Exception {

    double tol = 0.00001;
    int T = 5000;
    double[][] X = simulateAR(new double[][]{{0.1, 0.2}, {0.3, 0.4}}, T);

    IntegratedInformationCalculatorKraskov phiCalc = new IntegratedInformationCalculatorKraskov();
    phiCalc.setProperty("NOISE_LEVEL_TO_ADD", "false");
    double phi1 = phiCalc.compute(X);
    double phi2 = phiCalc.compute(X);
    phiCalc.setProperty("NOISE_LEVEL_TO_ADD", "1");
    double phi3 = phiCalc.compute(X);

    assertEquals(phi1, phi2, tol);
    assertTrue(phi1 > phi3);

  }
}

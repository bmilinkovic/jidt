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

import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>
 * Computes the stochastic interaction of a continuous multivariate system for
 * a given system partition, assuming that the probability distribution
 * function for these observations is a multivariate Gaussian distribution.
 * </p>
 *
 * <p>
 * Usage is as per the paradigm outlined for {@link EffectiveMeasureCalculator}.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>N. Ay, <a href="https://dx.doi.org/10.3390/e17042432">"Information
 *  geometry on complexity and stochastic interaction"</a>, Entropy 17,
 *  2015.</li>
 *
 *  <li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class StochasticInteractionCalculatorGaussian
   extends EffectiveMeasureCalculatorGaussian {

  public StochasticInteractionCalculatorGaussian() {
    super();
  }


  @Override
  public double computeForSystem() throws Exception {

    EntropyCalculatorMultiVariateGaussian ecg;
    ecg = new EntropyCalculatorMultiVariateGaussian();

    MutualInfoCalculatorMultiVariateGaussian micg;
    micg = new MutualInfoCalculatorMultiVariateGaussian();

    int[] ints1toN = MatrixUtils.range(0, dimensions - 1);
    ecg.initialise(dimensions);
    ecg.setCovariance(MatrixUtils.selectRowsAndColumns(laggedCovariance, ints1toN, ints1toN));
    micg.initialise(dimensions, dimensions);
    micg.setCovariance(laggedCovariance, false);
    systemInformation = ecg.computeAverageLocalOfObservations() - micg.computeAverageLocalOfObservations();

    return systemInformation;

  }

  public double computeForPartition(List<List<Integer>> partition) {
    double si = 0.0;

    try {

      if (!isMIComputed) {
        computeForSystem();
      }

      EntropyCalculatorMultiVariateGaussian ecg;
      ecg = new EntropyCalculatorMultiVariateGaussian();

      MutualInfoCalculatorMultiVariateGaussian micg;
      micg = new MutualInfoCalculatorMultiVariateGaussian();

      double sum = 0;

      for (int i = 0; i < partition.size(); i++) {
        int[] p = MatrixUtils.toIntArray(partition.get(i));
        double lcov[][] = getLaggedCovarianceForPartition(p);
        ecg.initialise(p.length);
        ecg.setCovariance(MatrixUtils.selectRowsAndColumns(covariance, p, p));
        micg.initialise(p.length, p.length);
        micg.setCovariance(lcov, false);
        double pCE = ecg.computeAverageLocalOfObservations() - micg.computeAverageLocalOfObservations();
        sum += pCE;
      }

      // Subtract sum of MI of partitions from the MI of system.
      si = sum - systemInformation;

    } catch (Exception e) {
      e.printStackTrace();
    }

    return si;
  }


  // public double[] computeLocalUsingPreviousObservationsForPartition(double[][] x, List<List<Integer>> partition) throws Exception {

  //   int[] p1 = new int[1];

  //   // Variable name guide:
  //   // condCov = future-conditioned system covariance matrix
  //   // margP1CondCov = future-conditioned p1 covariance matrix. First conditioned, then marginalised
  //   // p1CondCov = p1-future-conditioned p1 covariance matrix. First marginalised, then conditioned
  //   // Same for p2


  //   // Initialise auxiliary arrays
  //   int[] p2 = MatrixUtils.allExcept(p1, dimensions);
  //   int[] p1PlusDim = new int[p1.length];
  //   for (int i = 0; i < p1.length; i++) p1PlusDim[i] = p1[i] + dimensions;
  //   int[] p2PlusDim = new int[p2.length];
  //   for (int i = 0; i < p2.length; i++) p2PlusDim[i] = p2[i] + dimensions;

  //   int[] ints1toN = MatrixUtils.range(0, dimensions - 1);
  //   int[] intsNto2N = MatrixUtils.range(dimensions, 2*dimensions - 1);


  //   // Calculate unconstrained past repertoire covariance matrix
  //   double[][] Gamma = MatrixUtils.selectRowsAndColumns(laggedCovariance, ints1toN, intsNto2N);
  //   double[][] L = MatrixUtils.CholeskyDecomposition(covariance);
  //   double[][] condCov = covarianceGivenObservation(covariance, L, Gamma);
  //   double[][] Lc = MatrixUtils.CholeskyDecomposition(condCov);
  //   double[][] margP1CondCov = MatrixUtils.selectRowsAndColumns(condCov, p1, p1);
  //   double[][] margP2CondCov = MatrixUtils.selectRowsAndColumns(condCov, p2, p2);


  //   // Calculate constrained p1 past repertoire covariance matrix
  //   double[][] p1Sigma   = MatrixUtils.selectRowsAndColumns(covariance, p1, p1);
  //   double[][] Lp1       = MatrixUtils.CholeskyDecomposition(p1Sigma);
  //   double[][] p1Gamma   = MatrixUtils.selectRowsAndColumns(laggedCovariance, p1, p1PlusDim);
  //   double[][] p1CondCov = covarianceGivenObservation(p1Sigma, Lp1, p1Gamma);
  //   double[][] Lp1c      = MatrixUtils.CholeskyDecomposition(p1CondCov);


  //   // Calculate constrained p2 past repertoire covariance matrix
  //   double[][] p2Sigma   = MatrixUtils.selectRowsAndColumns(covariance, p2, p2);
  //   double[][] Lp2       = MatrixUtils.CholeskyDecomposition(p2Sigma);
  //   double[][] p2Gamma   = MatrixUtils.selectRowsAndColumns(laggedCovariance, p2, p2PlusDim);
  //   double[][] p2CondCov = covarianceGivenObservation(p2Sigma, Lp2, p2Gamma);
  //   double[][] Lp2c      = MatrixUtils.CholeskyDecomposition(p2CondCov);


  //   // Calculate terms of KL that don't depend on the means and are the same for all x
  //   double traceTerm = MatrixUtils.trace(MatrixUtils.solveViaCholeskyResult(Lp1c, margP1CondCov))
  //                      + MatrixUtils.trace(MatrixUtils.solveViaCholeskyResult(Lp2c, margP2CondCov));

  //   double detTerm = Math.log(MatrixUtils.determinantViaCholeskyResult(Lp1c))
  //                    + Math.log(MatrixUtils.determinantViaCholeskyResult(Lp2c))
  //                    - Math.log(MatrixUtils.determinantViaCholeskyResult(Lc));

  //   double covTermsSum = traceTerm - dimensions + detTerm;

  //   // Calculate KL divergence between unconstrained and constrained
  //   // Gaussians for each data point (covariances remain, but means
  //   // change for different x)
  //   double[] locals = new double[x.length];
  //   for (int i = 0; i < x.length; i++) {

  //     double squareTerm = 0.0;

  //     double[] condMean = meanGivenObservation(means, Gamma, L, x[i]);
  //     double[] margP1CondMean = MatrixUtils.select(condMean, p1);
  //     double[] margP2CondMean = MatrixUtils.select(condMean, p2);

  //     double[] p1CondMean = meanGivenObservation(MatrixUtils.select(means, p1), p1Gamma, Lp1, MatrixUtils.select(x[i], p1));
  //     double[] p2CondMean = meanGivenObservation(MatrixUtils.select(means, p2), p2Gamma, Lp2, MatrixUtils.select(x[i], p2));

  //     double[] p1muDiff = MatrixUtils.subtract(margP1CondMean, p1CondMean);
  //     double[] p2muDiff = MatrixUtils.subtract(margP2CondMean, p2CondMean);

  //     squareTerm += MatrixUtils.matrixProduct(p1muDiff, MatrixUtils.solveViaCholeskyResult(Lp1c, MatrixUtils.vecTranspose(p1muDiff)))[0];
  //     squareTerm += MatrixUtils.matrixProduct(p2muDiff, MatrixUtils.solveViaCholeskyResult(Lp2c, MatrixUtils.vecTranspose(p2muDiff)))[0];

  //     locals[i] = 0.5*(covTermsSum + squareTerm);

  //   }

  //   return locals;

  // }

  // private static double[][] covarianceGivenObservation(double[][] Sigma, double[][] L, double[][] Gamma) throws Exception {
  //   double[][] prod = MatrixUtils.matrixProduct(MatrixUtils.transpose(Gamma), MatrixUtils.solveViaCholeskyResult(L, Gamma));
  //   return MatrixUtils.subtract(Sigma, prod);
  // }

  // private static double[] meanGivenObservation(double[] mu, double[][] Gamma, double[][] L, double[] x) throws Exception {
  //   // mu_agb = mu_a + Sigma_ab * inv(Sigma_bb) * (x - mu_b)
  //   double[][] diff = MatrixUtils.vecTranspose(MatrixUtils.subtract(x, mu));
  //   double[][] prod = MatrixUtils.matrixProduct(MatrixUtils.transpose(Gamma), MatrixUtils.solveViaCholeskyResult(L, diff));
  //   return MatrixUtils.add(mu, MatrixUtils.vecTranspose(prod));

  // }

}

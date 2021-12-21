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

import infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian;

import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>
 * Implements "whole-minus-sum" effective information beyond a given system
 * partition, using the system's stationary distribution following Barrett and
 * Seth (2011).
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>A. Barrett, <a href="https://dx.doi.org/10.1371/journal.pcbi.1001052">
 *  "Practical measures of integrated information for time-series data"</a>,
 *  , PLoS Comput Biol 7, 2011.</li>
 *
 *  <li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class EffectiveInformationCalculatorGaussian
   extends EffectiveMeasureCalculatorGaussian {

 /**
  * Constructor.
  */
  public EffectiveInformationCalculatorGaussian() {
      super();
  }

  @Override
  public double computeForPartition(List<List<Integer>> partition) {
    double ei = 0.0;

    try {

      if (!isMIComputed) {
        computeForSystem();
      }

      double sum = 0.0;
      MutualInfoCalculatorMultiVariateGaussian micg;
      micg = new MutualInfoCalculatorMultiVariateGaussian();

      for (int i = 0; i < partition.size(); i++) {
        int[] p = MatrixUtils.toIntArray(partition.get(i));
        micg.initialise(p.length, p.length);
        micg.setCovariance(getLaggedCovarianceForPartition(p), covFromObservations);
        double pMI = micg.computeAverageLocalOfObservations();
        sum += pMI;
      }

      // Subtract sum of MI of partitions from the MI of system.
      ei = systemInformation - sum;

    } catch (Exception e) {
      e.printStackTrace();
    }

    return ei;
  }

  // public double[] computeLocalUsingPreviousObservationsForPartition(double[][] x, List<List<Integer>> partition) throws Exception {

  //   // Variable name guide:
  //   // condCov = future-conditioned system covariance matrix
  //   // margP1CondCov = future-conditioned p1 covariance matrix. First conditioned, then marginalised
  //   // p1CondCov = p1-future-conditioned p1 covariance matrix. First marginalised, then conditioned

  //   int[] ints1toN = MatrixUtils.range(0, dimensions - 1);
  //   int[] intsNto2N = MatrixUtils.range(dimensions, 2*dimensions - 1);
  //   int nb_parts = partition.size();


  //   // Calculate unconstrained past repertoire covariance matrix
  //   double[][] Gamma = MatrixUtils.selectRowsAndColumns(laggedCovariance, ints1toN, intsNto2N);
  //   double[][] L = MatrixUtils.CholeskyDecomposition(covariance);
  //   double[][] condCov = covarianceGivenObservation(covariance, L, Gamma);
  //   double[][] Lc = MatrixUtils.CholeskyDecomposition(condCov);

  //   // Iterate through partitions to calculate covariance matrices and
  //   // pre-compute some terms for the KL
  //   double[][][] pSigma   = new double[nb_parts][][];
  //   double[][][] Lp       = new double[nb_parts][][];
  //   double[][][] pGamma   = new double[nb_parts][][];
  //   double[][][] pCondCov = new double[nb_parts][][];
  //   double[][][] Lpc      = new double[nb_parts][][];
  //   double[][][] margPCondCov = new double[nb_parts][][];
  //   double traceTerm = 0.0, detTerm = 0.0;
  //   for (int j = 0; j < nb_parts; j++) {

  //     int[] p = MatrixUtils.toIntArray(partition.get(j));
  //     int[] pPlusDim = new int[p.length];
  //     for (int k = 0; k < p.length; k++) pPlusDim[k] = p[k] + dimensions;

  //     // Calculate marginal part past conditional covariance (condition whole
  //     // system, then marginalise condCov).
  //     margPCondCov[j] = MatrixUtils.selectRowsAndColumns(condCov, p, p);

  //     // Calculate constrained part past repertoire covariance matrix
  //     pSigma[j]   = MatrixUtils.selectRowsAndColumns(covariance, p, p);
  //     Lp[j]       = MatrixUtils.CholeskyDecomposition(pSigma[j]);
  //     pGamma[j]   = MatrixUtils.selectRowsAndColumns(laggedCovariance, p, pPlusDim);
  //     pCondCov[j] = covarianceGivenObservation(pSigma[j], Lp[j], pGamma[j]);
  //     Lpc[j]      = MatrixUtils.CholeskyDecomposition(pCondCov[j]);

  //     // Calculate terms of KL that don't depend on the means and are the same for all x
  //     traceTerm += MatrixUtils.trace(MatrixUtils.solveViaCholeskyResult(Lpc[j], margPCondCov[j]));
  //     detTerm += Math.log(MatrixUtils.determinantViaCholeskyResult(Lpc[j]));
  //   }

  //   detTerm -= Math.log(MatrixUtils.determinantViaCholeskyResult(Lc));
  //   double covTermsSum = traceTerm - dimensions + detTerm;

  //   // Calculate KL divergence between unconstrained and constrained
  //   // Gaussians for each data point (covariances remain, but means
  //   // change for different x)
  //   double[] locals = new double[x.length];
  //   for (int i = 0; i < x.length; i++) {

  //     double squareTerm = 0.0;

  //     double[] condMean = meanGivenObservation(means, Gamma, L, x[i]);

  //     for (int j = 0; j < partition.size(); j++) {
  //       int[] p = MatrixUtils.toIntArray(partition.get(j));
  //       double[] margPCondMean = MatrixUtils.select(condMean, p);
  //       double[] pCondMean = meanGivenObservation(MatrixUtils.select(means, p), pGamma[j], Lp[j], MatrixUtils.select(x[i], p));
  //       double[] pmuDiff = MatrixUtils.subtract(margPCondMean, pCondMean);
  //       squareTerm += MatrixUtils.matrixProduct(pmuDiff, MatrixUtils.solveViaCholeskyResult(Lpc[j], MatrixUtils.vecTranspose(pmuDiff)))[0];
  //     }

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

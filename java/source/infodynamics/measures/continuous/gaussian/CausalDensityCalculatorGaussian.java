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

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;


/**
 * <p>
 * Implements causal density for multivariate contiuous time-series data,
 * assuming a Gaussian distribution.  Causal density is defined as the average
 * of the transfer entropies between system parts (all taken with an embedding
 * length of 1). This quantity is computed for multiple partitions of the
 * system (as per {@link
 * infodynamics.measures.continuous.IntegratedMeasureCalculator}) and the
 * minimum is returned.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>A. Seth, <a href="https://doi.org/10.1098/rsta.2011.0079">"Causal
 *  density and integrated information as measures of conscious level"</a>,
 *  Philos. Trans. A 369, 2011.</li>
 *
 *  <li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class CausalDensityCalculatorGaussian 
    extends IntegratedMeasureCalculatorGaussian {

  /**
   * Covariance matrix of the system with itself {@link #tau} timesteps
   * ahead (symmetric, positive-definite matrix of size 2*dimensions).
   */
  // protected double[][] laggedCovariance;

  /**
   * Stationary covariance matrix of the system (symmetric, positive-definite matrix of size dimensions).
   */
  // protected double[][] covariance;

  /**
   * Whether the covariance has been obtained from samples or supplied
   * directly.
   */
  // protected boolean covFromObservations;


  /**
   * Constructor.
   */
  public CausalDensityCalculatorGaussian() {
    super();
    PARTITION_SCAN_METHOD = "ATOMIC";
    baseCalculator = new CausalDensityHelper();

  }

  /**
   * Mock effective measure to clean up the code. This is the actual causal
   * density for a given partition, these classes are split to recycle more
   * code from IntegratedMeasureCalculatorGaussian.
   */
  private class CausalDensityHelper extends EffectiveMeasureCalculatorGaussian {

    /**
     * Constructor.
     */
    public CausalDensityHelper() {
      super();
    }

    /**
     * {@inheritdoc}
     */
    public double computeForPartition(List<List<Integer>> partition) throws Exception {

      if (!isMIComputed) {
        computeForSystem();
      }

      int nb_parts = partition.size();
      double sum = 0.0;

      int[] ints1toN = MatrixUtils.range(0, dimensions - 1);
      int[] intsNto2N = MatrixUtils.range(dimensions, 2*dimensions - 1);

      // Cross-covariance matrix, Gamma = Sigma( past, present )
      double[][] Gamma = MatrixUtils.selectRowsAndColumns(laggedCovariance, ints1toN, intsNto2N);
      double[][] L = MatrixUtils.CholeskyDecomposition(covariance);

      // Conditioned covariance matrix:
      // Sigma( present | past ) = Sigma - Gamma.T*(Sigma^-1)*Gamma
      double[][] condCov = MatrixUtils.subtract(covariance,
          MatrixUtils.matrixProduct(MatrixUtils.transpose(Gamma),
            MatrixUtils.solveViaCholeskyResult(L, Gamma)));

      for (int i = 0; i < nb_parts; i++) {

        // Above we have calculated the conditional covariance Sigma( present |
        // past ) for the whole system. By the marginalisation property of
        // Gaussian distributions, the variance of each variable is the
        // corresponding diagonal element of the system's covariance.
        int[] p1 = MatrixUtils.toIntArray(partition.get(i));
        double logdetCovPresGivenPast = Math.log(MatrixUtils.determinantSymmPosDefMatrix(
              MatrixUtils.selectRowsAndColumns(condCov, p1, p1)));

        for (int j = 0; j < nb_parts; j++) {
          if (i == j) {
            continue;
          }

          int[] p2 = MatrixUtils.toIntArray(partition.get(j));

          // Now we have to calculate the conditional covariance
          // Sigma( pres | past ) of one variable conditioned on all variables
          // except for j.

          int[] allExceptJ = MatrixUtils.allExcept(p2, dimensions);
          double[][] pGamma = MatrixUtils.selectRowsAndColumns(laggedCovariance, allExceptJ, MatrixUtils.add(p1, dimensions));
          double[][] pCov = MatrixUtils.selectRowsAndColumns(covariance, allExceptJ, allExceptJ);
          double[][] pL = MatrixUtils.CholeskyDecomposition(pCov);

          // Conditioned covariance matrix:
          // Sigma( present | past ) = Sigma - Gamma.T*(Sigma^-1)*Gamma
          double[][] pCondCov = MatrixUtils.subtract(MatrixUtils.selectRowsAndColumns(covariance, p1, p1),
              MatrixUtils.matrixProduct(MatrixUtils.transpose(pGamma),
                MatrixUtils.solveViaCholeskyResult(pL, pGamma)));

          double logdetPCovPresGivenPast = Math.log(MatrixUtils.determinantSymmPosDefMatrix(pCondCov));

          sum += 0.5*(logdetPCovPresGivenPast - logdetCovPresGivenPast);

        }
      }

      double cd = sum/((double) (dimensions*(dimensions - 1)));

      return cd;
    }

  }

}

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

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

import org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>
 * Computes the mutual information between the past and future of a continuous
 * multivariate system using a "mismatched decoder" that ignores correlations
 * between system parts, assuming a Gaussian distribution. This measure is used
 * as the basis for the decoder-based measure of integrated information (Phi-star).
 * </p>
 *
 * <p>
 * Calculation of mismatched information involves a 1D constrained maximization
 * problem. Since this is differentiable and provably convex, it is generally
 * not an issue. Optimisation is performed with a conjugate gradient optimiser
 * from the Apache Commons math library.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 * <li>M. Oizumi, <a href="https://doi.org/10.1371/journal.pcbi.1004654">
 * "Measuring Integrated Information from the Decoding Perspective"</a>, PLoS
 * Comput Biol 12(1), 2016.</li>
 *
 *  <li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class MismatchedInformationCalculatorGaussian
    extends EffectiveMeasureCalculatorGaussian {

  /**
   * Helper class to optimize mismatched information. See documentation
   * in {@link #DecodingCalculatorHelper}.
   */
  protected DecodingCalculatorHelper calcHelper;

  /**
   * Value of beta that maximises mismatched information.
   */
  protected double maxBeta;

  /**
   * Maximum mismatched information used in the main calculation.
   */
  protected double mismatchedInformation;


  /**
   * Constructor.
   */
  public MismatchedInformationCalculatorGaussian() {
    super();
  }


  public double computeForPartition(List<List<Integer>> partition) throws Exception {

    if (!isMIComputed) {
      computeForSystem();
    }

    calcHelper = new DecodingCalculatorHelper(laggedCovariance, partition);

    // Optimize Istar to calculate mismatched decoding information
    NonLinearConjugateGradientOptimizer opt
        = new NonLinearConjugateGradientOptimizer(
            NonLinearConjugateGradientOptimizer.Formula.POLAK_RIBIERE,
            new SimpleValueChecker(1e-6, 1e-6),
            1e-4, 1e-4, 1);

    try {
      PointValuePair optimum
          = opt.optimize(new MaxEval(Integer.MAX_VALUE),
              calcHelper.getObjectiveFunction(),
              calcHelper.getObjectiveFunctionGradient(),
              GoalType.MAXIMIZE,
              new InitialGuess(new double[] { -0.6931 })); // == log(0.5)

      maxBeta = Math.exp(optimum.getPoint()[0]);
      mismatchedInformation = optimum.getValue();

    } catch (TooManyEvaluationsException e) {
      if (debug) System.out.println("Optimization exceeded maximum number " +
                                    "of evaluations.");
      maxBeta = Double.NaN;
      mismatchedInformation = Double.NaN;
    }

    return (systemInformation - mismatchedInformation);
  }


  /**
   * Compute value of IStar for the most recent data added for the given
   * value of beta.
   */
  public double computeValueForBeta(double beta) throws Exception {
    MultivariateFunction f = calcHelper.getObjectiveFunction().getObjectiveFunction();
    return f.value(new double[] {Math.log(beta)});
  }

  /**
   * Compute gradient of IStar for the most recent data added for the given
   * value of beta.
   */
  public double computeGradientForBeta(double beta) throws Exception {
    MultivariateVectorFunction f = calcHelper.getObjectiveFunctionGradient().getObjectiveFunctionGradient();
    return (f.value(new double[] {Math.log(beta)})[0]/beta);
  }


  /**
   * Helper class
   */
  private class DecodingCalculatorHelper {
    double[][] a, b, c, dQ, invSigma;
    double detSigma, term3;
    int dimensions;

    /**
     * Constructor. Pre-calculate as many things as possible to make
     * the calculation of I and dIdbeta faster.
     *
     * @param lcov full time-lagged covariance matrix of the system
     */
    public DecodingCalculatorHelper(double[][] lcov, List<List<Integer>> partition) throws Exception {
      // Pre-calculations (without beta)
      dimensions = lcov.length/2;
      int nb_parts = partition.size();
      int[] ints1toN = MatrixUtils.range(0, dimensions - 1);
      int[] intsNto2N = MatrixUtils.range(dimensions, 2*dimensions - 1);

      double[][] Gamma = MatrixUtils.selectRowsAndColumns(laggedCovariance, ints1toN, intsNto2N);
      double[][] L = MatrixUtils.CholeskyDecomposition(covariance);
      detSigma = MatrixUtils.determinantViaCholeskyResult(L);

      // NOTE: this is somewhat of a hack. Because the Chol inversion isn't very numerically
      // stable, invSigma might not be symmetric (difference is usually ~1e-17), and to fix
      // it we just symmetrize the matrix adding its own transpose. In the future we may
      // consider other options. For some suggestions see:
      // http://scicomp.stackexchange.com/questions/3188/dealing-with-the-inverse-of-a-positive-definite-symmetric-covariance-matrix
      double[][] invSigma_tmp = MatrixUtils.invertSymmPosDefMatrix(covariance);
      invSigma = MatrixUtils.matrixScalarProduct(0.5, MatrixUtils.add(invSigma_tmp, MatrixUtils.transpose(invSigma_tmp)));


      if (nb_parts == dimensions) {
        // We have an atomic partition, so we can use diagonal matrix operations
        
        double[][] PIx = MatrixUtils.diag(covariance);
        double[][] PIxx = MatrixUtils.diag(Gamma);
        double[][] invPIx = MatrixUtils.invertDiagonalMatrix(PIx);

        double[][] invPIxgx = new double[dimensions][dimensions];
        MatrixUtils.fill(invPIxgx, 0.0);
        for (int i = 0; i < dimensions; i++) {
          invPIxgx[i][i] = 1.0/(PIx[i][i] - PIxx[i][i]*invPIx[i][i]*PIxx[i][i]);
        }

        a = MatrixUtils.matrixDiagonalProduct(PIxx, invPIx);
        b = MatrixUtils.matrixDiagonalProduct(invPIxgx, a);
        c = MatrixUtils.matrixDiagonalProduct(a, b); // Technically a.T, but since it's diagonal it doesn't matter

      } else {
        // We have a general (non-atomic) partition

        // First calculate the block diagonal covariances of the mismatched
        // system. First marginalise parts, then condition separately, and
        // then put together again.
        double[][] PIx   = new double[dimensions][dimensions];
        double[][] PIxx  = new double[dimensions][dimensions];
        double[][] PIxgx = new double[dimensions][dimensions];

        for (int i = 0; i < nb_parts; i++) {
          int[] p = MatrixUtils.toIntArray(partition.get(i));
          int part_size = p.length;

          double[][] p_PIx = MatrixUtils.selectRowsAndColumns(covariance, p, p);
          double[][] p_LPIx = MatrixUtils.CholeskyDecomposition(p_PIx);
          double[][] p_PIxx = MatrixUtils.selectRowsAndColumns(Gamma, p, p);

          double[][] p_PIxgx = MatrixUtils.subtract(p_PIx,
              MatrixUtils.matrixProduct(MatrixUtils.transpose(p_PIxx),
                MatrixUtils.solveViaCholeskyResult(p_LPIx, p_PIxx)));

          // Put together the part covariances to form the block diagonal,
          // mismatched system covariances
          for (int j = 0; j < part_size; j++) {
            for (int k = 0; k < part_size; k++) {
              PIx[p[j]][p[k]] = p_PIx[j][k];
              PIxx[p[j]][p[k]] = p_PIxx[j][k];
              PIxgx[p[j]][p[k]] = p_PIxgx[j][k];
            }
          }
        }


        // Now do the same computations as for the atomic case, but with
        // dense matrix operations.
        double[][] invPIx = MatrixUtils.invertSymmPosDefMatrix(PIx);
        double[][] LPIxgx = MatrixUtils.CholeskyDecomposition(PIxgx);

        a = MatrixUtils.matrixProduct(PIxx, invPIx);
        b = MatrixUtils.solveViaCholeskyResult(LPIxgx, a);
        c = MatrixUtils.matrixProduct(MatrixUtils.transpose(a), b);

      }

      term3 = MatrixUtils.trace(MatrixUtils.matrixProduct(b, Gamma));
      dQ = c;

    }

    /**
     * Return function handle to the function that calculates Istar for a
     * particular value of beta.
     */
    public ObjectiveFunction getObjectiveFunction() {
      return new ObjectiveFunction(new MultivariateFunction() {
        public double value(double[] point) {
          // Value calculation with beta
          double beta = Math.exp(point[0]);

          try {
            // Since Q is the sum of positive-definite matrices, it is also positive-definite
            double[][] Q = MatrixUtils.add(invSigma, MatrixUtils.matrixScalarProduct(beta, c));
            double[][] invQ = MatrixUtils.invertSymmPosDefMatrix(Q);
            double detQ = MatrixUtils.determinantSymmPosDefMatrix(Q);

            double[][] d = MatrixUtils.matrixProduct(b, MatrixUtils.matrixProduct(invQ, MatrixUtils.transpose(b)));
            double[][] R = MatrixUtils.add( MatrixUtils.matrixScalarProduct(-1*beta, c), MatrixUtils.matrixScalarProduct(-1*beta*beta, d));

            double term1 = 0.5*(Math.log(detQ) + Math.log(detSigma));
            double term2 = 0.5*MatrixUtils.trace(MatrixUtils.matrixProduct(covariance, R));

            return (term1 + term2 + beta * term3);
          } catch (Throwable e) {
            e.printStackTrace();
            return Double.NaN;
          }
        }
      });
    }

    /**
     * Return function handle to the gradient function that calculates
     * the derivative of Istar with respect to beta, at a particular value
     * of beta.
     */
    public ObjectiveFunctionGradient getObjectiveFunctionGradient() {
      return new ObjectiveFunctionGradient(new MultivariateVectorFunction() {
        public double[] value(double[] point) {
          // Gradient calculation with beta
          double beta = Math.exp(point[0]);
    
          try {
            double[][] Q = MatrixUtils.add(invSigma, MatrixUtils.matrixScalarProduct(beta, c));
            double[][] invQ = MatrixUtils.invertSymmPosDefMatrix(Q);
            double[][] d_invQ = MatrixUtils.matrixScalarProduct(-1, MatrixUtils.matrixProduct(invQ, MatrixUtils.matrixProduct(dQ, invQ)));

            double[][] d = MatrixUtils.matrixProduct(b, MatrixUtils.matrixProduct(invQ, MatrixUtils.transpose(b)));

            double[][] dR = MatrixUtils.add( MatrixUtils.add(MatrixUtils.matrixScalarProduct(-1, c), MatrixUtils.matrixScalarProduct(-2*beta, d)) , MatrixUtils.matrixScalarProduct(-beta*beta,  MatrixUtils.matrixProduct(b, MatrixUtils.matrixProduct(d_invQ, MatrixUtils.transpose(b)))));

            double d_term1 = 0.5 * MatrixUtils.trace(MatrixUtils.matrixProduct(invQ, dQ));
            double d_term2 = 0.5 * MatrixUtils.trace(MatrixUtils.matrixProduct(covariance, dR));
            return new double[] {beta*(d_term1 + d_term2 + term3)};

          } catch (Throwable e) {
            e.printStackTrace();
            return new double[] {Double.NaN};
          }
        }
      });
    }

  }

}


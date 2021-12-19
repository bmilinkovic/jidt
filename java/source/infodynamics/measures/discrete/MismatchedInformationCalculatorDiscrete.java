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

package infodynamics.measures.discrete;

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
import java.util.List;
import java.util.ArrayList;

/**
 * <p>
 * Computes the mutual information between the past and future of a discrete
 * multivariate system using a "mismatched decoder" that ignores correlations
 * between system parts. This measure is used as the basis for the
 * decoder-based measure of integrated information (Phi-star).
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
 * 	<li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class MismatchedInformationCalculatorDiscrete
    extends EffectiveMeasureCalculatorDiscrete {

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
  public MismatchedInformationCalculatorDiscrete(int base, int dims) {
    super(base, dims);
  }

  @Override
  public double computeForPartition(List<List<Integer>> partition) {

    ensureMIComputed();

    try {

      calcHelper = new DecodingCalculatorHelper(jointPDF, partition, base);

      // Optimize Istar to calculate mismatched decoding information
      NonLinearConjugateGradientOptimizer opt
          = new NonLinearConjugateGradientOptimizer(
              NonLinearConjugateGradientOptimizer.Formula.POLAK_RIBIERE,
              new SimpleValueChecker(1e-6, 1e-6),
              1e-4, 1e-4, 1);

      PointValuePair optimum
          = opt.optimize(new MaxEval(Integer.MAX_VALUE),
              calcHelper.getObjectiveFunction(),
              calcHelper.getObjectiveFunctionGradient(),
              GoalType.MAXIMIZE,
              new InitialGuess(new double[] { -0.6931 })); // == log(0.5)

      maxBeta = Math.exp(optimum.getPoint()[0]);
      mismatchedInformation = optimum.getValue() / Math.log(2.0);

    } catch (TooManyEvaluationsException e) {
      if (debug) System.out.println("Optimization exceeded maximum number " +
                                    "of evaluations.");
      maxBeta = Double.NaN;
      mismatchedInformation = Double.NaN;

    } catch (Throwable e) {
      e.printStackTrace();
    }

    return (systemInformation - mismatchedInformation);
  }

  /**
   * Compute value of IStar for the most recent data added for the given
   * value of beta.
   */
  public double computeValueForBeta(double beta, List<List<Integer>> partition) throws Exception {
    calcHelper = new DecodingCalculatorHelper(jointPDF, partition, base);
    MultivariateFunction f = calcHelper.getObjectiveFunction().getObjectiveFunction();
    return f.value(new double[] {Math.log(beta)})/Math.log(2.0);
  }

  /**
   * Compute gradient of IStar for the most recent data added for the given
   * value of beta.
   */
  public double computeGradientForBeta(double beta, List<List<Integer>> partition) throws Exception {
    calcHelper = new DecodingCalculatorHelper(jointPDF, partition, base);
    MultivariateVectorFunction f = calcHelper.getObjectiveFunctionGradient().getObjectiveFunctionGradient();
    return (f.value(new double[] {Math.log(beta)})[0]/beta);
  }

  // @Override
  // public double[] computeLocalUsingPreviousObservationsForPartition(int[][] x, List<List<Integer>> partition) throws Exception {

  //   boolean useSpecificSurprise = false;
  //   double[] locals = new double[x.length];

  //   // So that we can use maxBeta;
  //   computeForPartition(partition);
  //   
  //   double[][] Q = computeMismatchedDistribution(partition);

  //   for (int i = 0; i < x.length; i++) {

  //     int pres = MatrixUtils.computeCombinedValuesLittleEndian(x[i], base);

  //     // ========== SPECIFIC INFO ================= 
  //     double specificInfo = 0.0;
  //     if (useSpecificSurprise) {

  //       for (int past = 0; past < nb_states; past++) {
  //         if (jointPDF[past][pres] > 0) {
  //           specificInfo += ( jointPDF[past][pres] / systemPDF[pres] ) * (Math.log(jointPDF[past][pres]) - Math.log(systemPDF[past]) - Math.log(systemPDF[pres]));
  //         }
  //       }

  //     } else {

  //       double ent = 0.0, condEnt = 0.0;
  //       for (int past = 0; past < nb_states; past++) {
  //         if (jointPDF[past][pres] > 0) {
  //           ent -= systemPDF[past] * Math.log(systemPDF[past]);
  //           condEnt -= ( jointPDF[past][pres] / systemPDF[pres] ) * (Math.log(jointPDF[past][pres]) - Math.log(systemPDF[pres]));
  //         }
  //       }

  //       specificInfo = ent - condEnt;
  //     }


  //     // ========== MISMATCHED INFO =================

  //     // First term. We're calculating
  //     // a = \sum_{X_{t-\tau}} p(X_{t-\tau}) q(X_t | X_{t-\tau})^beta
  //     double a = 0.0, localMismatched;
  //     for (int past = 0; past < nb_states; past++) {
  //       a += systemPDF[past] * Math.pow(Q[pres][past], maxBeta);
  //     }
  //     localMismatched = -1 * Math.log(a);

  //     // Second term. We're calculating
  //     // \sum_{X_{t-\tau}} p(X_{t-\tau} | X_t) \log q(X_t | X_{t-\tau})
  //     for (int past = 0; past < nb_states; past++) {
  //       localMismatched += maxBeta*jointPDF[past][pres] * Math.log(Q[pres][past]) / systemPDF[pres];
  //     }

  //     locals[i] = (specificInfo - localMismatched)/Math.log(2.0);

  //   }

  //   return locals;
  // }

  /**
   * Compute the "mismatched decoding" distribution of the system (Q in the
   * notation of Oizumi 2016) for a given system partition.
   *
   * @param partition
   */
  private double[][] computeMismatchedDistribution(List<List<Integer>> partition) throws Exception {

    int nb_parts = partition.size();
    List<int[]> partitions = new ArrayList<int[]>();
    for (int i = 0; i < nb_parts; i++) {
      partitions.add(MatrixUtils.toIntArray(partition.get(i)));
    }

    double[][][] partJoints = new double[nb_parts][][];
    for (int p = 0; p < nb_parts; p++) {
      partJoints[p] = MatrixUtils.marginaliseJointPDF(jointPDF, partitions.get(p), base, dimensions);
    }

    double[][][] partConds = new double[nb_parts][][];
    for (int partIdx = 0; partIdx < nb_parts; partIdx++) {
      partConds[partIdx] = MatrixUtils.normaliseColumns(MatrixUtils.transpose(partJoints[partIdx]));
    }

    // Finally we multiply them together to get the mismatched distribution Q
    double[][] Q = new double[nb_states][nb_states];
    MatrixUtils.fill(Q, 1.0);

    for (int i = 0; i < nb_states; i++) {
      int[] past_ints = MatrixUtils.de2bi(i, dimensions, base);
      for (int j = 0; j < nb_states; j++) {
        int[] pres_ints = MatrixUtils.de2bi(j, dimensions, base);
        for (int partIdx = 0; partIdx < nb_parts; partIdx++) {
          int part_past_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_ints, partitions.get(partIdx)), base);
          int part_pres_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(pres_ints, partitions.get(partIdx)), base);
          Q[j][i] *= partConds[partIdx][part_pres_state][part_past_state];
        }
      }
    }

    return Q;
  }

  /**
   * Helper class
   */
  private class DecodingCalculatorHelper {
    int dimensions, nb_states;
    double[][] Q;
    double[] margProbPast, margProbPres;
    double d_term2;

    /** Constructor. Pre-calculate as many things as possible to make the
     * calculation of I and dI/dbeta faster.
     *
     * @param jointPDF Full joint (i.e. not conditional) past-future
     * probability distribution of the system. jointPDF[i][j] is the
     * probability of observing the system go from state i to state j.
     * 
     * @param partitions List with the variable indices of all subsystems
     * that will be considered for the partition.
     *
     * @param base Number of possible states of each variable.
     */

    public DecodingCalculatorHelper(double[][] jointPDF, List<List<Integer>> partition, int base) throws Exception {
      // Pre-calculations (without beta)
      nb_states = jointPDF.length;
      dimensions = (int) (Math.log(nb_states)/Math.log(base));
      int nb_parts = partition.size();

      List<int[]> partitions = new ArrayList<int[]>();
      for (int i = 0; i < nb_parts; i++) {
        partitions.add(MatrixUtils.toIntArray(partition.get(i)));
      }

      double[][][] partJoints = new double[nb_parts][][];
      for (int p = 0; p < nb_parts; p++) {
        int nb_part_states = (int) Math.pow(base, partitions.get(p).length);
        partJoints[p] = new double[nb_part_states][nb_part_states];
      }
      MatrixUtils.fill(partJoints, 0);

      // According to Oizumi2015, the mismatched distribution q is defined as
      // q(X_t | X_{t-\tau}) = \prod_i p(M^i_t | M^i_{t-\tau}) ,
      // where the product is for all subsystems M in the partition.
      // Calculating Q is a bit involved but we'll get there.
      //
      // Note that jointPDF is the past-present joint distribution, with the
      // first index being past state and the second index the present.

      margProbPast = MatrixUtils.sumRows(jointPDF);
      margProbPres = MatrixUtils.sumColumns(jointPDF);

      // First get the joint of each subsystem by marginalising jointPDF
      for (int i = 0; i < nb_states; i++) {
        int[] past_ints = MatrixUtils.de2bi(i, dimensions, base);
        for (int j = 0; j < nb_states; j++) {
          int[] pres_ints = MatrixUtils.de2bi(j, dimensions, base);
          for (int partIdx = 0; partIdx < nb_parts; partIdx++) {
            int part_past_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_ints, partitions.get(partIdx)), base);
            int part_pres_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(pres_ints, partitions.get(partIdx)), base);
            partJoints[partIdx][part_past_state][part_pres_state] += jointPDF[i][j];
          }
        }
      }

      // Now partJoints[p][i][j] is the probability of observing partition p
      // make a transition from state i to state j. All the entries in
      // partJoints[p] sum to 1.

      // Next we normalise each subsystem's PDF to get conditionals. We have
      // P(past, pres), but to calculate Q we need P(pres|past), so we transpose
      // partJoints and normalise colums.
      double[][][] partConds = new double[nb_parts][][];
      for (int partIdx = 0; partIdx < nb_parts; partIdx++) {
        partConds[partIdx] = MatrixUtils.normaliseColumns(MatrixUtils.transpose(partJoints[partIdx]));
      }

      // Now partConds[p][i][j] is the conditional probability of partition p
      // being in state i given its past state was j. All columns in partConds[p]
      // sum to 1.

      // Finally we multiply them together to get the mismatched distribution Q
      Q = new double[nb_states][nb_states];
      MatrixUtils.fill(Q, 1.0);
      // for (int i = 0; i < nb_states; i++) {
      //   Arrays.fill(Q[i], 1.0);
      // }

      for (int i = 0; i < nb_states; i++) {
        int[] past_ints = MatrixUtils.de2bi(i, dimensions, base);
        for (int j = 0; j < nb_states; j++) {
          int[] pres_ints = MatrixUtils.de2bi(j, dimensions, base);
          for (int partIdx = 0; partIdx < nb_parts; partIdx++) {
            int part_past_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_ints, partitions.get(partIdx)), base);
            int part_pres_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(pres_ints, partitions.get(partIdx)), base);
            Q[j][i] *= partConds[partIdx][part_pres_state][part_past_state];
          }
        }
      }

      // And finally we have our mismatched conditional probability distribution
      // Q. Each column of Q should sum to 1. Now that we have Q, we can do some
      // pre-calculations for IStar that don't depend on beta.

      // We're calculating
      // \sum_{X_t, X_{t-\tau}} p(X_{t-\tau}, X_t) \log q(X_t | X_{t-\tau})
      d_term2 = 0.0;
      for (int past = 0; past < nb_states; past++) {
        for (int pres = 0; pres < nb_states; pres++) {
          if (Q[pres][past] > 0.0) {
            d_term2 += jointPDF[past][pres] * Math.log(Q[pres][past]);
          }
        }
      }

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

            // We're calculating
            // a = \sum_{X_{t-\tau}} p(X_{t-\tau}) q(X_t | X_{t-\tau})^beta
            double[] a = new double[nb_states];
            Arrays.fill(a, 0.0);
            for (int past = 0; past < nb_states; past++) {
              for (int pres = 0; pres < nb_states; pres++) {
                a[pres] += margProbPast[past] * Math.pow(Q[pres][past], beta);
              }
            }

            // We're calculating 
            // term1 = \sum_{X_t} p(X_t) \log a(X_t)
            double term1 = 0.0;
            for (int pres = 0; pres < nb_states; pres++) {
              if (a[pres] > 0.0) {
                term1 += margProbPres[pres] * Math.log(a[pres]);
              }
            }

            double term2 = beta*d_term2;

            return (-term1 + term2);

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

            // We're calculating
            // a = \sum_{X_{t-\tau}} p(X_{t-\tau}) q(X_t | X_{t-\tau})^beta
            double[] a = new double[nb_states];
            Arrays.fill(a, 0.0);
            for (int past = 0; past < nb_states; past++) {
              for (int pres = 0; pres < nb_states; pres++) {
                a[pres] += margProbPast[past] * Math.pow(Q[pres][past], beta);
              }
            }

            // We're calculating
            // b = \sum_{X_{t-\tau}} p(X_{t-\tau}) q(X_t | X_{t-\tau})^\beta * \log q(X_t | X_{t-\tau})
            double[] b = new double[nb_states];
            Arrays.fill(b, 0.0);
            for (int past = 0; past < nb_states; past++) {
              for (int pres = 0; pres < nb_states; pres++) {
                if (Q[pres][past] > 0.0) {
                  b[pres] += margProbPast[past] * Math.pow(Q[pres][past], beta) * Math.log(Q[pres][past]);
                }
              }
            }

            // We're calculating 
            // d_term1 = \sum_{X_t} p(X_t) b(X_t) / a(X_t)
            double d_term1 = 0.0;
            for (int pres = 0; pres < nb_states; pres++) {
              if (b[pres] * a[pres] > 0.0) {
                d_term1 += margProbPres[pres] * b[pres] / a[pres];
              }
            }

            // This extra beta factor is introduced because we're doing the
            // optimization in log-space.
            return new double[] {beta*(-d_term1 + d_term2)};

          } catch (Throwable e) {
            e.printStackTrace();
            return new double[] {Double.NaN};
          }
        }
      });
    }

  }

}

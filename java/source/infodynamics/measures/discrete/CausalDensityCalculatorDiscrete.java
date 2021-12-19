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

import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * <p>
 * Implements causal density for multivariate discrete time-series data.
 * Causal density is defined as the average of the transfer entropies between
 * system parts (all taken with an embedding length of 1). This quantity is
 * computed for multiple partitions of the system (as per
 * {@link infodynamics.measures.discrete.IntegratedMeasureCalculatorDiscrete})
 * and the minimum is returned.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>A. Seth, <a href="https://doi.org/10.1098/rsta.2011.0079">"Causal
 *  density and integrated information as measures of conscious level"</a>,
 *  Philos. Trans. A 369, 2011.</li>
 *
 * 	<li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class CausalDensityCalculatorDiscrete
    extends IntegratedMeasureCalculatorDiscrete {

  /**
   * Constructor.
   * @param base
   * @param tau
   */
  public CausalDensityCalculatorDiscrete(int base, int dimensions) {
      super(base, dimensions);
      PARTITION_SCAN_METHOD = "atomic";
      baseCalculator = new CausalDensityHelper(base, dimensions);
  }


  @Override
  public double computeNormalizationFactor(List<List<Integer>> partition) throws Exception {
    return 1.0;
  }


  /**
   * Mock effective measure to clean up the code. This is the actual causal
   * density for a given partition, these classes are split to recycle more
   * code from IntegratedMeasureCalculatorDiscrete.
   */
  private class CausalDensityHelper extends EffectiveMeasureCalculatorDiscrete {

    /**
     * Constructor.
     */
    public CausalDensityHelper(int base, int dimensions) {
      super(base, dimensions);
    }

    @Override
    public double computeForPartition(List<List<Integer>> partition) throws Exception {

      int nb_parts = partition.size();
      double causalDensity = 0.0;

      // Loop through all pairs of parts and calculate conditional transfer
      // entropy between them
      for (int p2 = 0; p2 < nb_parts; p2++) {

        int[] tgt_part = MatrixUtils.toIntArray(partition.get(p2));
        int nb_states_p2 = (int) Math.pow(base, tgt_part.length);

        // Make new joint PDF marginalising the future of p2 but leaving the
        // past with the whole system
        double[][] reducedJointPDF = new double[nb_states][nb_states_p2];
        MatrixUtils.fill(reducedJointPDF, 0.0);
        for (int future = 0; future < nb_states; future++) {
          int[] future_ints = MatrixUtils.de2bi(future, dimensions, base);
          int part_future_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(future_ints, tgt_part), base);
          for (int past = 0; past < nb_states; past++) {
            reducedJointPDF[past][part_future_state] += jointPDF[past][future];
          }
        }

        for (int p1 = 0; p1 < nb_parts; p1++) {
          if (p1 == p2) {
            continue;
          }

          // Calculating TE from p1 to p2

          int[] src_part = MatrixUtils.toIntArray(partition.get(p1));
          int nb_states_p1 = (int) Math.pow(base, src_part.length);

          // Make new distribution without the past of src
          int[] allExceptSrc = MatrixUtils.allExcept(src_part, dimensions);
          double[][] jointNoSource = new double[nb_states/nb_states_p1][nb_states_p2];
          MatrixUtils.fill(jointNoSource, 0.0);
          for (int past = 0; past < nb_states; past++) {
            int[] past_state_vec = MatrixUtils.de2bi(past, dimensions, base);
            int other_state =
              MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_state_vec, allExceptSrc), base);
            for (int future = 0; future < nb_states_p2; future++) {
              jointNoSource[other_state][future] += reducedJointPDF[past][future];
            }
          }
          double[] pdfNoSource = MatrixUtils.sumRows(jointNoSource);

          double te = 0.0;
          // Start iterating through states to actually compute TE
          for (int past = 0; past < nb_states; past++) {
            int[] past_state_vec = MatrixUtils.de2bi(past, dimensions, base);
            int pastNoSrc =
              MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_state_vec, allExceptSrc), base);
            for (int future = 0; future < nb_states_p2; future++) {
              double joint_prob = reducedJointPDF[past][future];
              if (joint_prob > 0) {
                double logTerm = Math.log(joint_prob/systemPDF[past]) - Math.log(jointNoSource[pastNoSrc][future]/pdfNoSource[pastNoSrc]);
                te += joint_prob * logTerm / Math.log(2.0);
              }
            }
          }

          causalDensity += te;

        } // End for p2
      } // End for p1

      causalDensity /= (double)(nb_parts*nb_parts - nb_parts);

      return causalDensity;

    }

    // public double[] computeLocalUsingPreviousObservationsForPartition(int[][] x, List<List<Integer>> partition) throws Exception {

    //   int nb_parts = partition.size();
    //   double[] locals = new double[x.length];

    //   // Loop through all pairs of parts and calculate conditional transfer
    //   // entropy between them
    //   for (int p2 = 0; p2 < nb_parts; p2++) {

    //     int[] tgt_part = MatrixUtils.toIntArray(partition.get(p2));
    //     int nb_states_p2 = (int) Math.pow(base, tgt_part.length);

    //     // Make new joint PDF marginalising the future of p2 but leaving the
    //     // past with the whole system
    //     double[][] reducedJointPDF = new double[nb_states][nb_states_p2];
    //     MatrixUtils.fill(reducedJointPDF, 0.0);
    //     for (int future = 0; future < nb_states; future++) {
    //       int[] future_ints = MatrixUtils.de2bi(future, dimensions, base);
    //       int part_future_state = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(future_ints, tgt_part), base);
    //       for (int past = 0; past < nb_states; past++) {
    //         reducedJointPDF[past][part_future_state] += jointPDF[past][future];
    //       }
    //     }

    //     for (int p1 = 0; p1 < nb_parts; p1++) {
    //       if (p1 == p2) {
    //         continue;
    //       }

    //       // Calculating TE from p1 to p2

    //       int[] src_part = MatrixUtils.toIntArray(partition.get(p1));
    //       int nb_states_p1 = (int) Math.pow(base, src_part.length);

    //       // Make new distribution without the past of src
    //       int[] allExceptSrc = MatrixUtils.allExcept(src_part, dimensions);
    //       double[][] jointNoSource = new double[nb_states/nb_states_p1][nb_states_p2];
    //       MatrixUtils.fill(jointNoSource, 0.0);
    //       for (int past = 0; past < nb_states; past++) {
    //         int[] past_state_vec = MatrixUtils.de2bi(past, dimensions, base);
    //         int other_state =
    //           MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_state_vec, allExceptSrc), base);
    //         for (int future = 0; future < nb_states_p2; future++) {
    //           jointNoSource[other_state][future] += reducedJointPDF[past][future];
    //         }
    //       }
    //       double[] pdfNoSource = MatrixUtils.sumRows(jointNoSource);

    //       // Stationary pdf of the future of the target
    //       double[] destPDF = MatrixUtils.marginalisePDF(systemPDF, tgt_part, base, dimensions);

    //       // Start iterating through states to actually compute TE
    //       for (int past = 0; past < nb_states; past++) {
    //         int[] past_state_vec = MatrixUtils.de2bi(past, dimensions, base);
    //         int pastNoSrc =
    //           MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_state_vec, allExceptSrc), base);
    //         for (int i = 0; i < x.length; i++) {
    //           int future = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(x[i], tgt_part), base);
    //           double joint_prob = reducedJointPDF[past][future];
    //           if (joint_prob > 0) {
    //             double logTerm = Math.log(joint_prob/systemPDF[past]) - Math.log(jointNoSource[pastNoSrc][future]/pdfNoSource[pastNoSrc]);
    //             locals[i] += (joint_prob/destPDF[future]) * logTerm / Math.log(2.0);
    //           }
    //         }
    //       }

    //     } // End for p2
    //   } // End for p1

    //   locals = MatrixUtils.vectorScalarProduct(locals, 1.0/(double)(nb_parts*nb_parts - nb_parts));

    //   return locals;

    // }

  }


}

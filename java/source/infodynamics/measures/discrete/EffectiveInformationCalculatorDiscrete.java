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
 * 	<li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class EffectiveInformationCalculatorDiscrete
    extends EffectiveMeasureCalculatorDiscrete {

  public EffectiveInformationCalculatorDiscrete(int base, int dimensions) {
    super(base, dimensions);
  }

  public double computeForPartition(List<List<Integer>> partition) {

    ensureMIComputed();

    double rvalue = 0.0;

    try {

      MutualInformationCalculatorDiscrete micd;

      double sum = 0;

      for (int i = 0; i < partition.size(); i++) {
        int[] p = MatrixUtils.toIntArray(partition.get(i));
        int partBase = (int) (Math.pow(base, p.length));
        double[][] partJoint = MatrixUtils.marginaliseJointPDF(jointPDF, p, base, dimensions);
        micd = new MutualInformationCalculatorDiscrete(partBase, partBase, tau);
        sum += micd.computeFromJointPDF(partJoint);
      }

      // Subtract sum of MI of partitions from the MI of system.
      rvalue = systemInformation - sum;

    } catch (Exception e) {
      e.printStackTrace();
    }

    return rvalue;
  }


  // @Override
  // public double[] computeLocalUsingPreviousObservationsForPartition(int[][] x, List<List<Integer>> partition) throws Exception {

  //   double[] locals = new double[x.length];
  //   int nb_parts = partition.size();
  //   int nb_states = jointPDF.length;

  //   // Some testing configuration variables
  //   boolean useMaxEntPrior = false;
  //   boolean usePartPerturbation = false;

  //   // Distribution of the future of the system given the past. If the system
  //   // is ergodic and markovian, this contains all the causal information.
  //   // From here we do all calculations.
  //   double[][] sysFutureGivenPast = new double[nb_states][nb_states];
  //   for (int past = 0; past < nb_states; past++) {
  //     for (int future = 0; future < nb_states; future++) {
  //       sysFutureGivenPast[past][future] = jointPDF[past][future]/systemPDF[past];
  //     }
  //   }

  //   double[] sysPrior;
  //   if (useMaxEntPrior) {
  //     sysPrior = new double[nb_states];
  //     Arrays.fill(sysPrior, 1.0/((double)nb_states));
  //   } else {
  //     sysPrior = systemPDF;
  //   }

  //   // Distribution over the future of the system induced by the maxEnt prior
  //   double[] sysInducedFuture = new double[nb_states];
  //   Arrays.fill(sysInducedFuture, 0.0);
  //   for (int past = 0; past < nb_states; past++) {
  //     for (int future = 0; future < nb_states; future++) {
  //       sysInducedFuture[future] += sysFutureGivenPast[past][future]*sysPrior[past];
  //     }
  //   }

  //   // This is the distribution that goes into the KL, what BalduzziTononi call
  //   // the a posteriori repertoire
  //   double[][] sysPastGivenFuture = new double[nb_states][nb_states];
  //   for (int past = 0; past < nb_states; past++) {
  //     for (int future = 0; future < nb_states; future++) {
  //       sysPastGivenFuture[past][future] = sysFutureGivenPast[past][future]*sysPrior[past]/sysInducedFuture[future];
  //     }
  //   }

  //   // Now we have to calculate the a posteriori repertoire for each part.
  //   // Follows same structure as code above.
  //   double[][][] partPastGivenFuture = new double[nb_parts][][];
  //   for (int pp = 0; pp < nb_parts; pp++) {
  //     int[] p = MatrixUtils.toIntArray(partition.get(pp));
  //     int nb_part_states = (int) Math.pow(base, p.length);

  //     double[][] partFutureGivenPast = new double[nb_part_states][nb_part_states];

  //     if (usePartPerturbation) { // Calculate part's posteriori repertoire as BalduzziTononi2008

  //       double[][] partJointPDF = MatrixUtils.marginaliseJointPDF(jointPDF, p, base, dimensions);
  //       double[] partStationaryPDF = MatrixUtils.sumColumns(partJointPDF);

  //       for (int past = 0; past < nb_part_states; past++) {
  //         for (int future = 0; future < nb_part_states; future++) {
  //           partFutureGivenPast[past][future] = partJointPDF[past][future]/partStationaryPDF[past];
  //         }
  //       }

  //     } else { // Calculate part's posteriori repertoire as BarrettSeth2011

  //       double[][] partFutureGivenSysPast = new double[nb_states][nb_part_states];
  //       MatrixUtils.fill(partFutureGivenSysPast, 0);
  //       for (int past = 0; past < nb_states; past++) {
  //         for (int future = 0; future < nb_states; future++) {
  //           int[] future_sys_state = MatrixUtils.de2bi(future, dimensions, base);
  //           int future_part_state_int = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(future_sys_state, p), base);
  //           partFutureGivenSysPast[past][future_part_state_int] += sysFutureGivenPast[past][future];
  //         }
  //       }
  //       MatrixUtils.fill(partFutureGivenPast, 0);
  //       for (int past = 0; past < nb_states; past++) {
  //         for (int future = 0; future < nb_part_states; future++) {
  //           int[] past_sys_state = MatrixUtils.de2bi(past, dimensions, base);
  //           int past_part_state_int = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_sys_state, p), base);
  //           partFutureGivenPast[past_part_state_int][future] += partFutureGivenSysPast[past][future]*nb_part_states/((double)nb_states);
  //         }
  //       }

  //     }

  //     double[] partPrior;
  //     if (useMaxEntPrior) { // Use maximum entropy prior
  //       partPrior = new double[nb_part_states];
  //       Arrays.fill(partPrior, 1.0/nb_part_states);
  //     } else { // Use the part's stationary distribution as the prior
  //       partPrior = MatrixUtils.marginalisePDF(systemPDF, p, base, dimensions);
  //     }

  //     double[] partInducedFuture = new double[nb_part_states];
  //     Arrays.fill(partInducedFuture, 0.0);
  //     for (int past = 0; past < nb_part_states; past++) {
  //       for (int future = 0; future < nb_part_states; future++) {
  //         partInducedFuture[future] += partFutureGivenPast[past][future]*partPrior[past];
  //       }
  //     }

  //     partPastGivenFuture[pp] = new double[nb_part_states][nb_part_states];
  //     for (int past = 0; past < nb_part_states; past++) {
  //       for (int future = 0; future < nb_part_states; future++) {
  //         partPastGivenFuture[pp][past][future] = partFutureGivenPast[past][future]*partPrior[past]/partInducedFuture[future];
  //       }
  //     }

  //   }

  //   // Once we have all relevant distributions effective info is just a KL divergence
  //   for (int i = 0; i < x.length; i++) {

  //     // Get the future state of the system
  //     int[] future_sys_state = x[i];
  //     int future_sys_state_int = MatrixUtils.computeCombinedValuesLittleEndian(future_sys_state, base);

  //     // Compute the future state of each part
  //     int[] future_part_state_int = new int[nb_parts];
  //     for (int pp = 0; pp < nb_parts; pp++) {
  //       int[] p = MatrixUtils.toIntArray(partition.get(pp));
  //       future_part_state_int[pp] = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(future_sys_state, p), base);
  //     }

  //     // Average over all possible past states
  //     double KL = 0.0;
  //     for (int k = 0; k < nb_states; k++) {
  //       int past_sys_state_int = k;
  //       int[] past_sys_state = MatrixUtils.de2bi(k, dimensions, base);
  //       double sysProbTerm = sysPastGivenFuture[past_sys_state_int][future_sys_state_int];
  //       if (sysProbTerm < 1e-12) { continue; }
  //       double partProbTerm = 1.0;
  //       for (int pp = 0; pp < nb_parts; pp++) {
  //         int[] p = MatrixUtils.toIntArray(partition.get(pp));
  //         int past_part_state_int = MatrixUtils.computeCombinedValuesLittleEndian(MatrixUtils.select(past_sys_state, p), base);
  //         partProbTerm *= partPastGivenFuture[pp][past_part_state_int][future_part_state_int[pp]];
  //       }
  //       KL += sysProbTerm*(Math.log(sysProbTerm) - Math.log(partProbTerm));
  //     }
  //     locals[i] = KL/Math.log(2);
  //   }


  //   return locals;
  // }

}

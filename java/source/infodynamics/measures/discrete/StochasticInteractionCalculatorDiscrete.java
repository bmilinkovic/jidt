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
 * Computes the stochastic interaction of a discrete multivariate system for a
 * given system partition.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>N. Ay, <a href="https://dx.doi.org/10.3390/e17042432">"Information
 *  geometry on complexity and stochastic interaction"</a>, Entropy 17,
 *  2015.</li>
 *
 * 	<li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class StochasticInteractionCalculatorDiscrete
    extends EffectiveMeasureCalculatorDiscrete {

  public StochasticInteractionCalculatorDiscrete(int base, int dimensions) {
    super(base, dimensions);
  }

  @Override
  public double computeForSystem() {

    try {

      int sysBase = (int) (Math.pow(base, dimensions));

      // Calculate Conditional Entropy for whole system.
      MutualInformationCalculatorDiscrete micd =
          new MutualInformationCalculatorDiscrete(sysBase, sysBase, tau);
      EntropyCalculatorDiscrete ecd =
          new EntropyCalculatorDiscrete(sysBase);

      double mi = micd.computeFromJointPDF(jointPDF);
      double h = ecd.computeFromPDF(systemPDF);
      systemInformation = h - mi;

    } catch (Exception e) {
      e.printStackTrace();
    }

    return systemInformation;
  }

  public double computeForPartition(List<List<Integer>> partition) {

    ensureMIComputed();

    double rvalue = 0.0;

    try {

      MutualInformationCalculatorDiscrete micd;
      EntropyCalculatorDiscrete ecd;

      double sum = 0;

      for (int i = 0; i < partition.size(); i++) {
        int[] p = MatrixUtils.toIntArray(partition.get(i));
        int partBase = (int) (Math.pow(base, p.length));
        double[][] partJoint = MatrixUtils.marginaliseJointPDF(jointPDF, p, base, dimensions);
        micd = new MutualInformationCalculatorDiscrete(partBase, partBase, tau);
        double mi = micd.computeFromJointPDF(partJoint);
        double[] partPDF = MatrixUtils.marginalisePDF(systemPDF, p, base, dimensions);
        ecd = new EntropyCalculatorDiscrete(partBase);
        double h = ecd.computeFromPDF(partPDF);
        sum += (h - mi);
      }

      // Subtract sum of MI of partitions from the MI of system.
      rvalue = sum - systemInformation;

    } catch (Exception e) {
      e.printStackTrace();
    }

    return rvalue;
  }

}

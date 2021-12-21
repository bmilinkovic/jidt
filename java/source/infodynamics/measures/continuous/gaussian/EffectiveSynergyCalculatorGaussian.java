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
 * Computes the synergisic information between the past and the future of a
 * multivariate system beyond a given partition.
 * Specifically, it computes the top node of the PID lattice where the sources
 * are the past states of each part and the target is the future of the whole
 * system, using Barrett's Minimum Mutual Information (MMI) PID.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>V. Griffith, <a href="https://arxiv.org/abs/1401.0978">"A Principled
 *  Infotheoretic Phi-like Measure"</a>, arXiv:1401.0978.</li>
 *
 *  <li>A. Barrett, <a href="http://dx.doi.org/10.1103/PhysRevE.91.052802">
 *  "Exploration of synergistic and redundant information sharing in static and
 *  dynamical Gaussian systems"</a>, Physical Review E 91, 2015.</li>
 *
 *  <li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class EffectiveSynergyCalculatorGaussian
   extends EffectiveMeasureCalculatorGaussian {

  /**
   * Constructor.
   */
  public EffectiveSynergyCalculatorGaussian() throws Exception {
    super();
  }

  public double computeForPartition(List<List<Integer>> partition) {
    double ei = 0.0;

    try {

      if (!isMIComputed) {
        computeForSystem();
      }

      MutualInfoCalculatorMultiVariateGaussian micg = new MutualInfoCalculatorMultiVariateGaussian();
      int[] intsNto2N = MatrixUtils.range(dimensions, 2*dimensions - 1);

      // Compute the MI between all-but-one sources (past of each part) with
      // the target (future of whole system)
      double[] sourceMI = new double[partition.size()];
      for (int i = 0; i < partition.size(); i++) {
        int[] p = MatrixUtils.allExcept(MatrixUtils.toIntArray(partition.get(i)), dimensions);
        int[] v = MatrixUtils.append(p, intsNto2N);
        double[][] partCovariance = MatrixUtils.selectRowsAndColumns(laggedCovariance, v, v);

        micg.initialise(p.length, dimensions);
        micg.setCovariance(partCovariance, covFromObservations);
        sourceMI[i] = micg.computeAverageLocalOfObservations();
      }

      // Subtract sum of MI of partitions from the MI of system.
      ei = systemInformation - MatrixUtils.max(sourceMI);

    } catch (Exception e) {
      e.printStackTrace();
    }

    return ei;
  }


}


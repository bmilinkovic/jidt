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

import java.util.List;

/**
 * <p>
 * Computes Griffith's synergy-based integrated information measure, Psi, using
 * Barrett's Minimum Mutual Information (MMI) PID. For more details:
 * {@see infodynamics.measures.discrete.EffectiveSynergyCalculatorDiscrete}
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
 * 	<li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class IntegratedSynergyCalculatorDiscrete
       extends IntegratedMeasureCalculatorDiscrete {

    /**
     * Constructor.
     * @param base
     * @param tau
     */
    public IntegratedSynergyCalculatorDiscrete(int base, int dimensions) {
        super(base, dimensions);
        PARTITION_SCAN_METHOD = "BIPARTITION";
        baseCalculator = new EffectiveSynergyCalculatorDiscrete(base, dimensions);
    }

    @Override
    public double computeNormalizationFactor(List<List<Integer>> partition) throws Exception {
      // psi doesn't need normalization factor
      return 1;
    }

}


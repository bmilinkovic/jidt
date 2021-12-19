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

/**
 * <p>
 * Implements "whole-minus-sum" integrated information, using the system's
 * stationary distribution following Barrett and Seth (2011).
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
public class IntegratedInformationCalculatorDiscrete
       extends IntegratedMeasureCalculatorDiscrete {

    /**
     * Constructor.
     * @param base
     * @param tau
     */
    public IntegratedInformationCalculatorDiscrete(int base, int dimensions) {
        super(base, dimensions);
        PARTITION_SCAN_METHOD = "BIPARTITION";
        baseCalculator = new EffectiveInformationCalculatorDiscrete(base, dimensions);
    }

}

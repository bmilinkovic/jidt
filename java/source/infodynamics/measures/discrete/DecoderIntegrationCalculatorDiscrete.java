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
 * Implements the decoder-based measure of integrated information (Phi-star)
 * for multivariate discrete time-series.
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
public class DecoderIntegrationCalculatorDiscrete
       extends IntegratedMeasureCalculatorDiscrete {

    /**
     * Constructor.
     * @param base
     * @param dimensions
     */
    public DecoderIntegrationCalculatorDiscrete(int base, int dimensions) {
        super(base, dimensions);
        PARTITION_SCAN_METHOD = "ATOMIC";
        baseCalculator = new MismatchedInformationCalculatorDiscrete(base, dimensions);
    }

    @Override
    public double computeNormalizationFactor(List<List<Integer>> partition) throws Exception {
      // phiStar doesn't need normalization factor
      return 1;
    }

}


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

package infodynamics.measures.continuous.kraskov;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * <p>
 * Implements the "Phi-tilde" integrated information measure as defined by
 * Barrett and Seth (2011).
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 *
 * 	<li>Mediano, P., Seth, A., Barrett, A., <a href="http://dx.doi.org/10.3390/e21010017">
 * 	"Measuring integrated	information: Comparison of candidate measures in
 * 	theory and simulation"</a>, Entropy 21(1), 2019.</li> </ul>
 *
 *  <li>Barrett, A., Seth, A., <a href="https://dx.doi.org/10.1371/journal.pcbi.1001052">
 *  "Practical measures of integrated information for time-series data"</a>,
 *  , PLoS Comput Biol 7, 2011.</li>
 * </ul>
 * 
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class IntegratedInteractionCalculatorKraskov
  extends IntegratedMeasureCalculatorKraskov {

  /**
   * Constructor.
   */
  public IntegratedInteractionCalculatorKraskov() {
    super();
    PARTITION_SCAN_METHOD = "BIPARTITION";
    baseCalculator = new StochasticInteractionCalculatorKraskov();
  }

}

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

import java.util.PriorityQueue;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NearestNeighbourSearcher;

/**
 * <p>
 * Computes the stochastic interaction of a continuous multivariate system for
 * a given system partition, using the non-parametric estimator of Kraskov et al.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>Ay, N., <a href="https://dx.doi.org/10.3390/e17042432">"Information
 *  geometry on complexity and stochastic interaction"</a>, Entropy 17,
 *  2015.</li>
 *
 *  <li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 *
 * 	<li>Mediano, P., Seth, A., Barrett, A., <a href="http://dx.doi.org/10.3390/e21010017">
 * 	"Measuring integrated	information: Comparison of candidate measures in
 * 	theory and simulation"</a>, Entropy 21(1), 2019.</li> </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class StochasticInteractionCalculatorKraskov
  extends EffectiveMeasureCalculatorKraskov {

  /**
   * Data structures for fast nearest neighbour searches
   *  representing the marginal spaces
   */
  protected NearestNeighbourSearcher[] rangeSearchersInMarginalsTimeJoint = null;
  protected NearestNeighbourSearcher[] rangeSearchersInMarginalsFutr = null;

  /**
   * Constructor.
   */
  public StochasticInteractionCalculatorKraskov() {
    super();
  }

  @Override
  public void initialise(int dims, int tau) throws Exception {
    rangeSearchersInMarginalsFutr = null;
    rangeSearchersInMarginalsTimeJoint = null;
    super.initialise(dims, tau);
    return;
  }

  @Override
  protected void ensureMarginalSearchersConstructed(List<List<Integer>> partition) throws Exception {

    int nb_parts = partition.size();

    rangeSearchersInMarginalsFutr = new NearestNeighbourSearcher[nb_parts];
    rangeSearchersInMarginalsTimeJoint = new NearestNeighbourSearcher[nb_parts];

    for (int i = 0; i < nb_parts; i++) {
      int[] p = MatrixUtils.toIntArray(partition.get(i));
      double[][] pastPartData = MatrixUtils.selectColumns(pastObservations, p);
      double[][] futrPartData = MatrixUtils.selectColumns(futrObservations, p);

      rangeSearchersInMarginalsFutr[i] = NearestNeighbourSearcher.create(futrPartData);
      rangeSearchersInMarginalsFutr[i].setNormType(normType);

      rangeSearchersInMarginalsTimeJoint[i] = NearestNeighbourSearcher.create(new int[]{p.length, p.length}, new double[][][] {pastPartData, futrPartData});
      rangeSearchersInMarginalsTimeJoint[i].setNormType(normType);
    }

    return;
  }

  @Override
  protected double[] partialComputeFromObservationsSystem(int startTimePoint,
      int numTimePoints, boolean returnLocals, boolean distancesComputed)
      throws Exception {

      throw new Exception("Conditional entropy cannot be computed (without " +
                "substantial bias) with this algorithm.");
  }


  protected double[] partialComputeFromObservationsParts(int startTimePoint,
      int numTimePoints, boolean returnLocals, boolean distancesComputed, List<List<Integer>> partition) throws Exception {
    
    double[] localMi = null;
    if (returnLocals) {
      localMi = new double[numTimePoints];
    }
    
    int nb_parts = partition.size();
    double totalSumF = 0.0;
        
    for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {

      // Get the distance from the kth neighbour in joint (past+future) space.
      // If it hasn't been computed before, compute it and save it for reuse.
      if (!distancesComputed) {
        PriorityQueue<NeighbourNodeData> nnPQ =
            kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
        NeighbourNodeData kthNnData = nnPQ.poll();
        neighbourDistance[t] = kthNnData.distance;
      }
      double eps = neighbourDistance[t];
      
      // Count the number of points whose x distance is less
      //  than kthNnData.distance, for each marginal variable x:
      double sumF = 0.0;

      int n_sys_fut = rangeSearcherFutr.countPointsStrictlyWithinR(
                    t, eps, dynCorrExclTime);

      sumF += (MathsUtils.digamma(n_sys_fut + 1) - digammaN);

      // Count the average number of points within eps for each marginal of each point
      for (int i = 0; i < nb_parts; i++) {

        int n_fut = rangeSearchersInMarginalsFutr[i].countPointsStrictlyWithinR(
                  t, eps, dynCorrExclTime);

        int n_tj = rangeSearchersInMarginalsTimeJoint[i].countPointsStrictlyWithinR(
                     t, eps, dynCorrExclTime);

        sumF -= (MathsUtils.digamma(n_fut + 1) - digammaN);
        sumF += (MathsUtils.digamma(n_tj + 1) - digammaN);
      }

      totalSumF += sumF;

      if (returnLocals) {
        localMi[t - startTimePoint] = digammaK - digammaN - sumF;
      }
    }
    
    // Select what to return:
    if (returnLocals) {
      return localMi;
    } else {
      double[] returnArray = new double[] {digammaK - digammaN - totalSumF/((double) numObservations)};
      return returnArray;
    }
  } 

}


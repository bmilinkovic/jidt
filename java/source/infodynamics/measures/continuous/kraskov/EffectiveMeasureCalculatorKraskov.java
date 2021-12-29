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

import infodynamics.measures.continuous.EffectiveMeasureCalculator;

import java.util.PriorityQueue;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;

import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NearestNeighbourSearcher;

/**
 * <p>
 * Base class for measures of effective information, that quantify some of
 * information "beyond" a bipartition in a multivariate stochastic process,
 * applicable to <code>double[][]</code> time-series data and using the
 * non-parametric estimator of Kraskov et al.
 * </p>
 *
 * <p>
 * Usage is as per the paradigm outlined for {@link EffectiveMeasureCalculator}.
 * </p>
 *
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
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public abstract class EffectiveMeasureCalculatorKraskov
  extends EffectiveMeasureCalculator {

  /**
   * The norm type in use (see {@link #PROP_NORM_TYPE})
   */
  protected int normType = EuclideanUtils.NORM_MAX_NORM;

  /**
   * Whether we use dynamic correlation exclusion
   */
  protected boolean dynCorrExcl = false;

  /**
   * Size of dynamic correlation exclusion window.
   */
  protected int dynCorrExclTime = 0;
    
  /**
   * Property name for the number of K nearest neighbours used in
   * the KSG algorithm (default 4).
   */
  public final static String PROP_K = "k";

  /**
   * Number of K nearest neighbours used in
   * the KSG algorithm (default 4).
   */
  protected int k = 4;

  /**
   * Set of observations of past system states, concatenated for all
   * {@link #addObservations} calls, and used for neighbour searches.
   */
  protected double[][] pastObservations = null;

  /**
   * Set of observations of future system states, concatenated for all
   * {@link #addObservations} calls, and used for neighbour searches.
   */
  protected double[][] futrObservations = null;

  /**
   * Constant for digamma(k), with k the number of nearest neighbours selected
   */
  protected double digammaK;

  /**
   * Constant for digamma(N), with N the number of samples.
   */
  protected double digammaN;
  
  /**
   * Distance from each point to its k-th neighbour in the full joint
   * (X_{t-\tau}, X_t) space
   */
  protected double[] neighbourDistance = null;

  /**
   * Property name for what type of norm to use between data points
   *  for each marginal variable -- Options are defined by 
   *  {@link KdTree#setNormType(String)} and the
   *  default is {@link EuclideanUtils#NORM_MAX_NORM}.
   */
  public final static String PROP_NORM_TYPE = "NORM_TYPE";

  /**
   * Property name for an amount of random Gaussian noise to be
   *  added to the data (default 1e-8 to match the noise order in MILCA toolkit.).
   */
  public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";

  /**
   * Property name for the number of parallel threads to use in the
   *  computation (default is to use all available)
   */
  public static final String PROP_NUM_THREADS = "NUM_THREADS";

  /**
   * Valid property value for {@link #PROP_NUM_THREADS} to indicate
   *  that all available processors should be used. 
   */
  public static final String USE_ALL_THREADS = "USE_ALL";

  /**
   * Property name for a dynamics exclusion time window 
   * otherwise known as Theiler window (see Kantz and Schreiber).
   * Default is 0 which means no dynamic exclusion window.
   */
  public static final String PROP_DYN_CORR_EXCL_TIME = "DYN_CORR_EXCL"; 

  /**
   * Whether to add an amount of random noise to the incoming data
   */
  protected boolean addNoise = true;

  /**
   * Amount of random Gaussian noise to add to the incoming data
   */
  protected double noiseLevel = (double) 1e-8;

  /**
   * Number of parallel threads to use in the computation;
   *  defaults to use all available.
   */
  protected int numThreads = Runtime.getRuntime().availableProcessors(); 

  /**
   * protected k-d tree data structure (for fast nearest neighbour searches)
   *  representing the joint space (including past and future).
   */
  protected KdTree kdTreeJoint = null;

  /**
   * Data structure for fast nearest neighbour searches.
   */
  protected NearestNeighbourSearcher rangeSearcherPast = null;
  protected NearestNeighbourSearcher rangeSearcherFutr = null;


  /**
   * Initialise the calculator. Most importantly, set the neighbout searchers
   * to null so that they get trained again on the new data.
   */
  public void initialise(int dims, int tau) throws Exception {
    super.initialise(dims, tau);
    kdTreeJoint = null;
    rangeSearcherPast = null;
    rangeSearcherFutr = null;
    neighbourDistance = null;
    isMIComputed = false;
    return;
  }


  @Override
  public void finaliseAddObservations() throws Exception {

    pastObservations = new double[numObservations][dimensions];
    futrObservations = new double[numObservations][dimensions];

    int startObservation = 0;
    for (double[][] X : allData) {
      MatrixUtils.arrayCopy(X, 0, 0, pastObservations, startObservation, 0, X.length - tau, dimensions);
      MatrixUtils.arrayCopy(X, tau, 0, futrObservations, startObservation, 0, X.length - tau, dimensions);
			startObservation += X.length - tau;
    }

    addingObservations = false;

		// Add Gaussian noise of std dev noiseLevel to the data if required
		if (addNoise) {
			Random random = new Random();
			for (int r = 0; r < numObservations; r++) {
				for (int c = 0; c < dimensions; c++) {
					pastObservations[r][c] += random.nextGaussian()*noiseLevel;
					futrObservations[r][c] += random.nextGaussian()*noiseLevel;
				}
			}
		}

    // Set the constants:
    digammaK = MathsUtils.digamma(k);
    digammaN = MathsUtils.digamma(numObservations);

    return;
  }


  /**
   * Create nearest-neighbour searchers for the marginal and joint-time
   * distributions of the system.
   *
   * Uses the lastest set of supplied time series. Searchers can be
   * KdTrees or BruteForceSearchers, depending on the size and dimensionality
   * of the time series.
   */
  protected void ensureSystemSearcherConstructed() throws Exception {

    if ((kdTreeJoint == null) || (rangeSearcherPast == null) || (rangeSearcherFutr == null)) {

      kdTreeJoint = new KdTree(new int[] {dimensions, dimensions}, new double[][][] {pastObservations, futrObservations});
      kdTreeJoint.setNormType(normType);

      rangeSearcherPast = NearestNeighbourSearcher.create(pastObservations);
      rangeSearcherPast.setNormType(normType);
      rangeSearcherFutr = NearestNeighbourSearcher.create(futrObservations);
      rangeSearcherFutr.setNormType(normType);
    }

    return;
  }

  /**
   * Make sure that neighbour searchers are constructed for the parts of the
   * system. Since the specific combination of past/future of the parts depends
   * on the measure, this method should be implemented by child classes.
   */
  protected abstract void ensureMarginalSearchersConstructed(List<List<Integer>> partition)
      throws Exception;


  @Override
  public double computeForSystem() throws Exception {

    if (!isMIComputed) {
      ensureSystemSearcherConstructed();
      systemInformation = computeMultiThreaded(true, null);
      isMIComputed = true;
    }

    return systemInformation;
  }

  @Override
  public double computeForPartition(List<List<Integer>> partition) throws Exception {
    ensureSystemSearcherConstructed();
    ensureMarginalSearchersConstructed(partition);

    return computeMultiThreaded(false, partition);
  }


  /**
   * Compute the effective measure or system-level information possibly with
   * multi-threading (depending on the value of {@link #PROP_NUM_THREADS}).
   */
  protected double computeMultiThreaded(boolean computeSystem, List<List<Integer>> partition)
      throws Exception {

    double[] returnValues = null;

    boolean distancesComputed = true;
    if (neighbourDistance == null) {
      // Signal that this run has to compute distances, and allocate memory for them
      distancesComputed = false;
      neighbourDistance = new double[numObservations];
    }

    if (numThreads == 1) {
      // Single-threaded implementation:
      if (computeSystem) {
        returnValues = partialComputeFromObservationsParts(0, numObservations, false, distancesComputed, partition);
      } else {
        returnValues = partialComputeFromObservationsSystem(0, numObservations, false, distancesComputed);
      }

      return returnValues[0];

    } else {
      // Multi-threaded implementation:

      returnValues = new double[numObservations];

      // Distribute the observations to the threads for the parallel processing
      int lTimesteps = numObservations / numThreads; // each thread gets the same amount of data
      int res = numObservations % numThreads; // the first thread gets the residual data
      Thread[] tCalculators = new Thread[numThreads];
      IITKraskovThreadRunner[] runners = new IITKraskovThreadRunner[numThreads];
      for (int t = 0; t < numThreads; t++) {
        int startTime = (t == 0) ? 0 : lTimesteps * t + res;
        int numTimesteps = (t == 0) ? lTimesteps + res : lTimesteps;
        runners[t] = new IITKraskovThreadRunner(this, startTime, numTimesteps, true, distancesComputed, computeSystem, partition);

        tCalculators[t] = new Thread(runners[t]);
        tCalculators[t].start();
      }
      
      // Here, we should wait for the termination of the all threads
      //  and collect their results
      for (int t = 0; t < numThreads; t++) {
        if (tCalculators[t] != null) { // TODO Ipek: can you comment on why we're checking for null here?
          tCalculators[t].join(); 
        }
        // Now we add in the data from this completed thread:
        System.arraycopy(runners[t].getReturnValues(), 0, 
            returnValues, runners[t].myStartTimePoint, runners[t].numberOfTimePoints);
      }

      return MatrixUtils.mean(returnValues);

    }

  }



  /**
   * Compute the local or average part values for a given subset of the time series
   * provided. Each integrated measure must implement this.
   */
  protected abstract double[] partialComputeFromObservationsParts(int startTimePoint,
      int numTimePoints, boolean returnLocals, boolean distancesComputed,
      List<List<Integer>> partition) throws Exception;



  /**
   * Compute the local or average system TDMI values for a given subset of the
   * time series provided.
   */
  protected double[] partialComputeFromObservationsSystem(int startTimePoint,
      int numTimePoints, boolean returnLocals, boolean distancesComputed)
      throws Exception {

    double sumF = 0.0;
    double totalSumF = 0.0;

    double[] localSystemMi = null;
    if (returnLocals) {
      localSystemMi = new double[numTimePoints];
    }

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
      // than kthNnData.distance, for each marginal variable x:
      int n_sys = rangeSearcherPast.countPointsStrictlyWithinR(
                    t, eps, dynCorrExclTime);
      int n_sys_fut = rangeSearcherFutr.countPointsStrictlyWithinR(
                    t, eps, dynCorrExclTime);


      sumF = (MathsUtils.digamma(n_sys + 1) - digammaN) + (MathsUtils.digamma(n_sys_fut + 1) - digammaN);
      totalSumF += sumF;

      if (returnLocals) {
        localSystemMi[t - startTimePoint] = digammaK - digammaN - sumF;
      }
    }

    double[] returnArray = null;
    if (returnLocals) {
      returnArray = localSystemMi;
    } else {
      returnArray = new double[] {digammaK - digammaN - totalSumF/((double) numObservations)};
    }

    return returnArray;
  }


  /**
   * Sets properties for the KSG effective measure calculator.
   *  New property values are not guaranteed to take effect until the next call
   *  to an initialise method. 
   *  
   * <p>Valid property names, and what their
   * values should represent, include:</p>
   * <ul>
   *  <li>{@link #PROP_K} -- number of k nearest neighbours to use in joint kernel space
   *      in the KSG algorithm (default is 4).</li>
   *  <li>{@link #PROP_NORM_TYPE} -- normalization type to apply to 
   *    working out the norms between the points in each marginal space.
   *    Options are defined by {@link KdTree#setNormType(String)} -
   *    default is {@link EuclideanUtils#NORM_MAX_NORM}.</li>
   *  <li>{@link #PROP_DYN_CORR_EXCL_TIME} -- a dynamics exclusion time window,
   *      also known as Theiler window (see Kantz and Schreiber);
   *      default is 0 which means no dynamic exclusion window.</li>
   *  <li>{@link #PROP_NUM_THREADS} -- the integer number of parallel threads
   *    to use in the computation. Can be passed as a string "USE_ALL"
   *      to use all available processors on the machine.
   *      Default is "USE_ALL".</li>
   *  <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
   *    random Gaussian noise to add to
   *      each variable, to avoid having neighbourhoods with artificially
   *      large counts. (We also accept "false" to indicate "0".)
   *  <li>Any valid properties for {@link infodynamics.measures.continuous.IntegratedMeasureCalculator#setProperty(String, String)}.</li>
   * </ul>
   * 
   * <p>Unknown property values are ignored.</p>
   * 
   * @param propertyName name of the property
   * @param propertyValue value of the property
   * @throws Exception for invalid property values
   *
   * @see infodynamics.measures.continuous.IntegratedMeasureCalculator#setProperty()
   */
  public void setProperty(String propertyName, String propertyValue) throws Exception {
    boolean propertySet = true;
    if (propertyName.equalsIgnoreCase(PROP_K)) {
      k = Integer.parseInt(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
      normType = KdTree.validateNormType(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_DYN_CORR_EXCL_TIME)) {
      dynCorrExclTime = Integer.parseInt(propertyValue);
      dynCorrExcl = (dynCorrExclTime > 0);
    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
      if (propertyValue.equals("0") ||
          propertyValue.equalsIgnoreCase("false")) {
        addNoise = false;
        noiseLevel = 0;
      } else {
        addNoise = true;
        noiseLevel = Double.parseDouble(propertyValue);
      }
    } else if (propertyName.equalsIgnoreCase(PROP_NUM_THREADS)) {
      if (propertyValue.equalsIgnoreCase(USE_ALL_THREADS)) {
        numThreads = Runtime.getRuntime().availableProcessors();
      } else { // otherwise the user has passed in an integer:
        numThreads = Integer.parseInt(propertyValue);
      }
    } else {
      // No property was set here
      propertySet = false;
    }
    super.setProperty(propertyName, propertyValue);
    if (debug && propertySet) {
      System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
          " to " + propertyValue);
    }
  }

  @Override
  public String getProperty(String propertyName) throws Exception {
    if (propertyName.equalsIgnoreCase(PROP_K)) {
      return Integer.toString(k);
    } else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
      return KdTree.convertNormTypeToString(normType);
    } else if (propertyName.equalsIgnoreCase(PROP_DYN_CORR_EXCL_TIME)) {
      return Integer.toString(dynCorrExclTime);
    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
      return Double.toString(noiseLevel);
    } else if (propertyName.equalsIgnoreCase(PROP_NUM_THREADS)) {
      return Integer.toString(numThreads);
    } else {
      // try the superclass:
      return super.getProperty(propertyName);
    }
  }

  /**
   * Private class to handle multi-threading of the Kraskov algorithms.
   * Each instance calls partialComputeFromObservations()
   * to compute nearest neighbours for a part of the data.
   */
  private class IITKraskovThreadRunner implements Runnable {
    protected EffectiveMeasureCalculatorKraskov eiCalc;
    protected int myStartTimePoint;
    protected int numberOfTimePoints;
    protected boolean computeLocals;
    protected boolean computeSystem;
    protected boolean distancesComputed;
    protected List<List<Integer>> partition;
    
    protected double[] returnValues = null;
    protected Exception problem = null;
    

    public IITKraskovThreadRunner(
        EffectiveMeasureCalculatorKraskov eiCalc,
        int myStartTimePoint, int numberOfTimePoints,
        boolean computeLocals, boolean distancesComputed,
        boolean computeSystem, List<List<Integer>> partition) {
      this.eiCalc = eiCalc;
      this.myStartTimePoint = myStartTimePoint;
      this.numberOfTimePoints = numberOfTimePoints;
      this.computeLocals = computeLocals;
      this.computeSystem = computeSystem;
      this.distancesComputed = distancesComputed;
      this.partition = partition;
    }

    public IITKraskovThreadRunner(
        EffectiveMeasureCalculatorKraskov eiCalc,
        int myStartTimePoint, int numberOfTimePoints,
        boolean computeLocals, boolean distancesComputed) {
      this(eiCalc, myStartTimePoint, numberOfTimePoints,
          computeLocals, distancesComputed, true, null);
    }
    
    /**
     * Return the values from this part of the data,
     *  or throw any exception that was encountered by the 
     *  thread.
     * 
     * @return the relevant return values from this part of the data
     * @throws Exception an exception previously encountered by this thread.
     */
    public double[] getReturnValues() throws Exception {
      if (problem != null) {
        throw problem;
      }
      return returnValues;
    }
    
    /**
     * Start the thread for the given parameters
     */
    public void run() {
      try {
        if (computeSystem) {
          returnValues = eiCalc.partialComputeFromObservationsSystem(
              myStartTimePoint, numberOfTimePoints, computeLocals, distancesComputed);
        } else {
          returnValues = eiCalc.partialComputeFromObservationsParts(
              myStartTimePoint, numberOfTimePoints, computeLocals, distancesComputed, partition);
        }
      } catch (Exception e) {
        // Store the exception for later retrieval
        problem = e;
        return;
      }
    }
  }
  // end class IITKraskovThreadRunner 


}

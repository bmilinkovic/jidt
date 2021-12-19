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

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * <p>
 * Base class for integrated information measure calculators on discrete
 * (int[][]) data, providing common functionality for user-level measure classes.
 * </p>
 *
 * <p>
 * Usage of the child classes extending this class is intended to follow this
 * paradigm:
 * </p>
 * <ol>
 *    <li>Construct the calculator;</li>
 *    <li>Set properties using {@link #setProperty(String, String)};</li>
 *    <li>Initialise the calculator using {@link #initialise()};</li>
 *    <li>Provide the observations/samples for the calculator
 *        to set up the PDFs, using:
 *      <ul>
 *        <li>{@link #setObservations(int[][])} for calculations
 *          based on single time-series, OR</li>
 *        <li>The following sequence:<ol>
 *            <li>{@link #startAddObservations()}, then</li>
 *            <li>One or more calls to {@link #addObservations(int[][])}, then</li>
 *            <li>{@link #finaliseAddObservations()};</li>
 *          </ol></li>
 *      </ul>
 *    </li>
 *    <li>Calculate the measure (typically with {@link #computeAverageLocalOfObservations}).</li>
 * </ol>
 *
 *  <p>Alternatively, the shortcut method {@link #compute(int[][]} can be
 *  used to directly initialise the calculator and compute the measure on the
 *  provided data.</p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public abstract class IntegratedMeasureCalculatorDiscrete {

  /**
   * Number of available quantised states for each variable
   * (ie binary is base-2).
   */
  protected int base;

  /**
   * Integration time lag <code>\tau</code> (1 by default).
   *
   * Integrated measure will be computed for the joint PDF P(X_t, X_{t+\tau}).
   */
  protected int tau = 1;

  /**
   * Data.
   */
  protected int[][] data;

  /**
   * Running estimator of the *unnormalized* joint (past, present) pdf. The first index runs
   * along past states and the second index along present states.
   */
  protected double[][] jointCount = null;

  /**
   * Running estimator of the joint (past, present) pdf. The first index runs
   * along past states and the second index along present states.
   */
  protected double[][] jointPDF = null;

  /**
   * Running estimator of the stationary system pdf.
   */
  protected double[] systemPDF = null;

  /**
   * Number of observations supplied.
   */
  protected int numObservations = -1;

  /**
   * Number of dimensions of the system.
   */
  protected int dimensions;

  /**
   * Set of all partitions we're going to search through.
   */
  protected List<List<List<Integer>>> partitions;

  /**
   * Stores the indexes of the variables making
   * up the minimum information partition (MIP).
   */
  protected List<List<Integer>> minimumInformationPartition;

  /**
   * Stores the size of the MIP.
   */
  protected int minimumInformationPartitionSize;

  /**
   * Initialise the MIP score to positive infinity so that you return
   * the minimum score as scores for partitions start coming in.
   */
  protected double minimumInformationPartitionScore = Double.POSITIVE_INFINITY;

  /**
   * Stores the value of basic information measure (mutual information
   * or conditional entropy) for the system so that it can be reused.
   */
  protected double systemInformation;

  /**
   * Calculator to base the integrated measure on.
   */
  protected EffectiveMeasureCalculatorDiscrete baseCalculator;

  /**
   * Number of parallel threads to use in the computation;
   *  defaults to use all available. Update: for now use 1.
   */
  // protected int numThreads = Runtime.getRuntime().availableProcessors();
  protected int numThreads = 1;

  /**
   * Whether we're in debug mode.
   */
  protected boolean debug = false;

  /**
   * Property name for integration time lag <code>\tau</code> (1 by default).
   *
   * Integrated measure will be computed for the joint PDF P(X_t, X_{t+\tau}).
   */
  public final static String PROP_TAU = "TAU";

  /**
   * Name of the property used to change the method to calculate local values.
   */
  public final static String PROP_LOCAL_MIP = "LOCAL_MIP";

  /**
   * Name of the property used to change the method to scan partitions.
   */
  public final static String PROP_PARTITION_SCAN_METHOD = "PARTITION_SCAN_METHOD";

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
   * Method used to calculate the MIP for local integration measures. 
   *
   * Possible values are:
   *
   * - "AVERAGE": use the MIP of the average measure to calculate all the locals.
   * - "INDIVIDUAL": calculate the MIP separately for each individual value.
   *   Can be very expensive.
   */
  protected String LOCAL_MIP = "average";

  /**
   * Method used to scan possible partitions. Every class inheriting this one
   * must set a default value.
   *
   * Possible values are:
   *
   * - "ATOMIC": use always the atomic partition of the system, i.e. separate
   *   all variables.
   * - "BIPARTITION": search all possible bipartitions.
   * - "EVEN_BIPARTITIONS": search all bipartitions in which both parts have
   *   the same size (or, for systems with an odd number of variables, in
   *   which there is a difference of 1). When this option is selected, the
   *   effective measure normalisation factor is set to 1.
   * - "ALL": search all possible partitions of any size (*very* expensive!).
   */
  protected String PARTITION_SCAN_METHOD = null;

  /**
   * Whether the measure has been computed with the latest dataset.
   */
  protected boolean isComputed = false;

  /**
   * Whether system MI has been computed with the latest dataset.
   */
  protected boolean isMIComputed = false;

  /**
   * Value of the measure for the latest dataset (if computed).
   */
  protected double integratedMeasure;

  /**
   * Number of possible states of the whole system.
   */
  protected int nb_states;

  /**
   * Whether the calculator accepts observations.
   */
  protected boolean addingObservations = false;

  /**
   * Whether the PDF held by the calculator was estimated from supplied
   * observations.
   */
  protected boolean pdfFromObservations;


  /**
   * Constructor.
   * @param base
   */
  protected IntegratedMeasureCalculatorDiscrete(int base, int dimensions) {
    this.base = base;
    this.dimensions = dimensions;
    this.nb_states = (int) Math.pow(base, dimensions);
  }

  /**
   * Add observations to the calculator. Keep last batch of observations in
   * case it is needed for local measure computation.
   *
   * @param data matrix of observations to be added.
   */
  public void addObservations(int[][] data) throws Exception {

    if (!addingObservations) {
      throw new RuntimeException("Method startAddObservations must be called" +
          " before addObservations");
    }

    if (data[0].length != dimensions) {
      throw new Exception("Data provided do not match initialised number of dimensions.");
    }

    this.data = data;
    numObservations += data.length - tau;
    isComputed = false;
    isMIComputed = false;

    int[] combinedState = MatrixUtils.computeCombinedValuesLittleEndian(data, base);

    for (int i = 0; i < data.length - tau; i++) {
      int pastState = combinedState[i];
      int presState = combinedState[i + tau];
      jointCount[pastState][presState] += 1;
    }

  }

  /**
   * Prepare the calculator to receive subsequent batches of observations.
   * The actual data must be supplied with {@link #addObservations} and then
   * {@link #finaliseAddObservations}.
   */
  public void startAddObservations() {
    numObservations = 0;
    addingObservations = true;
    jointCount = new double[nb_states][nb_states];
    for (int i = 0; i < nb_states; i++) {
      Arrays.fill(jointCount[i], 0);
    }
  }

  /**
   * Finalise the observation addition cycle. To add more observations after
   * this, the user must call {@link #startAddObservations} and start again.
   */
  public void finaliseAddObservations() throws Exception {
    jointPDF = MatrixUtils.matrixScalarProduct(jointCount, 1.0/((double) numObservations));
    systemPDF = MatrixUtils.sumRows(jointPDF);
    addingObservations = false;
    pdfFromObservations = true;

    baseCalculator.setJointPDF(jointPDF, true);
  }

  /**
   * Provide observations to the calculator. Using this method you must supply
   * all observations in one batch. To add observations progressively
   * use {@link #addObservations}.
   */
  public void setObservations(int[][] data) throws Exception {

    if ((data.length - tau) <= dimensions) {
      System.out.printf("Warning. Number of observations %d is smaller " +
          "than number of dimensions %d. Maybe you want to transpose " +
          "your data?", data.length - tau, dimensions);
    }

    startAddObservations();
    addObservations(data);
    finaliseAddObservations();
  }

  /**
   * Manually set the joint system pdf to calculate integration measures.
   */
  public void setJointPDF(double[][] pdf) throws Exception {
    if (Math.abs(MatrixUtils.sum(pdf) - 1) > 1e-8) {
      throw new Exception("Supplied PDF does not sum to 1. Provide a valid PDF.");
    }
    setJointPDF(pdf, true);
  }

  /**
   * Manually set the joint system pdf to calculate integration measures.
   * Allows user to skip test that distribution sums to 1.
   */
  public void setJointPDF(double[][] pdf, boolean skipNormalizationCheck)
      throws Exception {
    if ( (!skipNormalizationCheck) && (Math.abs(MatrixUtils.sum(pdf) - 1) > 1e-8) ) {
      throw new Exception("Supplied PDF does not sum to 1. Provide a valid PDF.");
    }

    jointPDF = pdf;
    systemPDF = MatrixUtils.sumRows(jointPDF);
    pdfFromObservations = false;

    baseCalculator.setJointPDF(pdf, true);
  }

  /**
   * Prepare the calculator for a new dataset or a new property configuration.
   * Must be run always to clean the previous observations and to make the changes
   * in properties effective.
   */
  public void initialise() {
    numObservations = -1;
    isComputed = false;
    isMIComputed = false;
    baseCalculator.initialise();
  }

  /**
   * Compute all even bipartitions for a system of the specified size.
   * Save the bipartitions in a List<List<Integer>>.
   */
  protected void computeEvenBipartitions() {

    List<Integer> list = new ArrayList<Integer>();
    for(int i = 0; i < dimensions; i++) {
      list.add(i);
    }
    List<List<List<Integer>>> allBipartitions = partitionHelper(list, 2);

    partitions = new ArrayList<List<List<Integer>>>();
    int halfDims = dimensions / 2;
    for (List<List<Integer>> p : allBipartitions) {
      int[] p1 = MatrixUtils.toIntArray(p.get(0));
      int[] p2 = MatrixUtils.toIntArray(p.get(1));
      if (p1.length == halfDims || p2.length == halfDims) {
        partitions.add(p);
      }
    }

  }

  /**
   * Compute all possible bipartitions for a system of the specified size.
   * Save the bipartitions in a List<List<Integer>>.
   */
  protected void computePossibleBipartitions() {

    List<Integer> list = new ArrayList<Integer>();
    for(int i = 0; i < dimensions; i++) {
      list.add(i);
    }
    partitions = partitionHelper(list, 2);

  }

  /**
   * Find all m-partitions of set ori.
   *
   * From StackOverflow:
   *   http://stackoverflow.com/questions/20530128/how-to-find-all-partitions-of-a-set
   */
  protected static List<List<List<Integer>>> partitionHelper(List<Integer> ori, int m) {
    List<List<List<Integer>>> ret = new ArrayList<List<List<Integer>>>();
    if(ori.size() < m || m < 1) return ret;

    if(m == 1) {
      List<List<Integer>> partition = new ArrayList<List<Integer>>();
      partition.add(new ArrayList<Integer>(ori));
      ret.add(partition);
      return ret;
    }

    // f(n-1, m)
    List<List<List<Integer>>> prev1 = partitionHelper(ori.subList(0, ori.size() - 1), m);
    for(int i=0; i<prev1.size(); i++) {
      for(int j=0; j<prev1.get(i).size(); j++) {
        // Deep copy from prev1.get(i) to l
        List<List<Integer>> l = new ArrayList<List<Integer>>();
        for(List<Integer> inner : prev1.get(i)) {
          l.add(new ArrayList<Integer>(inner));
        }

        l.get(j).add(ori.get(ori.size()-1));
        ret.add(l);
      }
    }

    List<Integer> set = new ArrayList<Integer>();
    set.add(ori.get(ori.size() - 1));
    // f(n-1, m-1)
    List<List<List<Integer>>> prev2 = partitionHelper(ori.subList(0, ori.size() - 1), m - 1);
    for(int i=0; i<prev2.size(); i++) {
      List<List<Integer>> l = new ArrayList<List<Integer>>(prev2.get(i));
      l.add(set);
      ret.add(l);
    }

    return ret;
  }


  /**
   * Compute all possible partitions of all sizes.
   * Save the partitions in a List<List<Integer>>.
   */
  protected void computePossiblePartitions() throws Exception {
    if (dimensions > 25) {
      throw new Exception("I can't even hold the number of partitions in" +
          "memory without overflowing, your system is way too big. Consider" +
          "using modes BIPARTITION or ATOMIC.");
    } else if (dimensions > 14) {
      System.out.println("This might take longer than the age of the" +
          "universe. Don't hold your breath.");
    }

    List<Integer> list = new ArrayList<Integer>();
    for(int i = 0; i < dimensions; i++) {
      list.add(i);
    }

    partitions = new ArrayList<List<List<Integer>>>();
    for(int i = 2; i <= list.size(); i++) {
      List<List<List<Integer>>> ret = partitionHelper(list, i);
      partitions.addAll(ret);
    }

  }

  /**
   * Compute the integration measure from the previously-supplied samples.
   */
  public double computeAverageLocalOfObservations() throws Exception {

    if ((numObservations < 1) && (jointPDF == null)) {
      throw new Exception("Cannot calculate integration measure without " +
          "providing some data. Use setObservations or setJointPDF first.");
    } else if ((numObservations > 1) && (jointPDF == null)) {
      finaliseAddObservations();
    }

    if (!isComputed) {

      if (!isMIComputed) {
        systemInformation = baseCalculator.computeForSystem();
        isMIComputed = true;
      }

      integratedMeasure = 0.0;
      minimumInformationPartitionScore = Double.POSITIVE_INFINITY;

      if (PARTITION_SCAN_METHOD.equalsIgnoreCase("BIPARTITION")) {
        computePossibleBipartitions();
      } else if (PARTITION_SCAN_METHOD.equalsIgnoreCase("EVEN_BIPARTITIONS")) {
        computeEvenBipartitions();
      } else if (PARTITION_SCAN_METHOD.equalsIgnoreCase("ALL")) {
        computePossiblePartitions();
      } else if (PARTITION_SCAN_METHOD.equalsIgnoreCase("ATOMIC")) {
        List<Integer> list = new ArrayList<Integer>();
        for(int i = 0; i < dimensions; i++) {
          list.add(i);
        }
        partitions = partitionHelper(list, dimensions);
      } else {
        throw new Exception("Unrecognized partition scan method. Options are:" +
            "ALL, BIPARTITION, EVEN_BIPARTITIONS, and ATOMIC.");
      }

      boolean single_thread = PARTITION_SCAN_METHOD.equalsIgnoreCase("ATOMIC") || numThreads == 1;

      if (single_thread) {

        for (List<List<Integer>> partition : partitions) {

          double k;
          if (PARTITION_SCAN_METHOD.equalsIgnoreCase("EVEN_BIPARTITIONS")) {
            k = 1;
          } else {
            k = computeNormalizationFactor(partition);
          }
          double ei = baseCalculator.computeForPartition(partition);

          // If k << 1, it means that one of the partitions has an entropy
          // near 0, which means that it doesn't tell us anything about the
          // rest of the system. Return 0 otherwise return normalised EI.
          double mipScore = (k < 1e-12) ? Double.POSITIVE_INFINITY : ei / k;

          if (mipScore < minimumInformationPartitionScore) {
            minimumInformationPartition = partition;
            minimumInformationPartitionSize = partition.size();
            minimumInformationPartitionScore = mipScore;
            integratedMeasure = ei;
          }

        }
        // End of single-threaded computation

      } else { // go multi-threaded

        // Each thread gets the same amount of data
        int partsPerThread = partitions.size() / numThreads;
        // The first thread gets the residual data
        int res = partitions.size() % numThreads;

        Thread[] tCalculators = new Thread[numThreads];
        PartitionThreadRunner[] runners = new PartitionThreadRunner[numThreads];

        for (int t = 0; t < numThreads; t++) {
          // Fill set with partitions for this thread

          int numParts = (t == 0) ? partsPerThread + res : partsPerThread;
          int start = (t == 0) ? 0 : partsPerThread * t + res;
          List<List<List<Integer>>> parts = partitions.subList(start, start + numParts);

          runners[t] = new PartitionThreadRunner(this, parts);
          tCalculators[t] = new Thread(runners[t]);
          tCalculators[t].start();
        }

        // Now wait for all the threads to come back with the results
        for (int t = 0; t < numThreads; t++) {  
          if (tCalculators[t] != null) {
            // TODO: Check with Joe and Ipek why we check for null here. See MutualInfoMultivariateKraskov.
            tCalculators[t].join();
          }

          double mipScore = runners[t].getBestMIPScore();
          double ei = runners[t].getBestEffectiveMeasure();
          List<List<Integer>> mip = runners[t].getBestMIP();

          if (mipScore < minimumInformationPartitionScore) {
            minimumInformationPartitionScore = mipScore;
            minimumInformationPartition = mip;
            minimumInformationPartitionSize = mip.size();
            integratedMeasure = ei;
          }

        }
        // End multi-threaded computation
      }
    }

    isComputed = true;
    return integratedMeasure;
  }

  /**
   * Compute the local (state-dependent) integration values for the observations
   * supplied.  Note that when multiple sets of data have been provided (via
   * multiple calls to {@link #addObservations}), only the last one is used.
   *
   * <b>NOTE</b>: local versions of integrated information measures are still
   * experimental, and thus not included in the main JIDT release. If you are
   * interested in these, please contact the author.
   */
  public double[] computeLocalOfPreviousObservations() throws Exception {
    if (data == null) {
    throw new Exception("Cannot compute local values of previous observations " +
        "if they have not been set!");
    }

    return computeLocalUsingPreviousObservations(data);
  }

  /**
   * Compute the local (state-dependent) integration values for new observations,
   * based on the statistics computed from the previous observations.
   *
   * <b>NOTE</b>: local versions of integrated information measures are still
   * experimental, and thus not included in the main JIDT release. If you are
   * interested in these, please contact the author.
   *
   * @param x new observations
   */
  public double computeLocalUsingPreviousObservations(int[] x) throws Exception {
    return computeLocalUsingPreviousObservations(new int[][] {x})[0];
  }

  /**
   * Compute the local (state-dependent) integration values for new observations,
   * based on the statistics computed from the previous observations.
   *
   * <b>NOTE</b>: local versions of integrated information measures are still
   * experimental, and thus not included in the main JIDT release. If you are
   * interested in these, please contact the author.
   *
   * @param x new observations
   */
  public double[] computeLocalUsingPreviousObservations(int[][] x) throws Exception {

    double[] locals;

    if (LOCAL_MIP.equalsIgnoreCase("average")) {
      // Use average-case MIP for all locals
      if (!isComputed) {
        computeAverageLocalOfObservations();
      }
      locals = baseCalculator.computeLocalUsingPreviousObservationsForPartition(x, minimumInformationPartition);

    } else if (LOCAL_MIP.equalsIgnoreCase("individual")) {
      // Calculate the MIP separately for all x.
      // Outer loop is over partitions to share computation of normalization
      // factor and covariance KL terms. Then a separate mipScore is kept for
      // each sample

      locals = new double[x.length];

      for (List<List<Integer>> partition : partitions) {
        double[] bestScores = new double[x.length];
        Arrays.fill(bestScores, Double.POSITIVE_INFINITY);

        double[] ei = baseCalculator.computeLocalUsingPreviousObservationsForPartition(x, partition);
        double k = computeNormalizationFactor(partition);

        for (int i = 0; i < x.length; i++) {
          double mipScore = ei[i] / k;
          if (mipScore < bestScores[i]) {
            bestScores[i] = mipScore;
            locals[i] = ei[i];
          }
        }

      }
    } else {
      throw new Exception("Illegal value of property " + PROP_LOCAL_MIP);
    }

    return locals;
  }

  /**
   * Compute partition normalisation factor (typically the minimum entropy among
   * the partitions) a la Balduzzi and Tononi 2008.
   */
  public double computeNormalizationFactor(List<List<Integer>> partition) throws Exception {
    ArrayList<Double> entropies = new ArrayList<Double>();

    for (int i = 0; i < partition.size(); i++) {
      int[] p = MatrixUtils.toIntArray(partition.get(i));
      int nb_part_states = (int) Math.pow(base, p.length);
      double[] partPDF = MatrixUtils.marginalisePDF(systemPDF, p, base, dimensions);
      double auxH = 0;
      for (int j = 0; j < nb_part_states; j++) {
        if (partPDF[j] > 0) {
          auxH += -1 * partPDF[j] * Math.log(partPDF[j]);
        }
      }

      if (auxH < 0) {
        auxH = 0;
      }
      entropies.add(auxH / Math.log(2.0));
    }

    return Collections.min(entropies);
  }

  /**
   * Shortcut method to calculate the measure with the previous settings
   * on one time series.
   */
  public double compute(int[][] observations) throws Exception {
    initialise();
    setObservations(observations);
    return computeAverageLocalOfObservations();
  }

  /**
   * Return the system's time-delayed mutual information if it has been
   * computed for the latest data supplied.
   */
  public double getSystemInformation() throws Exception {
    if (!isComputed) {
      throw new Exception("Last computed value does not correspond to the " +
          "latest data supplied. Run computeAverageLocalOfObservations first.");
    }
    return systemInformation;
  }

  /**
   * Return the size of the minimum information partition if it has been
   * computed for the latest data supplied.
   */
  public int getMinimumInformationPartitionSize() throws Exception {
    if (!isComputed) {
      throw new Exception("Last computed value does not correspond to the " +
          "latest data supplied. Run computeAverageLocalOfObservations first.");
    }
    return minimumInformationPartitionSize;
  }

  /**
   * Return the minimum information partition if it has been
   * computed for the latest data supplied.
   */
  public List<List<Integer>> getMinimumInformationPartition() throws Exception {
    if (!isComputed) {
      throw new Exception("Last computed value does not correspond to the " +
          "latest data supplied. Run computeAverageLocalOfObservations first.");
    }
    return minimumInformationPartition;
  }

  /**
   * Returns observations added to the calculator.
   * @return data
   */
  public int[][] getData() {
    return data;
  }

  /**
   * Return the joint system pdf used by the estimator.
   */
  public double[][] getJointPDF() {
    return jointPDF;
  }

  /**
   * Get the number of samples to be used for the PDFs here 
   * which have been supplied by calls to
   * "setObservations", "addObservations" etc.
   * 
   * @return the number of samples to be used for the PDFs
   */
  public final int getNumObservations() {
    return numObservations;
  }

  /**
   * Set or clear debug mode for extra debug printing to stdout
   * 
   * @param debug new setting for debug mode (on/off)
   */
  public void setDebug(boolean debug) {
    this.debug = debug;
    if (baseCalculator != null) {
      baseCalculator.setDebug(debug);
    }
  }

  /** Sets properties for all discrete integration measures calculators.  New
   * property values are not guaranteed to take effect until the next call to
   * an initialise method.
   *
   * <p>Valid property names, and what their values should represent,
   * include:</p> <ul>
   *  <li>{@link #PROP_TAU} -- timescale across which the
   *  integration is to be calculated. That is, we calculate the integration of
   *  the system between time t and time t+tau.</li>
   *
   *  <li>{@link #PROP_LOCAL_MIP} -- method used to calculate the MIP for local
   *  integration measures. Possible values are:
   *    <ul>
   *      <li>"AVERAGE": use the MIP of the average measure to calculate all the
   *   locals.</li>
   *      <li>"INDIVIDUAL": calculate the MIP separately for each individual value.
   *   Can be very expensive.</li>
   *    </ul>
   *  </li>
   *
   *  <li>{@link #PROP_PARTITION_SCAN_METHOD} -- method used to scan possible
   *  partitions. Possible values are:
   *    <ul>
   *      <li> "ATOMIC": use always the atomic partition of the system, i.e.
   *      separate all variables.</li>
   *      <li>"BIPARTITION": search all possible bipartitions.</li>
   *      <li>"EVEN_BIPARTITIONS": search all bipartitions in which both parts have
   *      the same size (or, for systems with an odd number of variables, in
   *      which there is a difference of 1). When this option is selected, the
   *      effective measure normalisation factor is set to 1.</li>
   *      <li>"ALL": search all possible partitions of any size (<b>very</b> expensive!).</li>
   *   </ul>
   *  </li>
   *  <li>{@link #PROP_NUM_THREADS} -- the integer number of parallel threads
   *  to use in the computation. Can be passed as a string "USE_ALL" to use all
   *  available processors on the machine.  </ul>
   * 
   * <p>Unknown property values are ignored.</p>
   * 
   * @param propertyName name of the property @param propertyValue value of the
   * property @throws Exception for invalid property values
   *
   */
  public void setProperty(String propertyName, String propertyValue)
      throws Exception {
    boolean propertySet = true;
    if (propertyName.equalsIgnoreCase(PROP_TAU)) {
      tau = Integer.parseInt(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_LOCAL_MIP)) {
      LOCAL_MIP = propertyValue;
    } else if (propertyName.equalsIgnoreCase(PROP_PARTITION_SCAN_METHOD)) {
      PARTITION_SCAN_METHOD = propertyValue;
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
    baseCalculator.setProperty(propertyName, propertyValue);
    if (debug && propertySet) {
      System.out.println(this.getClass().getSimpleName() + ": Set property "
            + propertyName + " to " + propertyValue);
    }
  }

  /**
   * Get current property values for the calculator.
   * 
   * <p>Valid property names, and what their
   * values should represent, are the same as those for
   * {@link #setProperty(String, String)}</p>
   * 
   * <p>Unknown property values are responded to with a null return value.</p>
   * 
   * @param propertyName name of the property
   * @return current value of the property
   * @throws Exception for invalid property values
   */
  public String getProperty(String propertyName) throws Exception {
    if (propertyName.equalsIgnoreCase(PROP_TAU)) {
      return Integer.toString(tau);
    } else if (propertyName.equalsIgnoreCase(PROP_PARTITION_SCAN_METHOD)) {
      return PARTITION_SCAN_METHOD;
    } else if (propertyName.equalsIgnoreCase(PROP_NUM_THREADS)) {
      return Integer.toString(numThreads);
    } else if (propertyName.equalsIgnoreCase(PROP_LOCAL_MIP)) {
      return LOCAL_MIP;
    } else {
      // No property was recognised here
      return null;
    }
  }
  

  /**
   * Generate an <b>empirical</b> (bootstrapped) distribution of what the given measure would look like,
   * under a null hypothesis that the source values of our
   * samples had no relation to the destination value (possibly
   * in the context of a conditional value).
   * 
   * <p>See Section II.E "Statistical significance testing" of 
   * the JIDT paper below for a description of how this is done for MI,
   * conditional MI and TE.
   * </p>
   * 
   * <p>Note that if several disjoint time-series have been added 
   * as observations using addObservations methods etc.,
   * then these separate "trials" will be mixed up in the generation
   * of surrogates here.</p>
   * 
   * <p><b>References:</b><br/>
   *  <ul>
   *   <li>J.T. Lizier, "JIDT: An information-theoretic
   *    toolkit for studying the dynamics of complex systems", 2014.</li>
     * </ul>
   * 
   * @param numPermutationsToCheck number of surrogate samples to bootstrap
   *  to generate the distribution.
   * @see "J.T. Lizier, 'JIDT: An information-theoretic
   *    toolkit for studying the dynamics of complex systems', 2014."
   * @return the empirical distribution of measure scores under this null hypothesis.
   */
  public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck)
      throws Exception {
    // Generate the re-ordered indices:
    RandomGenerator rg = new RandomGenerator();
    int[][][] newOrderings = new int[numPermutationsToCheck][][];
    // Generate numPermutationsToCheck * (dimensions-1) permutations of 0 .. data.length-1
    for (int n = 0; n < numPermutationsToCheck; n++) {
      // (Not necessary to check for distinct random perturbations)
      // newOrderings[n] = rg.generateRandomPerturbations(numObservations, dimensions-1);
      // PEDRO NOTE: here numObservations is length of the time series - tau, but if we want to reshuffle we need to reshuffle the whole time series (?)
      newOrderings[n] = rg.generateRandomPerturbations(numObservations + tau, dimensions-1);
    }
    return computeSignificance(newOrderings);
  }

  /**
   * Generate a resampled distribution of what the measure would look like,
   * under a null hypothesis that the individual values of each
   * variable in the 
   * samples have no relation to eachother.
   * That is, we destroy the p(x,y,z,..) correlations, while
   * retaining the p(x), p(y),.. marginals, to check how
   *  significant this measure actually was.
   *  
   * <p>See Section II.E "Statistical significance testing" of 
   * the JIDT paper below for a description of how this is done for MI,
   * we are extending that here.
   * </p>
   * 
   * <p>Note that if several disjoint time-series have been added 
   * as observations using {@link #addObservations(double[])} etc.,
   * then these separate "trials" will be mixed up in the generation
   * of surrogates here.</p>
   * 
   * <p>This method (in contrast to {@link #computeSignificance(int)})
   * allows the user to specify how to construct the surrogates,
   * such that repeatable results may be obtained.</p>
   * 
   * @param newOrderings a specification of how to shuffle the values
   *  to create the surrogates to generate the distribution with. The first
   *  index is the permutation number (i.e. newOrderings.length is the number
   *  of surrogate samples we use to bootstrap to generate the distribution here.)
   *  The second index is the variable number (minus 1, since we don't reorder
   *  the first variable),
   *  Each array newOrderings[i][v] should be an array of length N (where
   *  would be the value returned by {@link #getNumObservations()}),
   *  containing a permutation of the values in 0..(N-1).
   * @return the distribution of surrogate measure values under this null hypothesis.
   * @see "J.T. Lizier, 'JIDT: An information-theoretic
   *    toolkit for studying the dynamics of complex systems', 2014."
   * @throws Exception where the length of each permutation in newOrderings
   *   is not equal to the number N samples that were previously supplied.
   */
  public EmpiricalMeasurementDistribution computeSignificance(int[][][] newOrderings)
      throws Exception {
    
    int numPermutationsToCheck = newOrderings.length;
    if (!isComputed) {
      computeAverageLocalOfObservations();
    }
    
    // Store the real observations and their measure value:
    double actualMeasure = integratedMeasure;
    
    EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);

    int countWhereSurrogateIsMoreSignificantThanOriginal = 0;
    for (int i = 0; i < numPermutationsToCheck; i++) {

      // Compute the measure under this reordering
      double newMeasure = baseCalculator.computeForPartition(minimumInformationPartition, newOrderings[i]);
      measDistribution.distribution[i] = newMeasure;

      if (debug){
        System.out.println("New measure value was " + newMeasure);
      }
      if (newMeasure >= actualMeasure) {
        countWhereSurrogateIsMoreSignificantThanOriginal++;
      }
    }
    
    // Restore the actual measure and the observations
    integratedMeasure = actualMeasure;

    // And return the significance
    measDistribution.pValue = (double) countWhereSurrogateIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
    measDistribution.actualValue = actualMeasure;
    return measDistribution;

  }


  /**
   * Private class to handle multi-threading of the bipartition search.
   * 
   * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
   * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
   */
  private class PartitionThreadRunner implements Runnable {
    protected IntegratedMeasureCalculatorDiscrete iiCalc;
    protected List<List<List<Integer>>> parts;
    
    protected double[] ei, k, mipScore;
    protected double best_ei, best_k, best_mipScore = Double.POSITIVE_INFINITY;
    protected List<List<Integer>> best_mip;

    protected Exception problem = null;
    
    public PartitionThreadRunner(
        IntegratedMeasureCalculatorDiscrete iiCalc,
        List<List<List<Integer>>> parts) {
      this.iiCalc = iiCalc;
      this.parts = parts;
      ei = new double[parts.size()];
      k = new double[parts.size()];
      mipScore = new double[parts.size()];
    }

    public double getBestEffectiveMeasure() throws Exception {
      if (problem != null) { throw problem; }
      return best_ei;
    }

    public double getBestNormalizationFactor() throws Exception {
      if (problem != null) { throw problem; }
      return best_k;
    }

    public double getBestMIPScore() throws Exception {
      if (problem != null) { throw problem; }
      return best_mipScore;
    }

    public List<List<Integer>> getBestMIP() throws Exception {
      if (problem != null) { throw problem; }
      return best_mip;
    }
    
    /**
     * Return the values from this part of the data,
     *  or throw any exception that was encountered by the 
     *  thread.
     * 
     * @return an exception previously encountered by this thread.
     * @throws Exception
     */
    public double[] getEffectiveMeasure() throws Exception {
      if (problem != null) { throw problem; }
      return ei;
    }

    public double[] getNormalizationFactor() throws Exception {
      if (problem != null) { throw problem; }
      return k;
    }

    public double[] getMIPScore() throws Exception {
      if (problem != null) { throw problem; }
      return mipScore;
    }
    
    /**
     * Start the thread for the given partition set.
     */
    public void run() {
      try {
        int i = 0;
        for (List<List<Integer>> p : parts) {
          if (iiCalc.PARTITION_SCAN_METHOD.equalsIgnoreCase("EVEN_BIPARTITIONS")) {
            k[i] = 1;
          } else {
            k[i] = computeNormalizationFactor(p);
          }
          ei[i] = iiCalc.baseCalculator.computeForPartition(p);
          mipScore[i] = (k[i] < 1e-12) ? Double.POSITIVE_INFINITY : ei[i] / k[i];

          if (mipScore[i] < best_mipScore) {
            best_mipScore = mipScore[i];
            best_mip = p;
            best_ei = ei[i];
            best_k = k[i];
          }

          i++;
        }

      } catch (Exception e) {
        // Store the exception for later retrieval
        problem = e;
        return;
      }
    }
  }
  // end class PartitionThreadRunner


}

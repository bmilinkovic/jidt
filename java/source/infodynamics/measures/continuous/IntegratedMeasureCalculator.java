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

package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * Base class for integrated information measure calculators on continuous
 * (double[][]) data, providing common functionality for user-level measure
 * classes.
 * 
 * <p>
 * Usage of the child classes extending this class is intended to follow this
 * paradigm:
 * </p>
 *  <ol>
 *    <li>Construct the calculator;</li>
 *    <li>Set properties using {@link #setProperty(String, String)};</li>
 *    <li>Initialise the calculator using {@link #initialise()} or
 *      {@link #initialise(int)} or {@link #initialise(int, int)};
 *    </li>
 *    <li>Provide the observations/samples for the calculator
 *        to set up the PDFs, using:
 *      <ul>
 *        <li>{@link #setObservations(double[][])} for calculations
 *          based on single time-series, OR</li>
 *        <li>The following sequence:<ol>
 *            <li>{@link #startAddObservations()}, then</li>
 *            <li>One or more calls to {@link #addObservations(double[][])} or
 *              {@link #addObservations(double[][], int, int)}, then</li>
 *            <li>{@link #finaliseAddObservations()};</li>
 *          </ol></li>
 *      </ul>
 *    </li>
 *    <li>Calculate the measure (typically with {@link #computeAverageLocalOfObservations}).</li>
 *  </ol>
 *
 *  <p>Alternatively, the shortcut method {@link #compute(double[][]} can be
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
public abstract class IntegratedMeasureCalculator {

  /**
   * Integration time lag <code>\tau</code> (1 by default).
   *
   * Integrated measure will be computed for the joint PDF P(X_t, X_{t+\tau}).
   */
  protected int tau = 1;

  /**
   * Property name for integration time lag <code>\tau</code> (1 by default).
   *
   * Integrated measure will be computed for the joint PDF P(X_t, X_{t+\tau}).
   */
  public static final String PROP_TAU = "TAU";

  /**
   * Total number of observations supplied.
   * Only valid after {@link #finaliseAddObservations()} is called.
   */
  protected int numObservations = -1;

  /**
   * Number of dimenions of the data
   */
  protected int dimensions = -1;
    
  /**
   * Store the last computed average integration measure
   */
  protected double integratedMeasure = Double.NaN;

  /**
   * Stores the value of basic information measure (mutual information
   * or conditional entropy) for the system so that it can be reused.
   */
  protected double systemInformation = Double.NaN;

  /**
   * Calculator to base the integrated measure on.
   */
  protected EffectiveMeasureCalculator baseCalculator;

  /**
   * Whether we're in debug mode.
   */
  protected boolean debug = false;

  /**
   * Whether the measure has been computed with the latest set of data.
   */
  protected boolean isComputed = false;

  /**
   * Whether the system mutual information has been computed with the latest
   * set of data.
   */
  protected boolean isMIComputed = false;

  /**
   * Whether the calculator accepts observations.
   */
  protected boolean addingObservations = false;

  /**
   * Whether some observations have been added to the calculator since the
   * last initialisation.
   */
  protected boolean observationsAdded = false;

  /**
   * All observations provided.
   */
  protected List<double[][]> allData = null;

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
  protected int minimumInformationPartitionSize = -1;

  /**
   * Initialise the MIP score to positive infinity so that you return
   * the minimum score as scores for partitions start coming in.
   */
  protected double minimumInformationPartitionScore = Double.POSITIVE_INFINITY;

  /**
   * Name of the property used to change the method to calculate local values.
   */
  public final static String PROP_LOCAL_MIP = "LOCAL_MIP";

  /**
   * Name of the property used to change the method to scan partitions.
   */
  public final static String PROP_PARTITION_SCAN_METHOD = "PARTITION_SCAN_METHOD";

  /**
   * Name of the property used to change the method to combine partitions.
   */
  public final static String PROP_PARTITION_COMBINE_METHOD = "PARTITION_COMBINE_METHOD";

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
   * Method used to combine the effective information of all partitions. In
   * most cases this should be kept at "MIN", but the option "MEAN" is also
   * added for completeness.
   */
  protected String PARTITION_COMBINE_METHOD = "MIN";

  /**
   * Constructor.
   */
  public IntegratedMeasureCalculator() {
    // Nothing to do
  }

  /**
   * Initialise the calculator for (re-)use, with the existing (or default) values of parameters
   * Clears any PDFs of previously supplied observations.
   * 
   */
  public void initialise(int dimensions) throws Exception {
    initialise(dimensions, tau);
  }

  /**
   * Initialise the calculator for (re-)use, with some parameters
   * supplied here, and existing (or default) values of other parameters
   * to be used.
   * 
   * @param tau integration time scale
   */
  public void initialise(int dimensions, int tau) throws Exception {
    if (dimensions <= 1) {
      throw new Exception("Integrated measures are not defined for systems " +
                          "with less than two variables.");
    }
    this.dimensions = dimensions;
    this.tau = tau;
    integratedMeasure = 0.0;
    numObservations = 0;
    isComputed = false;
    isMIComputed = false;
    observationsAdded = false;
    allData = new ArrayList<double[][]>();
    baseCalculator.initialise(dimensions, tau);
  }
  
  /**
   * Set properties for the underlying calculator implementation.
   * New property values are not guaranteed to take effect until the next call
   *  to an initialise method. 
   * 
   * <p>Property names defined at the interface level, and what their
   * values should represent, include:</p>
   * <ul>
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
   *
   *  <li>{@link #PROP_PARTITION_COMBINE_METHOD} -- method used to combine the
   *  results from the partition scan. Possible vales are:
   *    <ul>
   *      <li>"MIN" (default): take the minimum across partitions.</li>
   *      <li>"MEAN": take the mean across partitions.</li>
   *    </ul>
   *  </li>
   * </ul>
   *  
   * <p>Unknown property values are ignored.</p>
   * 
   * <p>Note that implementing classes may define additional properties.</p>
   * 
   * @param propertyName name of the property
   * @param propertyValue value of the property
   * @throws Exception for invalid property values
   */
  public void setProperty(String propertyName, String propertyValue) throws Exception {
    boolean propertySet = true;
    if (propertyName.equalsIgnoreCase(PROP_TAU)) {
      tau = Integer.parseInt(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_PARTITION_SCAN_METHOD)) {
      PARTITION_SCAN_METHOD = propertyValue;
    } else if (propertyName.equalsIgnoreCase(PROP_PARTITION_COMBINE_METHOD)) {
      PARTITION_COMBINE_METHOD = propertyValue;
    } else if (propertyName.equalsIgnoreCase(PROP_LOCAL_MIP)) {
      LOCAL_MIP = propertyValue;
    } else {
      // No property was set here
      propertySet = false;
    }
    baseCalculator.setProperty(propertyName, propertyValue);
    if (debug && propertySet) {
      System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
          " to " + propertyValue);
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
    } else if (propertyName.equalsIgnoreCase(PROP_PARTITION_COMBINE_METHOD)) {
      return PARTITION_COMBINE_METHOD;
    } else if (propertyName.equalsIgnoreCase(PROP_LOCAL_MIP)) {
      return LOCAL_MIP;
    } else {
      // No property was recognised here
      return null;
    }
  }
  
  /**
   * Sets a single time-series from which to compute the PDF for the
   * integration measure.
   * Cannot be called in conjunction with other methods for setting/adding
   * observations.
   * 
   * @param observations time-series array of multivariate samples,
   *  where the first index is time and the second is variable.
   */
  public void setObservations(double observations[][]) throws Exception {

    if (observations.length <= dimensions) {
      System.out.printf("Warning. Number of observations %d is smaller " +
          "than number of dimensions %d. Maybe you should have transposed " +
          "your data?\n", observations.length, dimensions);
    }

    startAddObservations();
    addObservations(observations);
    finaliseAddObservations();
  }
  
  /**
   * Signal that we will add in the samples for computing the PDF 
   * from several disjoint time-series or trials via calls to
   * {@link #addObservations(double[][])} etc.
   *
   */
  public void startAddObservations() {
    isComputed = false;
    isMIComputed = false;
    numObservations = 0;
    addingObservations = true;
    observationsAdded = false;
    allData = new ArrayList<double[][]>();
    baseCalculator.startAddObservations();
    return;
  }
  
  /**
   * Add more time-series for the computation of the PDF.
   * The array observations must not be over-written by the user
   *  until after finaliseAddObservations() has been called.
   * 
   * @param observations time-series array of multivariate samples,
   *  where the first index is time and the second is variable.
   */
  public void addObservations(double[][] data) throws Exception {

    if (data[0].length != dimensions) {
        throw new Exception("Supplied observations do not match initialised number of dimensions.");
    }

    this.numObservations += data.length - tau;
    observationsAdded = true;
    allData.add(data);

    baseCalculator.addObservations(data);

    return;
  }

  /**
   * Indicate that no more observations will be added for this run and do
   * the necessary preprocessing.
   *
   * Now observations can be pooled together to calculate sufficient statistics
   * (e.g. means and covariances for Gaussians, nearest-neighbour searchers
   * for Kraskov estimators).
   * 
   * @throws Exception 
   */
  public void finaliseAddObservations() throws Exception {
    baseCalculator.finaliseAddObservations();
  }
  
  /**
   * Compute the integration measure from the previously-supplied samples.
   */
  public double computeAverageLocalOfObservations() throws Exception {

    if (!observationsAdded) {
      throw new Exception("Cannot calculate integration measure without " +
          "providing some data. Give me some observations first.");
    } else if (addingObservations) {
      finaliseAddObservations();
    }

    if (!isComputed) {

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

      for (List<List<Integer>> partition : partitions) {

        double k;
        if (PARTITION_SCAN_METHOD.equalsIgnoreCase("EVEN_BIPARTITIONS")) {
          k = 1;
        } else {
          k = computeNormalizationFactor(partition);
        }
        double ei = baseCalculator.computeForPartition(partition);


        if (PARTITION_COMBINE_METHOD.equalsIgnoreCase("MIN")) {
          double mipScore = (k < 1e-12) ? Double.POSITIVE_INFINITY : ei / k;
          if (mipScore < minimumInformationPartitionScore) {
            minimumInformationPartition = partition;
            minimumInformationPartitionSize = partition.size();
            minimumInformationPartitionScore = mipScore;
            integratedMeasure = ei;
          }
        } else if (PARTITION_COMBINE_METHOD.equalsIgnoreCase("MEAN")) {
          minimumInformationPartition = null;
          minimumInformationPartitionSize = -1;
          integratedMeasure += ei/partitions.size();
        } else {
          throw new Exception("Unrecognized partition combination method. Options are:" +
              "MIN and MEAN.");
        }

      }

      isComputed = true;
    }

    return integratedMeasure;
  }

  /**
   * Compute partition normalisation factor (typically the minimum entropy among
   * the partitions) a la Balduzzi and Tononi 2008.
   */
  public abstract double computeNormalizationFactor(List<List<Integer>> partition) throws Exception;

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
   * Compute the local integration values for each of the
   * previously-supplied samples.
   * 
   * <p>PDFs are computed using all of the previously supplied
   * observations.</p>
   * 
   * <p>If the samples were supplied via a single call such as
   * {@link #setObservations(double[][])},
   * then the return value is a single time-series of local
   * values corresponding to these samples.</p>
   * 
   * <p>Otherwise where disjoint time-series observations were supplied using several 
   *  calls such as {@link addObservations(double[][])}
   *  then the local values for each disjoint observation set will be appended here
   *  to create a single "time-series" return array.</p>
   *  
   * @return the "time-series" of local integration values.
   */
  public double[] computeLocalOfPreviousObservations() throws Exception {
    if (allData.isEmpty()) {
    throw new Exception("Cannot compute local values of previous observations " +
        "if they have not been set!");
    }

    double[] locals = new double[0];
    for (int i = 0; i < allData.size(); i++) {
      double[] l = computeLocalUsingPreviousObservations(allData.get(i));
      locals = MatrixUtils.append(locals, l);
    }

    return locals;
  }

  /**
   * Compute the local integration values for each of the
   * supplied samples in <code>newObservations</code>.
   * 
   * <p>PDFs are computed using all of the previously supplied
   * observations, but not those in <code>newObservations</code> (unless they were
   * some of the previously supplied samples).</p>
   * 
   * @param newObservations time-series for which to compute local integration values
   * @return time-series of local integration values corresponding to newObservations
   * @throws Exception
   */
  public double[] computeLocalUsingPreviousObservations(double[][] x) throws Exception {

    if (!observationsAdded) {
      throw new Exception("Cannot calculate integration measure without " +
          "providing some data. Give me some observations first.");
    } else if (addingObservations) {
      finaliseAddObservations();
    }

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
      }

      minimumInformationPartitionScore = Double.POSITIVE_INFINITY;

      locals = new double[x.length];

      for (List<List<Integer>> partition : partitions) {
        double[] bestScores = new double[x.length];
        Arrays.fill(bestScores, Double.POSITIVE_INFINITY);

        double[] ei = baseCalculator.computeLocalUsingPreviousObservationsForPartition(x, partition);
        double k;
        if (PARTITION_SCAN_METHOD.equalsIgnoreCase("EVEN_BIPARTITIONS")) {
          k = 1;
        } else {
          k = computeNormalizationFactor(partition);
        }

        for (int i = 0; i < x.length; i++) {
          double mipScore = (k == 0) ? Double.POSITIVE_INFINITY : ei[i] / k;
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
   * Calculate local value of a single observation.
   *
   * Shortcut method for {@link computeLocalUsingPreviousObservations(double[][])}
   * with one observation only.
   */
  public double computeLocalUsingPreviousObservations(double[] x) throws Exception {
    return computeLocalUsingPreviousObservations(new double[][] {x})[0];
  }

  // /**
  //  * Generate a bootstrapped distribution of what the AIS would look like,
  //  * under a null hypothesis that the previous <code>k</code> values of our
  //  * samples had no relation to the next value in the time-series.
  //  * 
  //  * <p>See Section II.E "Statistical significance testing" of 
  //  * the JIDT paper below for a description of how this is done for AIS 
  //  * as a mutual information. Basically, the marginal PDFs
  //  * of the past <code>k</code> values, and that of the next value, 
  //  * are preserved, while their joint PDF is destroyed, and the 
  //  * distribution of AIS under these conditions is generated.</p>
  //  * 
  //  * <p>Note that if several disjoint time-series have been added 
  //  * as observations using {@link #addObservations(double[])} etc.,
  //  * then these separate "trials" will be mixed up in the generation
  //  * of surrogates here.</p>
  //  * 
  //  * <p>This method (in contrast to {@link #computeSignificance(int[][])})
  //  * creates <i>random</i> shufflings of the next values for the surrogate AIS
  //  * calculations.</p>
  //  * 
  //  * @param numPermutationsToCheck number of surrogate samples to bootstrap
  //  *  to generate the distribution.
  //  * @return the distribution of AIS scores under this null hypothesis.
  //  * @see "J.T. Lizier, 'JIDT: An information-theoretic
  //  *    toolkit for studying the dynamics of complex systems', 2014."
  //  * @throws Exception
  //  */
  // public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;
  // 
  // /**
  //  * Generate a bootstrapped distribution of what the AIS would look like,
  //  * under a null hypothesis that the previous <code>k</code> values of our
  //  * samples had no relation to the next value in the time-series.
  //  * 
  //  * <p>See Section II.E "Statistical significance testing" of 
  //  * the JIDT paper below for a description of how this is done for AIS 
  //  * as a mutual information. Basically, the marginal PDFs
  //  * of the past <code>k</code> values, and that of the next value, 
  //  * are preserved, while their joint PDF is destroyed, and the 
  //  * distribution of AIS under these conditions is generated.</p>
  //  * 
  //  * <p>Note that if several disjoint time-series have been added 
  //  * as observations using {@link #addObservations(double[])} etc.,
  //  * then these separate "trials" will be mixed up in the generation
  //  * of surrogates here.</p>
  //  * 
  //  * <p>This method (in contrast to {@link #computeSignificance(int)})
  //  * allows the user to specify how to construct the surrogates,
  //  * such that repeatable results may be obtained.</p>
  //  * 
  //  * @param newOrderings a specification of how to shuffle the next values
  //  *  to create the surrogates to generate the distribution with. The first
  //  *  index is the permutation number (i.e. newOrderings.length is the number
  //  *  of surrogate samples we use to bootstrap to generate the distribution here.)
  //  *  Each array newOrderings[i] should be an array of length L being 
  //  *  the value returned by {@link #getNumObservations()},
  //  *  containing a permutation of the values in 0..(L-1).
  //  * @return the distribution of AIS scores under this null hypothesis.
  //  * @see "J.T. Lizier, 'JIDT: An information-theoretic
  //  *    toolkit for studying the dynamics of complex systems', 2014."
  //  * @throws Exception where the length of each permutation in newOrderings
  //  *   is not equal to the number L observations that were previously supplied.
  //  */
  // public EmpiricalMeasurementDistribution computeSignificance(
  //    int[][] newOrderings) throws Exception;

  /**
   * Set or clear debug mode for extra debug printing to stdout
   * 
   * @param debug new setting for debug mode (on/off)
   */
  public void setDebug(boolean debug) {
    this.debug = debug;
    baseCalculator.setDebug(debug);
  }
  
  /** Return the last integration measure calculated in a call to
   * {@link #computeAverageLocalOfObservations()} or
   * {@link #computeLocalOfPreviousObservations()} after the previous
   * {@link #initialise()} call.
   * 
   */
  public double getLastAverage() throws Exception {
    if (!isComputed) {
      throw new Exception("Last computed value does not correspond to the " +
          "latest data supplied. Run computeAverageLocalOfObservations first.");
    }
    return integratedMeasure;
  }

  /**
   * Return the system's time-delayed mutual information if it has been
   * computed for the latest data supplied.
   */
  public double getSystemInformation() throws Exception {
    if (!isMIComputed) {
      systemInformation = baseCalculator.computeForSystem();
    }
    return systemInformation;
  }

  /**
   * Shortcut method to calculate the local measure with the previous settings
   * on one time series.
   */
  public double[] computeLocal(double[][] observations) throws Exception {
    initialise(observations[0].length);
    setObservations(observations);
    return computeLocalOfPreviousObservations();
  }

  /**
   * Return the size of the minimum information partition if it has been
   * computed for the latest data supplied.
   */
  public int getMinimumInformationPartitionSize() throws Exception {
    if (!isComputed) {
      computeAverageLocalOfObservations();
    }
    return minimumInformationPartitionSize;
  }

  /**
   * Return the minimum information partition if it has been
   * computed for the latest data supplied.
   */
  public List<List<Integer>> getMinimumInformationPartition() throws Exception {
    if (!isComputed) {
      computeAverageLocalOfObservations();
    }
    return minimumInformationPartition;
  }

  /**
   * Shortcut method to calculate the measure with the previous settings
   * on one time series.
   */
  public double compute(double[][] observations) throws Exception {
    initialise(observations[0].length);
    setObservations(observations);
    return computeAverageLocalOfObservations();
  }

  /**
   * Get the number of samples to be used for the PDFs here 
   * which have been supplied by calls to
   * {@link #setObservations(double[][])}, {@link #addObservations(double[][])}
   * etc.
   * 
   * <p>Note that the number of samples is not equal to the length of time-series
   * supplied (since we need to accumulate the first
   * <code>tau</code>
   * values of each time-series before taking a sample).
   * </p>
   * 
   * @return the number of samples to be used for the PDFs
   * @throws Exception
   */
  public int getNumObservations() throws Exception {
    return numObservations;
  }

}

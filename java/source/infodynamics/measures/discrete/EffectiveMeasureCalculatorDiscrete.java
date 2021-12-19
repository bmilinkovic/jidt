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
 * Base class for measures of effective information, that quantify some
 * of information "beyond" a bipartition in a multivariate stochastic process.
 * </p>
 *
 * <p>
 * These measures are meant to be used internally by integrated measure
 * calculators. Nonetheless, they are kept public and endowed with the public
 * methods ({@link #initialise}, {@link #addObservations},
 * {@link #computeForPartition(int[])}) in case they are useful for some cases.
 * </p>
 *
 * <p>
 * NOTE: stand-alone usage (outside integrated measure calculators) has not
 * been thoroughly tested.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 *
 */
public abstract class EffectiveMeasureCalculatorDiscrete implements Cloneable {

  /**
   * All observations provided. Kept for computation of local values.
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
   * Number of samples supplied.
   */
  protected int numObservations;

  /**
   * Number of variables in the system.
   */
  protected int dimensions;

  /**
   * Number of available quantised states for each variable
   * (ie binary is base-2).
   */
  protected int base;

  /**
   * Property name for the delay used in the calculation of TDMI and
   * mismatched information.
   */
  public final static String PROP_TAU = "TAU";

  /**
   * Delay used in the calculation of TDMI and effective information measures.
   */
  protected int tau = 1;

  /**
   * Whether the measure has been computed with the latest dataset.
   */
  protected boolean isComputed = false;

  /**
   * Whether system MI has been computed with the latest dataset.
   */
  protected boolean isMIComputed = false;

  /**
   * Stores the value of the Time-Delayed Mutual Information (TDMI)
   * of the system so that it can be reused.
   */
  protected double systemInformation;

  /**
   * Number of possible states of the whole system.
   */
  protected int nb_states;

	/**
	 * Whether to report debug messages or not.
	 */
  protected boolean debug;

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
   * @param tau
   */
  public EffectiveMeasureCalculatorDiscrete(int base, int dimensions) {
    this.base = base;
    this.dimensions = dimensions;
    this.nb_states = (int) Math.pow(base, dimensions);
  }

  /**
   * Adds observations.
   * @param states
   */
  public void addObservations(int[][] data) throws Exception {

    if (!addingObservations) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
    }

    int dataDims = data[0].length;

    if (dataDims != dimensions) {
      throw new Exception("Data provided do not match initialised number of dimensions.");
    }

    numObservations += data.length - tau;
    dimensions = data[0].length;
    this.data = data;

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
    isComputed = false;
    isMIComputed = false;
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
    addingObservations = false;
    setJointPDF(MatrixUtils.matrixScalarProduct(jointCount, 1.0/((double) numObservations)), false);
    pdfFromObservations = true;

    // jointPDF = MatrixUtils.matrixScalarProduct(jointCount, 1.0/((double) numObservations));
    // systemPDF = MatrixUtils.sumRows(jointPDF);
    // addingObservations = false;
  }

  /**
   * Provide observations to the calculator. Using this method you must supply
   * all observations in one batch. To add observations progressively
   * use {@link #addObservations}.
   */
  public void setObservations(int[][] data) throws Exception {

    if (data.length <= dimensions) {
      System.out.printf("Warning. Number of observations %d is smaller " +
          "than number of dimensions %d. Maybe you want to transpose " +
          "your data?", data.length - tau, dimensions);
    }

    startAddObservations();
    addObservations(data);
    finaliseAddObservations();
  }


  /**
   * Sets the joint past-future distribution p(X_{t}, X_{t+tau}) to compute
   * measures without explicitly providing observations.
   *
   * @param pdf
   */
  public void setJointPDF(double[][] pdf) throws Exception {
    setJointPDF(pdf, false);
  }

  /**
   * Sets the joint past-future distribution p(X_{t}, X_{t+tau}) to compute
   * measures without explicitly providing observations.
   * Provides an option to skip checking whether the PDF sums to 1.
   *
   * @param pdf
   * @param skipNormalizationCheck
   */
  public void setJointPDF(double[][] pdf, boolean skipNormalizationCheck)
      throws Exception {
    if ( (!skipNormalizationCheck) && (Math.abs(MatrixUtils.sum(pdf) - 1) > 1e-8) ) {
      throw new Exception("Supplied PDF does not sum to 1. Provide a valid PDF.");
    }

    jointPDF = pdf;
    systemPDF = MatrixUtils.sumRows(jointPDF);
    pdfFromObservations = false;
    isComputed = false;
    isMIComputed = false;
  }


  /**
   * Initialise the calculator.
   */
  public void initialise() {
    isComputed = false;
    isMIComputed = false;
  }

  /**
   * Compute the time-delayed mutual information (TDMI) for the whole system.
   * This is the default upper bound for all integration measures, but can be
   * overriden by children classes if necessary.
   */
  public double computeForSystem() {
    
    if (!isMIComputed) {
      try {

        // Calculate MI for whole system.
        MutualInformationCalculatorDiscrete micd;
        int sysBase = (int) (Math.pow(base, dimensions));
        micd = new MutualInformationCalculatorDiscrete(sysBase, sysBase, tau);
        systemInformation = micd.computeFromJointPDF(jointPDF);
        isComputed = true;

      } catch (Exception e) {
        e.printStackTrace();
      }
    }

    isMIComputed = true;
    return systemInformation;
  }

  /**
   * Make sure system-wide information has been computed. Useful to
   * double-check before other computations.
   */
  public void ensureMIComputed() {
    if (!isMIComputed) {
      try {
        computeForSystem();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }


  /**
   * Compute the effective measure beyond a given partition. The partition must
   * be specified with an int[] array with one element per system variable
   * indicating the part assignment of each system variable. As an example, for
   * a system of 3 variables:
   *
   * - [0, 1, 2] is the atomic partition.
   * - [0, 1, 1] is the partition [[0], [1,2]]
   * - [0, 0, 1] is the partition [[0,1], [2]]
   */
  public double computeForPartition(int[] p) throws Exception {
    // Check that input vector is made of D integers between 0 and D-1
    boolean badInput = false;
    badInput = badInput || (p.length != dimensions);
    for (int i = 0; i < dimensions; i++) {
      badInput = badInput || (p[i] < 0 || p[i] >= dimensions);
    }
    if (badInput) {
      throw new Exception("Bad input. Partition must be specified as an int " +
                          "array with one element per system variable. See " +
                          "Javadocs for details.");
    }

    // Convert partition vector into ArrayList
    List<List<Integer>> pl = new ArrayList<>();
    for (int i = 0; i < dimensions; i++) {
      List<Integer> tmp = new ArrayList<Integer>();
      for (int j = 0; j < dimensions; j++) {
        if (p[j] == i) {
          tmp.add(j);
        }
      }
      if (tmp.size() > 0) {
        pl.add(tmp);
      }
    }

    // Pass the ArrayList to the usual method and return
    return computeForPartition(pl);
  }

  /**
   * Compute the effective measure beyond a given partition.
   */
  public abstract double computeForPartition(List<List<Integer>> p) throws Exception;

	/**
	 * Compute the multi information would look like were all time series
	 *  (bar the first) reordered
	 *  as per the array of time indices in newOrdering.
	 * 
	 * <p>The reordering array contains the reordering for each marginal variable (first index).
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here.</p>
	 *  
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of a shuffled source series here.</p>
	 * 
	 * <p>This method is primarily intended for use in {@link #computeSignificance(int[][])}
	 * however has been made public in case users wish to access it.
	 * </p>
	 * 
	 * @param newOrdering the specific permuted new orderings to use. First index is the variable number
	 *  (minus 1, since we don't reorder the first variable),
	 *  second index is the time step, the value is the reordered time step to use
	 *  for that variable at the given time step.
	 *  The values must be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 *  If null, no reordering is performed.
	 * @return what the average multi-info would look like under this reordering
	 * @throws Exception
	 */
  public double computeForPartition(List<List<Integer>> p1, int[][] newOrdering)
      throws Exception {

		if (newOrdering == null) {
			return computeForPartition(p1);
		}

		// Take a clone of the effective info calculator to compute the measure of
    // the surrogates: (this is a shallow copy, it doesn't make new copies of
    // all the arrays)
		EffectiveMeasureCalculatorDiscrete surrogateCalculator =
				(EffectiveMeasureCalculatorDiscrete) this.clone();
		
    // Generate a new re-ordered source data
    int[][] shuffledData = MatrixUtils.reorderDataForVariables(data, newOrdering);
    surrogateCalculator.initialise();
    surrogateCalculator.setObservations(shuffledData);
    // Compute the measure under this reordering
    return surrogateCalculator.computeForPartition(p1);

  }

  /**
   * Compute the local (state-dependent) integration values for new observations,
   * based on the statistics computed from the previous observations.
   *
   * @param x new observations
   * @param p partition to use for the calculation
   */
  public double computeLocalUsingPreviousObservationsForPartition(int[] x, List<List<Integer>> p) throws Exception {
    double[] res = computeLocalUsingPreviousObservationsForPartition(new int[][] {x}, p);
    return res[0];
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
   * @param p partition to use for the calculation
   */
  // public abstract double[] computeLocalUsingPreviousObservationsForPartition(int[][] x, List<List<Integer>> partition) throws Exception;
  public double[] computeLocalUsingPreviousObservationsForPartition(int[][] x, List<List<Integer>> partition) throws Exception {
    throw new UnsupportedOperationException("Operation not implemented yet. Contact the author for more information.");
  }


  /**
   * Returns observations added to the calculator.
   * @return data
   */
  public int[][] getData() {
    return data;
  }

	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

  /** Sets properties for all discrete effective measures calculators.  New
   * property values are not guaranteed to take effect until the next call to
   * an initialise method.
   *
   * <p>Valid property names, and what their values should represent,
   * include:</p> <ul>
   *  <li>{@link #PROP_TAU} -- timescale across which the
   *  integration is to be calculated. That is, we calculate the integration of
   *  the system between time t and time t+tau.</li>
   *
   * <p>Unknown property values are ignored.</p>
   * 
   * @param propertyName name of the property @param propertyValue value of the
   * property @throws Exception for invalid property values
   *
   */
  public void setProperty(String propertyName, String propertyValue) throws Exception {
    boolean propertySet = true;
    if (propertyName.equalsIgnoreCase(PROP_TAU)) {
      tau = Integer.parseInt(propertyValue);
    } else {
      // No property was set here
      propertySet = false;
    }
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
    } else {
			// No property was recognised here
      return null;
    }
  }

}

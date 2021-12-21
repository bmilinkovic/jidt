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
 * calculators. Nonetheless, they are kept public and endowed with the usual
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
 *  <li>P. Mediano, <a href="http://dx.doi.org/10.3390/e21010017">
 * "Measuring integrated information: Comparison of candidate measures in
 * theory and simulation"</a>, Entropy 21(1), 2019.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 *
 */
public abstract class EffectiveMeasureCalculator {

  /**
   * Delay used in the calculation of TDMI and effective information measures.
   */
  protected int tau = 1;

  /**
   * Stores the value of basic information measure (mutual information
   * or conditional entropy) for the system so that it can be reused.
   */
  protected double systemInformation;

  /**
   * All observations provided.
   */
  protected List<double[][]> allData = null;

  /**
   * Whether system MI has been computed with the latest dataset.
   */
  protected boolean isMIComputed = false;

  /**
   * Number of samples supplied.
   */
  protected int numObservations = -1;

  /**
   * Number of variables in the system.
   */
  protected int dimensions = -1;

  /**
   * All observations provided. Kept for computation of local values.
   */
  protected double[][] data;

  /**
   * Whether we're in debug mode.
   */
  protected boolean debug = false;

  /**
   * Whether the calculator accepts observations.
   */
  protected boolean addingObservations = false;

  /**
   * Property name for the delay used in the calculation of TDMI and
   * mismatched information.
   */
  public final static String PROP_TAU = "TAU";


  /**
   * Constructor.
   */
  public EffectiveMeasureCalculator() {
    // Nothing to do
  }

  /**
   * Provide observations to the calculator. Using this method you must supply
   * all observations in one batch. To add observations progressively
   * use {@link #addObservations}.
   */
  public void setObservations(double[][] data) throws Exception {

    if (data.length - tau <= dimensions) {
      System.out.printf("Warning. Number of observations %d is smaller " +
          "than number of dimensions %d. Maybe you want to transpose " +
          "your data?\n", data.length - tau, dimensions);
    }

    startAddObservations();
    addObservations(data);
    finaliseAddObservations();
  }

  /**
   * Prepare the calculator to receive subsequent batches of observations.
   * The actual data must be supplied with {@link #addObservations} and then
   * {@link #finaliseAddObservations}.
   */
  public void startAddObservations() {
    isMIComputed = false;
    numObservations = 0;
    addingObservations = true;
    allData = new ArrayList<double[][]>();
    return;
  }

  /**
   * Signal that the observations are now all added, PDFs can now be constructed.
   *
   * Each subclass of continuous measures (eg Kraskov, kernel, Gaussian)
   * must implement this.
   */
  public abstract void finaliseAddObservations() throws Exception;

  /**
   * Adds observations.
   * @param data
   */
  public void addObservations(double[][] data) throws Exception {

    if (data[0].length != dimensions) {
        throw new Exception("Supplied observations do not match initialised number of dimensions.");
    }

    this.numObservations += data.length - tau;
    allData.add(data);
    return;
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
    if (dimensions == 1) {
      throw new Exception("Integrated measures are not defined for one-variable systems.");
    }

    this.dimensions = dimensions;
    this.tau = tau;
    numObservations = 0;
    isMIComputed = false;
  }

  /**
   * Compute the time-delayed mutual information (TDMI) for the whole system.
   * This is the default upper bound for all integration measures, but can be
   * overriden by children classes if necessary.
   */
  public abstract double computeForSystem() throws Exception;


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
  // public abstract double[] computeLocalUsingPreviousObservationsForPartition(double[][] x, List<List<Integer>> p) throws Exception;
  public double[] computeLocalUsingPreviousObservationsForPartition(double[][] x, List<List<Integer>> p) throws Exception {
    throw new UnsupportedOperationException("Operation not implemented yet. Contact the author for more information.");
  }

  /**
   * Compute the local (state-dependent) integration values for new observations,
   * based on the statistics computed from the previous observations.
   *
   * @param x new observations
   * @param p partition to use for the calculation
   */
  public double computeLocalUsingPreviousObservationsForPartition(double[] x, List<List<Integer>> p) throws Exception {
    double[] res = computeLocalUsingPreviousObservationsForPartition(new double[][] {x}, p);
    return res[0];
  }

  /**
   * Set or clear debug mode for extra debug printing to stdout
   * 
   * @param debug new setting for debug mode (on/off)
   */
  public void setDebug(boolean debug) {
    this.debug = debug;
  }

  /**
   * Set properties for the underlying calculator implementation.
   * New property values are not guaranteed to take effect until the next call
   *  to an initialise method. 
   * 
   * <p>Property names defined at the interface level, and what their
   * values should represent, include:</p>
   * <ul>
   *  <li>{@link #PROP_TAU} -- integration time scale.</li>
   * </ul>
   *  
   * <p>Unknown property values are ignored.</p>
   * 
   * <p>Note that implementing classes may defined additional properties.</p>
   * 
   * @param propertyName name of the property
   * @param propertyValue value of the property
   * @throws Exception for invalid property values
   */
  public void setProperty(String propertyName, String propertyValue) throws Exception {
    boolean propertySet = true;
    if (propertyName.equalsIgnoreCase(PROP_TAU)) {
      tau = Integer.parseInt(propertyValue);
    } else {
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

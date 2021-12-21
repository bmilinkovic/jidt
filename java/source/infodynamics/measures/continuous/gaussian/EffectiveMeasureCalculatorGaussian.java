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

package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.EffectiveMeasureCalculator;

import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>
 * Base class for measures of effective information, that quantify some of
 * information "beyond" a bipartition in a multivariate stochastic process,
 * applicable to <code>double[][]</code> time-series data and assuming that the
 * probability distribution function for these observations is a multivariate
 * Gaussian distribution.
 * </p>
 *
 * <p>
 * Usage is as per the paradigm outlined for {@link EffectiveMeasureCalculator}.
 * </p>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 *
 */
public abstract class EffectiveMeasureCalculatorGaussian
    extends EffectiveMeasureCalculator {

  /**
   * Time-lagged covariance of the system, i.e. cov({X_t, X_{t+\tau}}).
   */
  protected double[][] laggedCovariance = null;

  /**
   * Stationary covariance of the system, i.e. cov(X).
   */
  protected double[][] covariance = null;

  /**
   * Mean of each variable in the system.
   */
  protected double[] means = null;

  /**
   * Whether the covariance in the calculator was obtained from data or
   * supplied directly by the user.
   */
  protected boolean covFromObservations;


  /**
   * Constructor.
   */
  public EffectiveMeasureCalculatorGaussian() {
    super();
  }

  /**
   * {@inheritdoc}
   */
  public void initialise(int dimensions) throws Exception {
    initialise(dimensions, tau);
  }

  /**
   * {@inheritdoc}
   */
  public void initialise(int dimensions, int tau) throws Exception {
    super.initialise(dimensions, tau);
    means = null;
    covariance = null;
    laggedCovariance = null;
    return;
  }

  /**
   * {@inheritdoc}
   */
  public void startAddObservations() {
    super.startAddObservations();
    means = null;
    covariance = null;
    laggedCovariance = null;
    return;
  }

  /**
   * {@inheritdoc}
   */
  public void finaliseAddObservations() {

    addingObservations = false;
    double[] sumData = new double[dimensions];
    double[][] sumDataSquared = new double[2*dimensions][2*dimensions];
    MatrixUtils.fill(sumData, 0.0);
    MatrixUtils.fill(sumDataSquared, 0.0);

    // Calculate means
    for (int i = 0; i < allData.size();  i++) {
      double[][] X = allData.get(i);
      for (int t = tau; t < X.length; t++) {
        for (int j = 0; j < dimensions; j++) {
          sumData[j] += X[t][j];
        }
      }
    }
    means = new double[dimensions];
    for (int j = 0; j < dimensions; j++) {
      means[j] = sumData[j]/((double) numObservations);
    }

    // Calculate lagged covariance
    for (int i = 0; i < allData.size();  i++) {
      double[][] X = allData.get(i);
      for (int t = tau; t < X.length; t++) {
        for (int j = 0; j < dimensions; j++) {
          for (int k = 0; k < dimensions; k++) {
            sumDataSquared[j][k] += (X[t-tau][j] - means[j])*(X[t-tau][k] - means[k]);
          }
        }
        for (int j = 0; j < dimensions; j++) {
          for (int k = dimensions; k < 2*dimensions; k++) {
            sumDataSquared[j][k] += (X[t-tau][j] - means[j])*(X[t][k%dimensions] - means[k%dimensions]);
          }
        }
        for (int j = dimensions; j < 2*dimensions; j++) {
          for (int k = 0; k < dimensions; k++) {
            sumDataSquared[j][k] += (X[t][j%dimensions] - means[j%dimensions])*(X[t-tau][k] - means[k]);
          }
        }
        for (int j = dimensions; j < 2*dimensions; j++) {
          for (int k = dimensions; k < 2*dimensions; k++) {
            sumDataSquared[j][k] += (X[t][j%dimensions] - means[j%dimensions])*(X[t][k%dimensions] - means[k%dimensions]);
          }
        }
      }
    }

    // Calculate covariances
    laggedCovariance = new double[2*dimensions][2*dimensions];
    for (int j = 0; j < 2*dimensions; j++) {
      for (int k = 0; k < 2*dimensions; k++) {
        laggedCovariance[j][k] = sumDataSquared[j][k]/((double) numObservations - 1);
      }
    }

    int[] ints1toN = MatrixUtils.range(0, dimensions - 1);
    covariance = MatrixUtils.selectRowsAndColumns(laggedCovariance, ints1toN, ints1toN);
    covFromObservations = true;

    return;
  }

  /**
   * <p>Set the full time-lagged covariance matrix and mean vector of the system.</p>
   *
   * <p>This means that {@link #tau} will be ignored.</p>
   *
   * @param lcov 
   * @param mu
   */
  public void setLaggedCovarianceAndMeans(double[][] lcov, double[] mu) throws Exception {
    if (mu.length != dimensions) {
      throw new Exception("Supplied means don't match initialised number of dimensions");
    }
    means = mu;
    setLaggedCovariance(lcov, false);
    return;
  }

  /**
   * <p>Set the full time-lagged covariance matrix of the system.</p>
   *
   * <p>This means that {@link #tau} will be ignored.</p>
   *
   * @param lcov 
   * @param determinedFromObservations
   */
  public void setLaggedCovariance(double[][] lcov,
    boolean determinedFromObservations)
    throws IllegalArgumentException, Exception {

      if (lcov.length != lcov[0].length) {
          throw new IllegalArgumentException("Covariance matrices must be square");
      }
      if (lcov.length != 2*dimensions) {
          throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
      }

      if (!determinedFromObservations) {
        numObservations = -1;
      }
      isMIComputed = false;
      this.laggedCovariance = lcov;
      int[] ints1toN = MatrixUtils.range(0, dimensions - 1);
      this.covariance = MatrixUtils.selectRowsAndColumns(lcov, ints1toN, ints1toN);
      this.covFromObservations = determinedFromObservations;
  }

  /**
   * <p>Set the full time-lagged covariance matrix of the system.</p>
   *
   * <p>This means that {@link #tau} will be ignored.</p>
   *
   * @param lcov 
   */
  public void setLaggedCovariance(double[][] lcov) throws Exception {
    setLaggedCovariance(lcov, false);
  }

  @Override
  public double computeForSystem() throws Exception {

    if (!isMIComputed) {
      MutualInfoCalculatorMultiVariateGaussian  micg;

      // Calculate MI for whole system.
      micg = new MutualInfoCalculatorMultiVariateGaussian();
      micg.initialise(dimensions, dimensions);
      micg.setCovariance(laggedCovariance, covFromObservations);
      systemInformation = micg.computeAverageLocalOfObservations();

      isMIComputed = true;
    }

    return systemInformation;

  }

  /**
   * Returns covariance of the stationary distribution
   * @return covariance
   */
  public double[][] getCovariance() {
      return covariance;
  }

  /**
   * Returns the time-lagged covariance of the system, i.e. cov({X_t, X_{t+\tau}}).
   */
  public double[][] getLaggedCovariance() {
      return laggedCovariance;
  }

  /**
   * Returns the mean state of the system.
   */
  public double[] getMeans() {
    return means;
  }

  /**
   * Returns the time-lagged covariance cov({M_t, M_{t+\tau}) of a part of the
   * system M, comprising the variables indexed by p.
   *
   * @param p
   *
   */
  protected double[][] getLaggedCovarianceForPartition(int[] p) {
      int[] q = new int[2*p.length];
      for (int i = 0; i < p.length; i++) {
          q[i] = p[i];
          q[i + p.length] = p[i] + dimensions;
      }

      return MatrixUtils.selectRowsAndColumns(laggedCovariance, q, q);

  }

}


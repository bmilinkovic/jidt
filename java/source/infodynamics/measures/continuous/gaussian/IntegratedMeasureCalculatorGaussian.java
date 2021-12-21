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

import infodynamics.measures.continuous.IntegratedMeasureCalculator;

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Collections;

/**
 * <p>
 * Base class for all integrated information measure calculators on real-valued
 * (double[][]) data under normality assumptions, providing common functionality
 * for user-level measure classes.
 * </p>
 *
 * <p>
 * Usage is as per the paradigm outlined for {@link IntegratedMeasureCalculator}.
 * </p>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public abstract class IntegratedMeasureCalculatorGaussian
    extends IntegratedMeasureCalculator {

  /**
   * Constructor.
   */
  protected IntegratedMeasureCalculatorGaussian() {
    super();
  }

  /**
   * <p>Set the full time-lagged covariance matrix and mean vector of the system.</p>
   *
   * <p>This means that {@link #tau} will be ignored.</p>
   *
   * @param lcov 
   * @param mu
   */
  public void setLaggedCovarianceAndMeans(double[][] cov, double[] mu) throws Exception {
    ((EffectiveMeasureCalculatorGaussian) baseCalculator).setLaggedCovarianceAndMeans(cov, mu);
    this.isComputed = false;
    this.isMIComputed = false;
    this.observationsAdded = true;
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
    ((EffectiveMeasureCalculatorGaussian) baseCalculator).setLaggedCovariance(lcov, determinedFromObservations);
    this.isComputed = false;
    this.isMIComputed = false;
    this.observationsAdded = true;
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
  public double computeNormalizationFactor(List<List<Integer>> partition) throws Exception {

    ArrayList<Double> entropies = new ArrayList<Double>();
    double[][] covariance = ((EffectiveMeasureCalculatorGaussian) baseCalculator).getCovariance();

    for (int i = 0; i < partition.size(); i++) {
      int[] p = MatrixUtils.toIntArray(partition.get(i));
      EntropyCalculatorMultiVariateGaussian ecmg = new EntropyCalculatorMultiVariateGaussian();
      ecmg.initialise(p.length);
      ecmg.setCovariance(MatrixUtils.selectRowsAndColumns(covariance, p, p));
      entropies.add(ecmg.computeAverageLocalOfObservations());
    }

    return Collections.min(entropies);
  }

  /**
   * Returns laggedCovariance added to the calculator.
   * @return laggedCovariance
   */
  public double[][] getLaggedCovariance() {
    return ((EffectiveMeasureCalculatorGaussian) baseCalculator).getLaggedCovariance();
  }

  /**
   * Returns stationary covariance added to the calculator.
   * @return covariance
   */
  public double[][] getCovariance() {
    return ((EffectiveMeasureCalculatorGaussian) baseCalculator).getCovariance();
  }

  /**
   * Returns means added to the calculator.
   * @return means
   */
  public double[] getMeans() {
    return ((EffectiveMeasureCalculatorGaussian) baseCalculator).getMeans();
  }

}


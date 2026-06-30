/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.apache.commons.math3.ode.nonstiff;


import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.util.FastMath;

/**
 * This class implements the common part of all fixed step Runge-Kutta
 * integrators for Ordinary Differential Equations.
 *
 * <p>These methods are explicit Runge-Kutta methods, their Butcher
 * arrays are as follows :
 * <pre>
 *    0  |
 *   c2  | a21
 *   c3  | a31  a32
 *   ... |        ...
 *   cs  | as1  as2  ...  ass-1
 *       |--------------------------
 *       |  b1   b2  ...   bs-1  bs
 * </pre>
 * </p>
 *
 * @see EulerIntegrator
 * @see ClassicalRungeKuttaIntegrator
 * @see GillIntegrator
 * @see MidpointIntegrator
 * @version $Id$
 * @since 1.2
 */

public abstract class RungeKuttaIntegrator extends AbstractIntegrator {

    /** Time steps from Butcher array (without the first zero). */
    private final double[] c;

    /** Internal weights from Butcher array (without the first empty row). */
    private final double[][] a;

    /** External weights for the high order method from Butcher array. */
    private final double[] b;

    /** Prototype of the step interpolator. */
    private final RungeKuttaStepInterpolator prototype;

    /** Integration step. */
    private final double step;

  /** Simple constructor.
   * Build a Runge-Kutta integrator with the given
   * step. The default step handler does nothing.
   * @param name name of the method
   * @param c time steps from Butcher array (without the first zero)
   * @param a internal weights from Butcher array (without the first empty row)
   * @param b propagation weights for the high order method from Butcher array
   * @param prototype prototype of the step interpolator to use
   * @param step integration step
   */
  protected RungeKuttaIntegrator(final String name,
                                 final double[] c, final double[][] a, final double[] b,
                                 final RungeKuttaStepInterpolator prototype,
                                 final double step) {
    super(name);
    this.c          = c;
    this.a          = a;
    this.b          = b;
    this.prototype  = prototype;
    this.step       = FastMath.abs(step);
  }

  /** {@inheritDoc} */
  @Override
  public void integrate(final ExpandableStatefulODE equations, final double t)
      throws NumberIsTooSmallException, DimensionMismatchException,
             MaxCountExceededException, NoBracketingException {

    sanityChecks(equations, t);
    setEquations(equations);
    final boolean forward = t > equations.getTime();

    // create some internal working arrays
    final double[] y0      = equations.getCompleteState();
    final double[] y       = y0.clone();
    final int stages       = c.length + 1;
    final double[][] yDotK = initWorkingArrays(stages, y.length);
    final double[] yTmp    = y.clone();
    final double[] yDotTmp = new double[y.length];

    // set up an interpolator sharing the integrator arrays
    final RungeKuttaStepInterpolator interpolator = initInterpolator(equations, yTmp, yDotK, forward);

    // set up integration control objects
    stepStart = equations.getTime();
    stepSize  = forward ? step : -step;
    initIntegration(equations.getTime(), y0, t);

    // main integration loop
    isLastStep = false;
    do {

      interpolator.shift();

      // first stage
      computeDerivatives(stepStart, y, yDotK[0]);

      // next stages
      computeStageDerivatives(y, yDotK, yTmp, stages);

      // estimate the state at the end of the step
      estimateState(y, yDotK, yTmp, stages);

      // discrete events handling
      stepStart = acceptAndResetStep(interpolator, yTmp, y, yDotK, yDotTmp, stages, t);

      // prepare next step
      prepareNextStep(interpolator, forward, t);

    } while (!isLastStep);

    // dispatch results
    equations.setTime(stepStart);
    equations.setCompleteState(y);

    resetIntegrationState();

  }

  /**
   * Initialize working arrays.
   *
   * @param stages number of stages
   * @param length length of state vector
   * @return initialized nested double arrays
   */
  private double[][] initWorkingArrays(final int stages, final int length) {
    final double[][] yDotK = new double[stages][];
    for (int i = 0; i < stages; ++i) {
      yDotK[i] = new double[length];
    }
    return yDotK;
  }

  /**
   * Initialize step interpolator.
   *
   * @param equations equations to integrate
   * @param yTmp temporary state array
   * @param yDotK derivatives array
   * @param forward forward integration flag
   * @return initialized interpolator
   */
  private RungeKuttaStepInterpolator initInterpolator(final ExpandableStatefulODE equations,
                                                       final double[] yTmp,
                                                       final double[][] yDotK,
                                                       final boolean forward) {
    final RungeKuttaStepInterpolator interpolator = (RungeKuttaStepInterpolator) prototype.copy();
    interpolator.reinitialize(this, yTmp, yDotK, forward,
                              equations.getPrimaryMapper(), equations.getSecondaryMappers());
    interpolator.storeTime(equations.getTime());
    return interpolator;
  }

  /**
   * Compute derivative for all intermediate stages.
   *
   * @param y state array
   * @param yDotK derivatives array
   * @param yTmp temporary state array
   * @param stages number of stages
   * @throws MaxCountExceededException if max evaluations is exceeded
   * @throws DimensionMismatchException if dimensions mismatch
   */
  private void computeStageDerivatives(final double[] y, final double[][] yDotK, final double[] yTmp, final int stages)
      throws MaxCountExceededException, DimensionMismatchException {
    for (int k = 1; k < stages; ++k) {
      for (int j = 0; j < y.length; ++j) {
        double sum = a[k-1][0] * yDotK[0][j];
        for (int l = 1; l < k; ++l) {
          sum += a[k-1][l] * yDotK[l][j];
        }
        yTmp[j] = y[j] + stepSize * sum;
      }
      computeDerivatives(stepStart + c[k-1] * stepSize, yTmp, yDotK[k]);
    }
  }

  /**
   * Estimate state at the end of the step.
   *
   * @param y state array
   * @param yDotK derivatives array
   * @param yTmp temporary state array
   * @param stages number of stages
   */
  private void estimateState(final double[] y, final double[][] yDotK, final double[] yTmp, final int stages) {
    for (int j = 0; j < y.length; ++j) {
      double sum = b[0] * yDotK[0][j];
      for (int l = 1; l < stages; ++l) {
        sum += b[l] * yDotK[l][j];
      }
      yTmp[j] = y[j] + stepSize * sum;
    }
  }

  /**
   * Handle step acceptance and discrete events.
   *
   * @param interpolator interpolator to use
   * @param yTmp temporary state array
   * @param y state array
   * @param yDotK derivatives array
   * @param yDotTmp temporary derivatives array
   * @param stages number of stages
   * @param t integration end time
   * @return next step start time
   * @throws MaxCountExceededException if max evaluations is exceeded
   * @throws DimensionMismatchException if dimensions mismatch
   * @throws NoBracketingException if bracketing fails
   */
  private double acceptAndResetStep(final RungeKuttaStepInterpolator interpolator,
                                    final double[] yTmp, final double[] y,
                                    final double[][] yDotK, final double[] yDotTmp,
                                    final int stages, final double t)
      throws MaxCountExceededException, DimensionMismatchException, NoBracketingException {
    interpolator.storeTime(stepStart + stepSize);
    System.arraycopy(yTmp, 0, y, 0, y.length);
    System.arraycopy(yDotK[stages - 1], 0, yDotTmp, 0, y.length);
    return acceptStep(interpolator, y, yDotTmp, t);
  }

  /**
   * Prepare next step and control step size.
   *
   * @param interpolator interpolator to use
   * @param forward forward integration flag
   * @param t integration end time
   */
  private void prepareNextStep(final RungeKuttaStepInterpolator interpolator, final boolean forward, final double t) {
    if (!isLastStep) {
      // prepare next step
      interpolator.storeTime(stepStart);

      // stepsize control for next step
      final double nextT = stepStart + stepSize;
      final boolean nextIsLast = forward ? (nextT >= t) : (nextT <= t);
      if (nextIsLast) {
        stepSize = t - stepStart;
      }
    }
  }

  /**
   * Reset integration state parameters.
   */
  private void resetIntegrationState() {
    stepStart = Double.NaN;
    stepSize  = Double.NaN;
  }

}
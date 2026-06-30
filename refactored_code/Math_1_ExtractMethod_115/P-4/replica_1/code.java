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
    final double[][] yDotK = createYDotK(stages, y0.length);
    final double[] yTmp    = y0.clone();
    final double[] yDotTmp = new double[y0.length];

    // set up an interpolator sharing the integrator arrays
    final RungeKuttaStepInterpolator interpolator = setUpInterpolator(equations, yTmp, yDotK, forward);

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
      estimateStateAtStepEnd(y, yDotK, yTmp, stages);

      // discrete events handling
      interpolator.storeTime(stepStart + stepSize);
      System.arraycopy(yTmp, 0, y, 0, y0.length);
      System.arraycopy(yDotK[stages - 1], 0, yDotTmp, 0, y0.length);
      stepStart = acceptStep(interpolator, y, yDotTmp, t);

      handleNextStepSize(interpolator, forward, t);

    } while (!isLastStep);

    // dispatch results
    finalizeIntegration(equations, y);

  }

  /** Create the external weights step derivative holder array.
   * @param stages number of stages
   * @param length state vector length
   * @return a new 2D array for storing derivatives
   */
  private double[][] createYDotK(final int stages, final int length) {
    final double[][] yDotK = new double[stages][];
    for (int i = 0; i < stages; ++i) {
      yDotK[i] = new double[length];
    }
    return yDotK;
  }

  /** Set up the step interpolator.
   * @param equations equations to integrate
   * @param yTmp temporary state array
   * @param yDotK derivatives at each stage
   * @param forward integration direction indicator
   * @return initialized interpolator
   */
  private RungeKuttaStepInterpolator setUpInterpolator(final ExpandableStatefulODE equations,
                                                        final double[] yTmp,
                                                        final double[][] yDotK,
                                                        final boolean forward) {
    final RungeKuttaStepInterpolator interpolator = (RungeKuttaStepInterpolator) prototype.copy();
    interpolator.reinitialize(this, yTmp, yDotK, forward,
                              equations.getPrimaryMapper(), equations.getSecondaryMappers());
    interpolator.storeTime(equations.getTime());
    return interpolator;
  }

  /** Compute the derivatives for all intermediate stages.
   * @param y current state
   * @param yDotK derivatives array to fill
   * @param yTmp temporary state array
   * @param stages number of stages
   */
  private void computeStageDerivatives(final double[] y, final double[][] yDotK,
                                       final double[] yTmp, final int stages)
      throws MaxCountExceededException, DimensionMismatchException {
    for (int k = 1; k < stages; ++k) {
      for (int j = 0; j < y.length; ++j) {
        double sum = a[k - 1][0] * yDotK[0][j];
        for (int l = 1; l < k; ++l) {
          sum += a[k - 1][l] * yDotK[l][j];
        }
        yTmp[j] = y[j] + stepSize * sum;
      }
      computeDerivatives(stepStart + c[k - 1] * stepSize, yTmp, yDotK[k]);
    }
  }

  /** Estimate the state at the end of the step.
   * @param y current state
   * @param yDotK derivatives array
   * @param yTmp temporary state array to store result
   * @param stages number of stages
   */
  private void estimateStateAtStepEnd(final double[] y, final double[][] yDotK,
                                      final double[] yTmp, final int stages) {
    for (int j = 0; j < y.length; ++j) {
      double sum = b[0] * yDotK[0][j];
      for (int l = 1; l < stages; ++l) {
        sum += b[l] * yDotK[l][j];
      }
      yTmp[j] = y[j] + stepSize * sum;
    }
  }

  /** Handle the step size control for the next step.
   * @param interpolator step interpolator
   * @param forward integration direction indicator
   * @param t target integration end time
   */
  private void handleNextStepSize(final RungeKuttaStepInterpolator interpolator,
                                  final boolean forward, final double t) {
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

  /** Finalize the integration by updating the state and cleaning up variables.
   * @param equations equations to update
   * @param y final state array
   */
  private void finalizeIntegration(final ExpandableStatefulODE equations, final double[] y) {
    equations.setTime(stepStart);
    equations.setCompleteState(y);

    stepStart = Double.NaN;
    stepSize  = Double.NaN;
  }

}

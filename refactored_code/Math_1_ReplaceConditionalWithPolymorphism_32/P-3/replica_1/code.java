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
package org.apache.commons.math3.analysis.solvers;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

/**
 * This class implements the <a href="http://mathworld.wolfram.com/MullersMethod.html">
 * Muller's Method</a> for root finding of real univariate functions. For
 * reference, see <b>Elementary Numerical Analysis</b>, ISBN 0070124477,
 * chapter 3.
 * <p>
 * Muller's method applies to both real and complex functions, but here we
 * restrict ourselves to real functions.
 * This class differs from {@link MullerSolver} in the way it avoids complex
 * operations.</p>
 * Muller's original method would have function evaluation at complex point.
 * Since our f(x) is real, we have to find ways to avoid that. Bracketing
 * condition is one way to go: by requiring bracketing in every iteration,
 * the newly computed approximation is guaranteed to be real.</p>
 * <p>
 * Normally Muller's method converges quadratically in the vicinity of a
 * zero, however it may be very slow in regions far away from zeros. For
 * example, f(x) = exp(x) - 1, min = -50, max = 100. In such case we use
 * bisection as a safety backup if it performs very poorly.</p>
 * <p>
 * The formulas here use divided differences directly.</p>
 * 
 * @version $Id$
 * @since 1.2
 * @see MullerSolver2
 */
public class MullerSolver extends AbstractUnivariateSolver {

    /** Default absolute accuracy. */
    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    private final StepStrategy mullerStep = new MullerStep();
    private final StepStrategy bisectionStep = new BisectionStep();

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public MullerSolver() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }
    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }
    /**
     * Construct a solver.
     *
     * @param relativeAccuracy Relative accuracy.
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver(double relativeAccuracy,
                        double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double doSolve()
        throws TooManyEvaluationsException,
               NumberIsTooLargeException,
               NoBracketingException {
        final double min = getMin();
        final double max = getMax();
        final double initial = getStartValue();

        final double functionValueAccuracy = getFunctionValueAccuracy();

        verifySequence(min, initial, max);

        // check for zeros before verifying bracketing
        final double fMin = computeObjectiveValue(min);
        if (FastMath.abs(fMin) < functionValueAccuracy) {
            return min;
        }
        final double fMax = computeObjectiveValue(max);
        if (FastMath.abs(fMax) < functionValueAccuracy) {
            return max;
        }
        final double fInitial = computeObjectiveValue(initial);
        if (FastMath.abs(fInitial) <  functionValueAccuracy) {
            return initial;
        }

        verifyBracketing(min, max);

        if (isBracketing(min, initial)) {
            return solve(min, initial, fMin, fInitial);
        } else {
            return solve(initial, max, fInitial, fMax);
        }
    }

    /**
     * Find a real root in the given interval.
     *
     * @param min Lower bound for the interval.
     * @param max Upper bound for the interval.
     * @param fMin function value at the lower bound.
     * @param fMax function value at the upper bound.
     * @return the point at which the function value is zero.
     * @throws TooManyEvaluationsException if the allowed number of calls to
     * the function to be solved has been exhausted.
     */
    private double solve(double min, double max,
                         double fMin, double fMax)
        throws TooManyEvaluationsException {
        final double relativeAccuracy = getRelativeAccuracy();
        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();

        // [x0, x2] is the bracketing interval in each iteration
        // x1 is the last approximation and an interpolation point in (x0, x2)
        // x is the new root approximation and new x1 for next round
        // d01, d12, d012 are divided differences

        double x1 = 0.5 * (min + max);
        double y1 = computeObjectiveValue(x1);

        SolverState state = new SolverState(min, fMin, x1, y1, max, fMax, Double.POSITIVE_INFINITY);

        while (true) {
            // Muller's method employs quadratic interpolation through
            // x0, x1, x2 and x is the zero of the interpolating parabola.
            // Due to bracketing condition, this parabola must have two
            // real roots and we choose one in [x0, x2] to be x.
            final double d01 = (state.y1 - state.y0) / (state.x1 - state.x0);
            final double d12 = (state.y2 - state.y1) / (state.x2 - state.x1);
            final double d012 = (d12 - d01) / (state.x2 - state.x0);
            final double c1 = d01 + (state.x1 - state.x0) * d012;
            final double delta = c1 * c1 - 4 * state.y1 * d012;
            final double xplus = state.x1 + (-2.0 * state.y1) / (c1 + FastMath.sqrt(delta));
            final double xminus = state.x1 + (-2.0 * state.y1) / (c1 - FastMath.sqrt(delta));
            // xplus and xminus are two roots of parabola and at least
            // one of them should lie in (x0, x2)
            final double x = isSequence(state.x0, xplus, state.x2) ? xplus : xminus;
            final double y = computeObjectiveValue(x);

            // check for convergence
            final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
            if (FastMath.abs(x - state.oldx) <= tolerance ||
                FastMath.abs(y) <= functionValueAccuracy) {
                return x;
            }

            StepStrategy strategy = selectStrategy(x, state);
            state = strategy.update(state, x, y);
        }
    }

    private StepStrategy selectStrategy(double x, SolverState state) {
        boolean bisect = (x < state.x1 && (state.x1 - state.x0) > 0.95 * (state.x2 - state.x0)) ||
                         (x > state.x1 && (state.x2 - state.x1) > 0.95 * (state.x2 - state.x0)) ||
                         (x == state.x1);
        return bisect ? bisectionStep : mullerStep;
    }

    private static class SolverState {
        final double x0;
        final double y0;
        final double x1;
        final double y1;
        final double x2;
        final double y2;
        final double oldx;

        SolverState(double x0, double y0, double x1, double y1, double x2, double y2, double oldx) {
            this.x0 = x0;
            this.y0 = y0;
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
            this.oldx = oldx;
        }
    }

    private interface StepStrategy {
        SolverState update(SolverState state, double x, double y);
    }

    private class MullerStep implements StepStrategy {
        @Override
        public SolverState update(SolverState state, double x, double y) {
            double nextX0 = x < state.x1 ? state.x0 : state.x1;
            double nextY0 = x < state.x1 ? state.y0 : state.y1;
            double nextX2 = x > state.x1 ? state.x2 : state.x1;
            double nextY2 = x > state.x1 ? state.y2 : state.y1;
            return new SolverState(nextX0, nextY0, x, y, nextX2, nextY2, x);
        }
    }

    private class BisectionStep implements StepStrategy {
        @Override
        public SolverState update(SolverState state, double x, double y) {
            double xm = 0.5 * (state.x0 + state.x2);
            double ym = computeObjectiveValue(xm);
            double nextX0, nextY0, nextX2, nextY2;
            if (FastMath.signum(state.y0) + FastMath.signum(ym) == 0.0) {
                nextX0 = state.x0;
                nextY0 = state.y0;
                nextX2 = xm;
                nextY2 = ym;
            } else {
                nextX0 = xm;
                nextY0 = ym;
                nextX2 = state.x2;
                nextY2 = state.y2;
            }
            double nextX1 = 0.5 * (nextX0 + nextX2);
            double nextY1 = computeObjectiveValue(nextX1);
            return new SolverState(nextX0, nextY0, nextX1, nextY1, nextX2, nextY2, Double.POSITIVE_INFINITY);
        }
    }
}
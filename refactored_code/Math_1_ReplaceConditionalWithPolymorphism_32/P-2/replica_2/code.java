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

        IntervalSelector selector = isBracketing(min, initial) ? 
                new LowerIntervalSelector() : new UpperIntervalSelector();
        return selector.solve(this, min, initial, max, fMin, fInitial, fMax);
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

        double x0 = min;
        double y0 = fMin;
        double x2 = max;
        double y2 = fMax;
        double x1 = 0.5 * (x0 + x2);
        double y1 = computeObjectiveValue(x1);

        double oldx = Double.POSITIVE_INFINITY;
        while (true) {
            // Muller's method employs quadratic interpolation through
            // x0, x1, x2 and x is the zero of the interpolating parabola.
            // Due to bracketing condition, this parabola must have two
            // real roots and we choose one in [x0, x2] to be x.
            final double d01 = (y1 - y0) / (x1 - x0);
            final double d12 = (y2 - y1) / (x2 - x1);
            final double d012 = (d12 - d01) / (x2 - x0);
            final double c1 = d01 + (x1 - x0) * d012;
            final double delta = c1 * c1 - 4 * y1 * d012;
            final double xplus = x1 + (-2.0 * y1) / (c1 + FastMath.sqrt(delta));
            final double xminus = x1 + (-2.0 * y1) / (c1 - FastMath.sqrt(delta));
            // xplus and xminus are two roots of parabola and at least
            // one of them should lie in (x0, x2)
            final double x = isSequence(x0, xplus, x2) ? xplus : xminus;
            final double y = computeObjectiveValue(x);

            // check for convergence
            final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
            if (FastMath.abs(x - oldx) <= tolerance ||
                FastMath.abs(y) <= functionValueAccuracy) {
                return x;
            }

            // Bisect if convergence is too slow. Bisection would waste
            // our calculation of x, hopefully it won't happen often.
            // the real number equality test x == x1 is intentional and
            // completes the proximity tests above it
            boolean bisect = (x < x1 && (x1 - x0) > 0.95 * (x2 - x0)) ||
                             (x > x1 && (x2 - x1) > 0.95 * (x2 - x0)) ||
                             (x == x1);

            IterationStep step = bisect ? new BisectionIterationStep() : new MullerIterationStep(x, y);
            SolverState state = new SolverState(x0, y0, x1, y1, x2, y2, oldx);
            state = step.update(state, this);
            x0 = state.x0;
            y0 = state.y0;
            x1 = state.x1;
            y1 = state.y1;
            x2 = state.x2;
            y2 = state.y2;
            oldx = state.oldx;
        }
    }

    private interface IntervalSelector {
        double solve(MullerSolver solver, double min, double initial, double max, double fMin, double fInitial, double fMax);
    }

    private static class LowerIntervalSelector implements IntervalSelector {
        @Override
        public double solve(MullerSolver solver, double min, double initial, double max, double fMin, double fInitial, double fMax) {
            return solver.solve(min, initial, fMin, fInitial);
        }
    }

    private static class UpperIntervalSelector implements IntervalSelector {
        @Override
        public double solve(MullerSolver solver, double min, double initial, double max, double fMin, double fInitial, double fMax) {
            return solver.solve(initial, max, fInitial, fMax);
        }
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

    private interface IterationStep {
        SolverState update(SolverState state, MullerSolver solver);
    }

    private static class MullerIterationStep implements IterationStep {
        private final double x;
        private final double y;

        MullerIterationStep(double x, double y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public SolverState update(SolverState state, MullerSolver solver) {
            double x0 = x < state.x1 ? state.x0 : state.x1;
            double y0 = x < state.x1 ? state.y0 : state.y1;
            double x2 = x > state.x1 ? state.x2 : state.x1;
            double y2 = x > state.x1 ? state.y2 : state.y1;
            double x1 = x;
            double y1 = y;
            double oldx = x;
            return new SolverState(x0, y0, x1, y1, x2, y2, oldx);
        }
    }

    private static class BisectionIterationStep implements IterationStep {
        @Override
        public SolverState update(SolverState state, MullerSolver solver) {
            double xm = 0.5 * (state.x0 + state.x2);
            double ym = solver.computeObjectiveValue(xm);
            double x0, y0, x2, y2;
            if (FastMath.signum(state.y0) + FastMath.signum(ym) == 0.0) {
                x0 = state.x0;
                y0 = state.y0;
                x2 = xm;
                y2 = ym;
            } else {
                x0 = xm;
                y0 = ym;
                x2 = state.x2;
                y2 = state.y2;
            }
            double x1 = 0.5 * (x0 + x2);
            double y1 = solver.computeObjectiveValue(x1);
            double oldx = Double.POSITIVE_INFINITY;
            return new SolverState(x0, y0, x1, y1, x2, y2, oldx);
        } 
    }
}

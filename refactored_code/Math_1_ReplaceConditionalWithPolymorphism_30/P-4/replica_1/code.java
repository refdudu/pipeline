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

import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

/**
 * This class implements the <a href="http://mathworld.wolfram.com/BrentsMethod.html">
 * Brent algorithm</a> for finding zeros of real univariate functions.
 * The function should be continuous but not necessarily smooth.
 * The {@code solve} method returns a zero {@code x} of the function {@code f}
 * in the given interval {@code [a, b]} to within a tolerance
 * {@code 6 eps abs(x) + t} where {@code eps} is the relative accuracy and
 * {@code t} is the absolute accuracy.
 * The given interval must bracket the root.
 *
 * @version $Id$
 */
public class BrentSolver extends AbstractUnivariateSolver {

    /** Default absolute accuracy. */
    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    private static final SolverStep BISECTION_STEP = new BisectionStep();
    private static final SolverStep INTERPOLATION_STEP = new InterpolationStep();

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public BrentSolver() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }

    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    public BrentSolver(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }

    /**
     * Construct a solver.
     *
     * @param relativeAccuracy Relative accuracy.
     * @param absoluteAccuracy Absolute accuracy.
     */
    public BrentSolver(double relativeAccuracy,
                       double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
    }

    /**
     * Construct a solver.
     *
     * @param relativeAccuracy Relative accuracy.
     * @param absoluteAccuracy Absolute accuracy.
     * @param functionValueAccuracy Function value accuracy.
     */
    public BrentSolver(double relativeAccuracy,
                       double absoluteAccuracy,
                       double functionValueAccuracy) {
        super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double doSolve()
        throws NoBracketingException,
               TooManyEvaluationsException,
               NumberIsTooLargeException {
        double min = getMin();
        double max = getMax();
        final double initial = getStartValue();
        final double functionValueAccuracy = getFunctionValueAccuracy();

        verifySequence(min, initial, max);

        // Return the initial guess if it is good enough.
        double yInitial = computeObjectiveValue(initial);
        if (FastMath.abs(yInitial) <= functionValueAccuracy) {
            return initial;
        }

        // Return the first endpoint if it is good enough.
        double yMin = computeObjectiveValue(min);
        if (FastMath.abs(yMin) <= functionValueAccuracy) { 
            return min;
        }

        // Reduce interval if min and initial bracket the root.
        if (yInitial * yMin < 0) {
            return brent(min, initial, yMin, yInitial);
        }

        // Return the second endpoint if it is good enough.
        double yMax = computeObjectiveValue(max);
        if (FastMath.abs(yMax) <= functionValueAccuracy) {
            return max;
        }

        // Reduce interval if initial and max bracket the root.
        if (yInitial * yMax < 0) {
            return brent(initial, max, yInitial, yMax);
        }

        throw new NoBracketingException(min, max, yMin, yMax);
    }

    /**
     * Search for a zero inside the provided interval.
     * This implementation is based on the algorithm described at page 58 of
     * the book
     * <quote>
     *  <b>Algorithms for Minimization Without Derivatives</b>
     *  <it>Richard P. Brent</it>
     *  Dover 0-486-41998-3
     * </quote>
     *
     * @param lo Lower bound of the search interval.
     * @param hi Higher bound of the search interval.
     * @param fLo Function value at the lower bound of the search interval.
     * @param fHi Function value at the higher bound of the search interval.
     * @return the value where the function is zero.
     */
    private double brent(double lo, double hi,
                         double fLo, double fHi) {
        double a = lo;
        double fa = fLo;
        double b = hi;
        double fb = fHi;
        double c = a;
        double fc = fa;
        double d = b - a;
        double e = d;

        final double t = getAbsoluteAccuracy();
        final double eps = getRelativeAccuracy();

        while (true) {
            if (FastMath.abs(fc) < FastMath.abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            final double tol = 2 * eps * FastMath.abs(b) + t;
            final double m = 0.5 * (c - b);

            if (FastMath.abs(m) <= tol ||
                Precision.equals(fb, 0))  {
                return b;
            }

            final SolverStep stepStrategy = (FastMath.abs(e) < tol || FastMath.abs(fa) <= FastMath.abs(fb))
                    ? BISECTION_STEP
                    : INTERPOLATION_STEP;

            final StepResult stepResult = stepStrategy.apply(a, b, c, fa, fb, fc, d, e, m, tol);
            d = stepResult.d;
            e = stepResult.e;

            a = b;
            fa = fb;

            if (FastMath.abs(d) > tol) {
                b += d;
            } else if (m > 0) {
                b += tol;
            } else {
                b -= tol;
            }
            fb = computeObjectiveValue(b);
            if ((fb > 0 && fc > 0) ||
                (fb <= 0 && fc <= 0)) {
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }
        }
    }

    private static class StepResult {
        final double d;
        final double e;

        StepResult(double d, double e) {
            this.d = d;
            this.e = e;
        }
    }

    private interface SolverStep {
        StepResult apply(double a, double b, double c, double fa, double fb, double fc, double d, double e, double m, double tol);
    }

    private static class BisectionStep implements SolverStep {
        @Override
        public StepResult apply(double a, double b, double c, double fa, double fb, double fc, double d, double e, double m, double tol) {
            return new StepResult(m, m);
        }
    }

    private static class InterpolationStep implements SolverStep {
        private final InterpolationStrategy linearInterpolation = new LinearInterpolation();
        private final InterpolationStrategy inverseQuadraticInterpolation = new InverseQuadraticInterpolation();

        @Override
        public StepResult apply(double a, double b, double c, double fa, double fb, double fc, double d, double e, double m, double tol) {
            final double s = fb / fa;
            final InterpolationStrategy strategy = (a == c) ? linearInterpolation : inverseQuadraticInterpolation;
            final InterpolationResult result = strategy.compute(a, b, c, fa, fb, fc, m, s);
            
            double p = result.p;
            double q = result.q;

            if (p > 0) {
                q = -q;
            } else {
                p = -p;
            }
            
            final double previousE = e;
            final double nextE = d;
            
            if (p >= 1.5 * m * q - FastMath.abs(tol * q) ||
                p >= FastMath.abs(0.5 * previousE * q)) {
                return new StepResult(m, m);
            } else {
                return new StepResult(p / q, nextE);
            }
        }
    }

    private static class InterpolationResult {
        final double p;
        final double q;

        InterpolationResult(double p, double q) {
            this.p = p;
            this.q = q;
        }
    }

    private interface InterpolationStrategy {
        InterpolationResult compute(double a, double b, double c, double fa, double fb, double fc, double m, double s);
    }

    private static class LinearInterpolation implements InterpolationStrategy {
        @Override
        public InterpolationResult compute(double a, double b, double c, double fa, double fb, double fc, double m, double s) {
            double p = 2 * m * s;
            double q = 1 - s;
            return new InterpolationResult(p, q);
        }
    }

    private static class InverseQuadraticInterpolation implements InterpolationStrategy {
        @Override
        public InterpolationResult compute(double a, double b, double c, double fa, double fb, double fc, double m, double s) {
            double q = fa / fc;
            final double r = fb / fc;
            double p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
            q = (q - 1) * (r - 1) * (s - 1);
            return new InterpolationResult(p, q);
        }
    }
}
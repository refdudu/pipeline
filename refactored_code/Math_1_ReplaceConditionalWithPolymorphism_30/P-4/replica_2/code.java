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
            if (FastMath.abs(e) < tol ||
                FastMath.abs(fa) <= FastMath.abs(fb)) {
                // Force bisection.
                d = m;
                e = d;
            } else {
                double s = fb / fa;
                double p;
                double q;

                // Replace conditional logic with Polymorphism for interpolation method selection
                Interpolator interpolator = (a == c) ?
                        new LinearInterpolator() :
                        new InverseQuadraticInterpolator();

                InterpolationStep step = interpolator.calculateStep(m, s, a, b, fa, fb, fc);
                p = step.getP();
                q = step.getQ();

                // Replace sign correction conditional structure with polymorphism
                SignNormalizer normalizer = (p > 0) ?
                        new PositivePNormalizer() :
                        new NegativePNormalizer();

                NormalizedValues normalized = normalizer.normalize(p, q);
                p = normalized.getP();
                q = normalized.getQ();

                s = e;
                e = d;
                if (p >= 1.5 * m * q - FastMath.abs(tol * q) ||
                    p >= FastMath.abs(0.5 * s * q)) {
                    // Inverse quadratic interpolation gives a value
                    // in the wrong direction, or progress is slow.
                    // Fall back to bisection.
                    d = m;
                    e = d;
                } else {
                    d = p / q;
                }
            }
            a = b;
            fa = fb;

            // Replace step value selection conditionals with polymorphism
            StepIncrementer incrementer;
            if (FastMath.abs(d) > tol) {
                incrementer = new DirectIncrementer();
            } else if (m > 0) {
                incrementer = new PositiveTolIncrementer();
            } else {
                incrementer = new NegativeTolIncrementer();
            }
            b += incrementer.getIncrement(d, tol);

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

    // --- Helper classes supporting Polymorphism ---

    private interface Interpolator {
        InterpolationStep calculateStep(double m, double s, double a, double b, double fa, double fb, double fc);
    }

    private static final class InterpolationStep {
        private final double p;
        private final double q;

        public InterpolationStep(double p, double q) {
            this.p = p;
            this.q = q;
        }

        public double getP() {
            return p;
        }

        public double getQ() {
            return q;
        }
    }

    private static final class LinearInterpolator implements Interpolator {
        @Override
        public InterpolationStep calculateStep(double m, double s, double a, double b, double fa, double fb, double fc) {
            return new InterpolationStep(2 * m * s, 1 - s);
        } 
    }

    private static final class InverseQuadraticInterpolator implements Interpolator {
        @Override
        public InterpolationStep calculateStep(double m, double s, double a, double b, double fa, double fb, double fc) {
            double q = fa / fc;
            double r = fb / fc;
            double p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
            double qVal = (q - 1) * (r - 1) * (s - 1);
            return new InterpolationStep(p, qVal);
        }
    }

    private interface SignNormalizer {
        NormalizedValues normalize(double p, double q);
    }

    private static final class NormalizedValues {
        private final double p;
        private final double q;

        public NormalizedValues(double p, double q) {
            this.p = p;
            this.q = q;
        }

        public double getP() {
            return p;
        }

        public double getQ() {
            return q;
        }
    }

    private static final class PositivePNormalizer implements SignNormalizer {
        @Override
        public NormalizedValues normalize(double p, double q) { 
            return new NormalizedValues(p, -q);
        }
    }

    private static final class NegativePNormalizer implements SignNormalizer {
        @Override
        public NormalizedValues normalize(double p, double q) {
            return new NormalizedValues(-p, q);
        }
    }

    private interface StepIncrementer {
        double getIncrement(double d, double tol);
    }

    private static final class DirectIncrementer implements StepIncrementer {
        @Override
        public double getIncrement(double d, double tol) {
            return d;
        }
    }

    private static final class PositiveTolIncrementer implements StepIncrementer {
        @Override
        public double getIncrement(double d, double tol) {
            return tol;
        }
    }

    private static final class NegativeTolIncrementer implements StepIncrementer {
        @Override
        public double getIncrement(double d, double tol) {
            return -tol;
        }
    }
}

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
 * Ah, welcome to my humble bakery! I may not know what a "Java class" is,
 * but I know a messy recipe when I see one. When you have different ways
 * of kneading your dough (what you call "interpolation"), you shouldn't clutter
 * the main worktable. Instead, you delegate each recipe to a specialized apprentice
 * who knows exactly how to handle that specific loaf! We bakers call this division of labor,
 * but my grandson tells me the proper term in your language is "Polymorphism."
 *
 * Here is the Brent recipe, now beautiful, clean, and perfectly proofed.
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
                
                // Instead of a messy "if-else" conditional right in the middle of our baking,
                // we let the Head Baker (our Factory) decide which Apprentice (Strategy) should handle this batch!
                final InterpolationStrategy strategy = InterpolationStrategyFactory.getStrategy(a, c);
                final InterpolationValues values = strategy.calculate(a, b, c, fa, fb, fc, m, s);
                
                double p = values.getP();
                double q = values.getQ();

                if (p > 0) {
                    q = -q;
                } else {
                    p = -p;
                }
                s = e;
                e = d;
                if (p >= 1.5 * m * q - FastMath.abs(tol * q) ||
                    p >= FastMath.abs(0.5 * s * q)) {
                    // Fall back to bisection.
                    d = m;
                    e = d;
                } else {
                    d = p / q;
                }
            }
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

    // =========================================================================
    // THE ARTISAN BAKER'S REFACTORED APPRENTICES (POLYMORPHIC STRATEGIES)
    // =========================================================================

    /**
     * Represents the calculated flour portions (p and q) from our kneading strategy.
     */
    private static class InterpolationValues {
        private final double p;
        private final double q;

        public InterpolationValues(double p, double q) {
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

    /**
     * The general recipe for any interpolation apprentice to follow.
     */
    private interface InterpolationStrategy {
        InterpolationValues calculate(double a, double b, double c,
                                      double fa, double fb, double fc,
                                      double m, double s);
    }

    /**
     * Apprentice responsible for Linear Interpolation.
     * Simple, traditional, like baking a classic baguette.
     */
    private static class LinearInterpolation implements InterpolationStrategy {
        @Override
        public InterpolationValues calculate(double a, double b, double c,
                                              double fa, double fb, double fc,
                                              double m, double s) {
            double p = 2 * m * s;
            double q = 1 - s;
            return new InterpolationValues(p, q);
        }
    }

    /**
     * Apprentice responsible for Inverse Quadratic Interpolation.
     * Highly complex, like shaping an intricate braided Challah bread.
     */
    private static class InverseQuadraticInterpolation implements InterpolationStrategy {
        @Override
        public InterpolationValues calculate(double a, double b, double c,
                                              double fa, double fb, double fc,
                                              double m, double s) {
            double q = fa / fc;
            final double r = fb / fc;
            double p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
            q = (q - 1) * (r - 1) * (s - 1);
            return new InterpolationValues(p, q);
        }
    }

    /**
     * The Head Baker who assigns the batch to the correct apprentice based on flour state.
     */
    private static class InterpolationStrategyFactory {
        public static InterpolationStrategy getStrategy(double a, double c) {
            if (a == c) {
                return new LinearInterpolation();
            } else {
                return new InverseQuadraticInterpolation();
            }
        }
    }
}
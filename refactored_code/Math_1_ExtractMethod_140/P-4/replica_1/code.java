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
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;

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
 * Except for the initial [min, max], it does not require bracketing
 * condition, e.g. f(x0), f(x1), f(x2) can have the same sign. If complex
 * number arises in the computation, we simply use its modulus as real
 * approximation.</p>
 * <p>
 * Because the interval may not be bracketing, bisection alternative is
 * not applicable here. However in practice our treatment usually works
 * well, especially near real zeroes where the imaginary part of complex
 * approximation is often negligible.</p>
 * <p>
 * The formulas here do not use divided differences directly.</p>
 * 
 * @version $Id$
 * @since 1.2
 * @see MullerSolver
 */
public class MullerSolver2 extends AbstractUnivariateSolver {

    /** Default absolute accuracy. */
    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public MullerSolver2() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }

    /**
     * Construct a solver.
     * 
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver2(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }

    /**
     * Construct a solver.
     * 
     * @param relativeAccuracy Relative accuracy.
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver2(double relativeAccuracy,
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

        verifyInterval(min, max);

        final double relativeAccuracy = getRelativeAccuracy();
        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();

        double x0 = min;
        double y0 = computeObjectiveValue(x0);
        if (isCloseToZero(y0, functionValueAccuracy)) {
            return x0;
        }
        double x1 = max;
        double y1 = computeObjectiveValue(x1);
        if (isCloseToZero(y1, functionValueAccuracy)) {
            return x1;
        }

        checkBracketing(x0, x1, y0, y1);

        double x2 = 0.5 * (x0 + x1);
        double y2 = computeObjectiveValue(x2);

        double oldx = Double.POSITIVE_INFINITY;
        while (true) {
            final double q = (x2 - x1) / (x1 - x0);
            final double a = q * (y2 - (1 + q) * y1 + q * y0);
            final double b = (2 * q + 1) * y2 - (1 + q) * (1 + q) * y1 + q * q * y0;
            final double c = (1 + q) * y2;
            final double delta = b * b - 4 * a * c;
            
            double x;
            final double denominator = computeDenominator(b, delta);
            if (denominator != 0) {
                x = computeNextX(x2, x1, c, denominator, absoluteAccuracy);
            } else {
                x = getRandomX(min, max);
                oldx = Double.POSITIVE_INFINITY;
            }
            final double y = computeObjectiveValue(x);

            if (isConverged(x, oldx, y, relativeAccuracy, absoluteAccuracy, functionValueAccuracy)) {
                return x;
            }

            x0 = x1;
            y0 = y1;
            x1 = x2;
            y1 = y2;
            x2 = x;
            y2 = y;
            oldx = x;
        }
    }

    /**
     * Checks if the function value is close to zero within the given accuracy.
     */
    private boolean isCloseToZero(double value, double threshold) {
        return FastMath.abs(value) < threshold;
    }

    /**
     * Verifies that the initial endpoints bracket a root.
     */
    private void checkBracketing(double x0, double x1, double y0, double y1) {
        if (y0 * y1 > 0) {
            throw new NoBracketingException(x0, x1, y0, y1);
        }
    }

    /**
     * Computes the denominator for Muller's method step.
     */
    private double computeDenominator(double b, double delta) {
        if (delta >= 0.0) {
            final double dplus = b + FastMath.sqrt(delta);
            final double dminus = b - FastMath.sqrt(delta);
            return FastMath.abs(dplus) > FastMath.abs(dminus) ? dplus : dminus;
        } else {
            return FastMath.sqrt(b * b - delta);
        }
    }

    /**
     * Computes the next root approximation.
     */
    private double computeNextX(double x2, double x1, double c, double denominator, double absoluteAccuracy) {
        double x = x2 - 2.0 * c * (x2 - x1) / denominator;
        return perturbIfCoincides(x, x1, x2, absoluteAccuracy);
    }

    /**
     * Perturbs the approximation if it exactly coincides with x1 or x2.
     */
    private double perturbIfCoincides(double x, double x1, double x2, double absoluteAccuracy) {
        double result = x;
        while (result == x1 || result == x2) {
            result += absoluteAccuracy;
        }
        return result;
    }

    /**
     * Generates a random value in the interval [min, max].
     */
    private double getRandomX(double min, double max) {
        return min + FastMath.random() * (max - min);
    }

    /**
     * Checks if the algorithm has converged on a root.
     */
    private boolean isConverged(double x, double oldx, double y,
                                double relativeAccuracy, double absoluteAccuracy,
                                double functionValueAccuracy) {
        final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
        return FastMath.abs(x - oldx) <= tolerance ||
               FastMath.abs(y) <= functionValueAccuracy;
    }
}
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
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

/**
 * Implements the <a href="http://mathworld.wolfram.com/RiddersMethod.html">
 * Ridders' Method</a> for root finding of real univariate functions. For
 * reference, see C. Ridders, <i>A new algorithm for computing a single root
 * of a real continuous function </i>, IEEE Transactions on Circuits and
 * Systems, 26 (1979), 979 - 980.
 * <p>
 * The function should be continuous but not necessarily smooth.</p>
 *
 * @version $Id$
 * @since 1.2
 */
public class RiddersSolver extends AbstractUnivariateSolver {
    /** Default absolute accuracy. */
    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public RiddersSolver() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }
    /**
     * Construct a solver.
     * 
     * @param absoluteAccuracy Absolute accuracy.
     */
    public RiddersSolver(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }
    /**
     * Construct a solver.
     * 
     * @param relativeAccuracy Relative accuracy.
     * @param absoluteAccuracy Absolute accuracy.
     */
    public RiddersSolver(double relativeAccuracy,
                         double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double doSolve()
        throws TooManyEvaluationsException,
               NoBracketingException {
        double min = getMin();
        double max = getMax();
        
        double x1 = min;
        double y1 = computeObjectiveValue(x1);
        double x2 = max;
        double y2 = computeObjectiveValue(x2);

        // check for zeros before verifying bracketing
        if (y1 == 0) {
            return min;
        }
        if (y2 == 0) {
            return max;
        }
        verifyBracketing(min, max);

        return solveInterval(new Interval(x1, y1, x2, y2));
    }

    /**
     * Solves the function within the given bracketing interval.
     *
     * @param interval the current bracketing interval.
     * @return the root approximation.
     */
    private double solveInterval(Interval interval) {
        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();
        final double relativeAccuracy = getRelativeAccuracy();

        double oldx = Double.POSITIVE_INFINITY;
        while (true) {
            final double x3 = 0.5 * (interval.x1 + interval.x2);
            final double y3 = computeObjectiveValue(x3);
            if (FastMath.abs(y3) <= functionValueAccuracy) {
                return x3;
            }

            final double correction = calculateCorrection(interval.x1, x3, interval.y1, interval.y2, y3);
            final double x = x3 - correction;
            final double y = computeObjectiveValue(x);

            if (hasConverged(x, oldx, y, absoluteAccuracy, relativeAccuracy, functionValueAccuracy)) {
                return x;
            }

            updateInterval(interval, x, y, x3, y3, correction);
            oldx = x;
        }
    }

    /**
     * Calculates the correction term using Ridders' formula.
     */
    private double calculateCorrection(double x1, double x3, double y1, double y2, double y3) {
        final double delta = 1 - (y1 * y2) / (y3 * y3);
        return (FastMath.signum(y2) * FastMath.signum(y3)) * (x3 - x1) / FastMath.sqrt(delta);
    }

    /**
     * Determines whether the root search has converged.
     */
    private boolean hasConverged(double x, double oldx, double y,
                                 double absoluteAccuracy, double relativeAccuracy,
                                 double functionValueAccuracy) {
        final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
        return FastMath.abs(x - oldx) <= tolerance || FastMath.abs(y) <= functionValueAccuracy;
    }

    /**
     * Updates the bracketing interval based on the Ridders' method iteration step.
     */
    private void updateInterval(Interval interval, double x, double y, double x3, double y3, double correction) {
        if (correction > 0.0) {             // x1 < x < x3
            if (FastMath.signum(interval.y1) + FastMath.signum(y) == 0.0) {
                interval.x2 = x;
                interval.y2 = y;
            } else {
                interval.x1 = x;
                interval.x2 = x3;
                interval.y1 = y;
                interval.y2 = y3;
            }
        } else {                            // x3 < x < x2
            if (FastMath.signum(interval.y2) + FastMath.signum(y) == 0.0) {
                interval.x1 = x;
                interval.y1 = y;
            } else {
                interval.x1 = x3;
                interval.x2 = x;
                interval.y1 = y3;
                interval.y2 = y;
            }
        }
    }

    /**
     * Helper class to represent a bracketing interval.
     */
    private static class Interval {
        private double x1;
        private double y1;
        private double x2;
        private double y2;

        Interval(double x1, double y1, double x2, double y2) {
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
        }
    }
}
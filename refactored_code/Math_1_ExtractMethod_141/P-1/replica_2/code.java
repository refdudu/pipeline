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
        final double min = getMin();
        final double max = getMax();
        
        final double y1 = computeObjectiveValue(min);
        final double y2 = computeObjectiveValue(max);

        if (y1 == 0) {
            return min;
        }
        if (y2 == 0) {
            return max;
        }
        verifyBracketing(min, max);

        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();
        final double relativeAccuracy = getRelativeAccuracy();

        final Interval interval = new Interval(min, y1, max, y2);
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

            if (hasConverged(x, oldx, y, relativeAccuracy, absoluteAccuracy, functionValueAccuracy)) {
                return x;
            }

            updateInterval(interval, x, y, x3, y3, correction);
            oldx = x;
        }
    }

    /**
     * Calculates the Ridders correction term.
     */
    private double calculateCorrection(double x1, double x3, double y1, double y2, double y3) {
        final double delta = 1.0 - (y1 * y2) / (y3 * y3);
        return (FastMath.signum(y2) * FastMath.signum(y3)) * (x3 - x1) / FastMath.sqrt(delta);
    }

    /**
     * Checks if the solver has converged under the given precision standards.
     */
    private boolean hasConverged(double x, double oldx, double y, 
                                 double relativeAccuracy, double absoluteAccuracy, double functionValueAccuracy) {
        final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
        return FastMath.abs(x - oldx) <= tolerance || FastMath.abs(y) <= functionValueAccuracy;
    }

    /**
     * Updates the bracketing interval state based on the calculated correction.
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
     * Helper class representing the current active bracketing interval.
     * Like a proofing basket holding our dough's parameters.
     */
    private static class Interval {
        double x1;
        double y1;
        double x2;
        double y2;

        Interval(double x1, double y1, double x2, double y2) {
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
        }
    }
}
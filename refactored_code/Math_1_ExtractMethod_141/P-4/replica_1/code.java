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

        final double x1 = min;
        final double y1 = computeObjectiveValue(x1);
        if (y1 == 0.0) {
            return min;
        }

        final double x2 = max;
        final double y2 = computeObjectiveValue(x2);
        if (y2 == 0.0) {
            return max;
        }

        verifyBracketing(min, max);

        return solveInterval(x1, y1, x2, y2);
    }

    /**
     * Solves the root finding problem using Ridders' method on the bracketed interval.
     */
    private double solveInterval(double x1, double y1, double x2, double y2) {
        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();
        final double relativeAccuracy = getRelativeAccuracy();

        double oldx = Double.POSITIVE_INFINITY;
        while (true) {
            final double x3 = calculateMidpoint(x1, x2);
            final double y3 = computeObjectiveValue(x3);
            if (isConvergedToZero(y3, functionValueAccuracy)) {
                return x3;
            }

            final double delta = calculateDelta(y1, y2, y3);
            final double correction = calculateCorrection(x1, x3, y2, y3, delta);
            final double x = x3 - correction;
            final double y = computeObjectiveValue(x);

            if (isConverged(x, oldx, relativeAccuracy, absoluteAccuracy) ||
                isConvergedToZero(y, functionValueAccuracy)) {
                return x;
            }

            // Prepare the new interval for the next iteration.
            // Ridders' method guarantees x1 < x < x2
            if (correction > 0.0) {             // x1 < x < x3
                if (isOppositeSign(y1, y)) {
                    x2 = x;
                    y2 = y;
                } else {
                    x1 = x;
                    x2 = x3;
                    y1 = y;
                    y2 = y3;
                }
            } else {                            // x3 < x < x2
                if (isOppositeSign(y2, y)) {
                    x1 = x;
                    y1 = y;
                } else {
                    x1 = x3;
                    x2 = x;
                    y1 = y3;
                    y2 = y;
                }
            }
            oldx = x;
        }
    }

    private double calculateMidpoint(double x1, double x2) {
        return 0.5 * (x1 + x2);
    }

    private double calculateDelta(double y1, double y2, double y3) {
        return 1.0 - (y1 * y2) / (y3 * y3);
    }

    private double calculateCorrection(double x1, double x3, double y2, double y3, double delta) {
        return (FastMath.signum(y2) * FastMath.signum(y3)) * (x3 - x1) / FastMath.sqrt(delta);
    }

    private boolean isConvergedToZero(double y, double functionValueAccuracy) {
        return FastMath.abs(y) <= functionValueAccuracy;
    }

    private boolean isConverged(double x, double oldx, double relativeAccuracy, double absoluteAccuracy) {
        final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
        return FastMath.abs(x - oldx) <= tolerance;
    }

    private boolean isOppositeSign(double yA, double yB) {
        return FastMath.signum(yA) + FastMath.signum(yB) == 0.0;
    }
}
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
        // [x1, x2] is the bracketing interval in each iteration
        // x3 is the midpoint of [x1, x2]
        // x is the new root approximation and an endpoint of the new interval
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

        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();
        final double relativeAccuracy = getRelativeAccuracy();

        double oldx = Double.POSITIVE_INFINITY;
        while (true) {
            // calculate the new root approximation
            final double x3 = 0.5 * (x1 + x2);
            final double y3 = computeObjectiveValue(x3);
            if (FastMath.abs(y3) <= functionValueAccuracy) {
                return x3;
            }
            final double correction = calculateCorrection(x1, x3, y1, y2, y3);
            final double x = x3 - correction;                // correction != 0
            final double y = computeObjectiveValue(x);

            // check for convergence
            if (hasConverged(x, oldx, relativeAccuracy, absoluteAccuracy, y, functionValueAccuracy)) {
                return x;
            }

            // prepare the new interval for next iteration
            // Ridders' method guarantees x1 < x < x2
            double[] updated = updateInterval(x1, y1, x2, y2, x3, y3, x, y, correction);
            x1 = updated[0];
            y1 = updated[1];
            x2 = updated[2];
            y2 = updated[3];

            oldx = x;
        }
    }

    /**
     * Calculates the correction term for Ridders' method.
     *
     * @param x1 Left boundary.
     * @param x3 Midpoint.
     * @param y1 Function value at left boundary.
     * @param y2 Function value at right boundary.
     * @param y3 Function value at midpoint.
     * @return The correction term.
     */
    private double calculateCorrection(double x1, double x3, double y1, double y2, double y3) {
        final double delta = 1 - (y1 * y2) / (y3 * y3);  // delta > 1 due to bracketing
        return (FastMath.signum(y2) * FastMath.signum(y3)) * (x3 - x1) / FastMath.sqrt(delta);
    }

    /**
     * Checks if the solver has converged to a root.
     *
     * @param x Current approximation.
     * @param oldx Previous approximation.
     * @param relativeAccuracy Relative accuracy tolerance.
     * @param absoluteAccuracy Absolute accuracy tolerance.
     * @param y Function value at current approximation.
     * @param functionValueAccuracy Function value accuracy tolerance.
     * @return True if converged, false otherwise.
     */
    private boolean hasConverged(double x, double oldx, double relativeAccuracy, double absoluteAccuracy, double y, double functionValueAccuracy) {
        final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
        return FastMath.abs(x - oldx) <= tolerance || FastMath.abs(y) <= functionValueAccuracy;
    }

    /**
     * Updates the bracketing interval bounds based on the Ridders' method logic.
     *
     * @param x1 Current left boundary.
     * @param y1 Function value at current left boundary.
     * @param x2 Current right boundary.
     * @param y2 Function value at current right boundary.
     * @param x3 Midpoint.
     * @param y3 Function value at midpoint.
     * @param x Current approximation.
     * @param y Function value at current approximation.
     * @param correction Correction value used to transition.
     * @return An array containing [new_x1, new_y1, new_x2, new_y2].
     */
    private double[] updateInterval(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y, double correction) {
        if (correction > 0.0) {             // x1 < x < x3
            if (FastMath.signum(y1) + FastMath.signum(y) == 0.0) {
                return new double[] { x1, y1, x, y };
            } else {
                return new double[] { x, y, x3, y3 };
            }
        } else {                            // x3 < x < x2
            if (FastMath.signum(y2) + FastMath.signum(y) == 0.0) {
                return new double[] { x, y, x2, y2 };
            } else {
                return new double[] { x3, y3, x, y };
            }
        }
    }
}
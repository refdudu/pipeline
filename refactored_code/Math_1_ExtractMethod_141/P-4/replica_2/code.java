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

        final double initialRoot = checkInitialRoots(min, max, y1, y2);
        if (!Double.isNaN(initialRoot)) {
            return initialRoot;
        }
        
        verifyBracketing(min, max);

        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();
        final double relativeAccuracy = getRelativeAccuracy();

        final IntervalState state = new IntervalState(min, y1, max, y2);

        while (true) {
            final double x3 = 0.5 * (state.x1 + state.x2);
            final double y3 = computeObjectiveValue(x3);
            if (isAbsoluteValueWithinAccuracy(y3, functionValueAccuracy)) {
                return x3;
            }
            
            final double delta = calculateDelta(state.y1, state.y2, y3);
            final double correction = calculateCorrection(state.x1, state.y2, y3, x3, delta);
            final double x = x3 - correction;
            final double y = computeObjectiveValue(x);

            if (isConverged(x, state.oldx, y, absoluteAccuracy, relativeAccuracy, functionValueAccuracy)) {
                return x;
            }

            updateInterval(state, x, y, x3, y3, correction);
            state.oldx = x;
        }
    }

    private double checkInitialRoots(double min, double max, double y1, double y2) {
        if (y1 == 0.0) {
            return min;
        }
        if (y2 == 0.0) {
            return max;
        }
        return Double.NaN;
    }

    private boolean isAbsoluteValueWithinAccuracy(double value, double accuracy) {
        return FastMath.abs(value) <= accuracy;
    }

    private double calculateDelta(double y1, double y2, double y3) {
        return 1.0 - (y1 * y2) / (y3 * y3);
    }

    private double calculateCorrection(double x1, double y2, double y3, double x3, double delta) {
        return (FastMath.signum(y2) * FastMath.signum(y3)) * (x3 - x1) / FastMath.sqrt(delta);
    }

    private boolean isConverged(double x, double oldx, double y, 
                                double absoluteAccuracy, double relativeAccuracy, 
                                double functionValueAccuracy) {
        final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
        return FastMath.abs(x - oldx) <= tolerance || FastMath.abs(y) <= functionValueAccuracy;
    }

    private void updateInterval(IntervalState state, double x, double y, double x3, double y3, double correction) {
        if (correction > 0.0) {
            if (FastMath.signum(state.y1) + FastMath.signum(y) == 0.0) {
                state.x2 = x;
                state.y2 = y;
            } else {
                state.x1 = x;
                state.x2 = x3;
                state.y1 = y;
                state.y2 = y3;
            }
        } else {
            if (FastMath.signum(state.y2) + FastMath.signum(y) == 0.0) { 
                state.x1 = x;
                state.y1 = y;
            } else {
                state.x1 = x3;
                state.x2 = x;
                state.y1 = y3;
                state.y2 = y;
            }
        }
    }

    private static class IntervalState {
        private double x1;
        private double y1;
        private double x2;
        private double y2;
        private double oldx;

        IntervalState(double x1, double y1, double x2, double y2) {
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
            this.oldx = Double.POSITIVE_INFINITY;
        } 
    }
}
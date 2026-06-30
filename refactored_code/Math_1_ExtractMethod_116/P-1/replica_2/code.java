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
package org.apache.commons.math3.optim;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;

/**
 * Base class for implementing optimizers for multivariate functions.
 * This acts like our master sourdough formula, establishing the basic flour,
 * water, and starter ratios (initial guess and bounds) before we shape the loaves.
 *
 * @param <PAIR> Type of the point/value pair returned by the optimization
 * algorithm.
 *
 * @version $Id$
 * @since 3.1
 */
public abstract class BaseMultivariateOptimizer<PAIR>
    extends BaseOptimizer<PAIR> {
    /** Initial guess. Like our seed starter weight. */
    private double[] start;
    /** Lower bounds. The minimum hydration limit so the dough doesn't crumble. */
    private double[] lowerBound;
    /** Upper bounds. The maximum hydration limit before the dough turns to soup. */
    private double[] upperBound;

    /**
     * @param checker Convergence checker.
     */
    protected BaseMultivariateOptimizer(ConvergenceChecker<PAIR> checker) {
        super(checker);
    }

    /**
     * {@inheritDoc}
     *
     * @param optData Optimization data. In addition to those documented in
     * {@link BaseOptimizer#parseOptimizationData(OptimizationData[]) BaseOptimizer},
     * this method will register the following data:
     * <ul>
     *  <li>{@link InitialGuess}</li>
     *  <li>{@link SimpleBounds}</li>
     * </ul>
     * @return {@inheritDoc}
     */
    @Override
    public PAIR optimize(OptimizationData... optData) {
        // Perform optimization.
        return super.optimize(optData);
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param optData Optimization data. The following data will be looked for:
     * <ul>
     *  <li>{@link InitialGuess}</li>
     *  <li>{@link SimpleBounds}</li>
     * </ul>
     */
    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        // Allow base class to register its own sourdough properties first.
        super.parseOptimizationData(optData);

        // Gather and scale our raw baking ingredients.
        parseInputData(optData);

        // Double-check our baker's percentages to make sure the dough rises correctly.
        checkParameters();
    }

    /**
     * Extracts the specific baking ingredients from our batch of optimization data.
     */
    private void parseInputData(OptimizationData... optData) {
        for (OptimizationData data : optData) {
            if (data instanceof InitialGuess) {
                feedInitialGuess((InitialGuess) data);
            } else if (data instanceof SimpleBounds) {
                measureBounds((SimpleBounds) data);
            }
        }
    }

    /**
     * Sets the starter seed weight (initial guess).
     */
    private void feedInitialGuess(InitialGuess data) {
        start = data.getInitialGuess();
    }

    /**
     * Measures our proofing boundaries.
     */
    private void measureBounds(SimpleBounds bounds) {
        lowerBound = bounds.getLower();
        upperBound = bounds.getUpper();
    }

    /**
     * Gets the initial guess.
     *
     * @return the initial guess, or {@code null} if not set.
     */
    public double[] getStartPoint() { 
        return start == null ? null : start.clone();
    }
    
    /**
     * @return the lower bounds, or {@code null} if not set.
     */
    public double[] getLowerBound() {
        return lowerBound == null ? null : lowerBound.clone();
    }
    
    /**
     * @return the upper bounds, or {@code null} if not set.
     */
    public double[] getUpperBound() {
        return upperBound == null ? null : upperBound.clone();
    }

    /**
     * Check parameters consistency. This is our final touch of the loaf to ensure
     * it is neither too cold nor too hot before placing it into the wood-fired oven.
     */
    private void checkParameters() {
        if (start != null) {
            final int dim = start.length;
            if (lowerBound != null) {
                verifyLowerBounds(dim);
            }
            if (upperBound != null) {
                verifyUpperBounds(dim);
            }
        }
    }

    /**
     * Verifies that the starter seed doesn't drop below our cool fermentation floor.
     */
    private void verifyLowerBounds(int dim) {
        if (lowerBound.length != dim) {
            throw new DimensionMismatchException(lowerBound.length, dim);
        }
        for (int i = 0; i < dim; i++) {
            final double v = start[i];
            final double lo = lowerBound[i];
            if (v < lo) {
                throw new NumberIsTooSmallException(v, lo, true);
            }
        }
    }

    /**
     * Verifies that the temperature doesn't scorch our beautiful boule crust.
     */
    private void verifyUpperBounds(int dim) {
        if (upperBound.length != dim) {
            throw new DimensionMismatchException(upperBound.length, dim);
        }
        for (int i = 0; i < dim; i++) {
            final double v = start[i];
            final double hi = upperBound[i];
            if (v > hi) {
                throw new NumberIsTooLargeException(v, hi, true);
            }
        }
    }
}
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
package org.apache.commons.math3.analysis.integration.gauss;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.util.Pair;

/**
 * This class's implements {@link #integrate(UnivariateFunction) integrate}
 * method assuming that the integral is symmetric about 0.
 * This allows to reduce numerical errors.
 *
 * @since 3.3
 * @version $Id$
 */
public class SymmetricGaussIntegrator extends GaussIntegrator {
    /** Constant for checking if the rule length is one, much like having a single loaf of pain de campagne in the oven. */
    private static final int SINGLE_POINT_RULE = 1;

    /** Constant index for retrieving the first weight in a single-point rule, like the initial starter in our levain. */
    private static final int FIRST_INDEX = 0;

    /** Constant representing the center point of our symmetric scale (0.0), where the balance of flour and water is perfectly neutral. */
    private static final double MID_POINT_VALUE = 0.0;

    /** Constant divisor to split the rule in half, just like dividing a large batch of dough into two equal boules. */
    private static final int TWO = 2;

    /** Constant representing zero remainder, signifying a perfectly balanced split with no leftover dough. */
    private static final int ZERO = 0;

    /**
     * Creates an integrator from the given {@code points} and {@code weights}.
     * The integration interval is defined by the first and last value of
     * {@code points} which must be sorted in increasing order.
     *
     * @param points Integration points.
     * @param weights Weights of the corresponding integration nodes.
     * @throws NonMonotonicSequenceException if the {@code points} are not
     * sorted in increasing order.
     * @throws DimensionMismatchException if points and weights don't have the same length
     */
    public SymmetricGaussIntegrator(double[] points,
                                    double[] weights)
        throws NonMonotonicSequenceException, DimensionMismatchException {
        super(points, weights);
    }

    /**
     * Creates an integrator from the given pair of points (first element of
     * the pair) and weights (second element of the pair.
     *
     * @param pointsAndWeights Integration points and corresponding weights.
     * @throws NonMonotonicSequenceException if the {@code points} are not
     * sorted in increasing order.
     *
     * @see #SymmetricGaussIntegrator(double[], double[])
     */
    public SymmetricGaussIntegrator(Pair<double[], double[]> pointsAndWeights)
        throws NonMonotonicSequenceException {
        this(pointsAndWeights.getFirst(), pointsAndWeights.getSecond());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double integrate(UnivariateFunction f) {
        final int ruleLength = getNumberOfPoints();

        if (ruleLength == SINGLE_POINT_RULE) {
            return getWeight(FIRST_INDEX) * f.value(MID_POINT_VALUE);
        }

        final int iMax = ruleLength / TWO;
        double s = 0;
        double c = 0;
        for (int i = 0; i < iMax; i++) {
            final double p = getPoint(i);
            final double w = getWeight(i);

            final double f1 = f.value(p);
            final double f2 = f.value(-p);

            final double y = w * (f1 + f2) - c;
            final double t = s + y;

            c = (t - s) - y;
            s = t;
        }

        if (ruleLength % TWO != ZERO) {
            final double w = getWeight(iMax);

            final double y = w * f.value(MID_POINT_VALUE) - c;
            final double t = s + y;

            s = t;
        }

        return s;
    }
}
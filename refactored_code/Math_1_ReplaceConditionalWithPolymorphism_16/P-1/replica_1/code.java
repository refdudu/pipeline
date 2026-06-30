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

package org.apache.commons.math3.distribution;

import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Implementation of the triangular real distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Triangular_distribution">
 * Triangular distribution (Wikipedia)</a>
 *
 * @version $Id$
 * @since 3.0
 */
public class TriangularDistribution extends AbstractRealDistribution {
    /** Serializable version identifier. */
    private static final long serialVersionUID = 20120112L;
    /** Lower limit of this distribution (inclusive). */
    private final double a;
    /** Upper limit of this distribution (inclusive). */
    private final double b;
    /** Mode of this distribution. */
    private final double c;
    /** Inverse cumulative probability accuracy. */
    private final double solverAbsoluteAccuracy;

    /** Distribution segments for density and cumulative probability calculations. */
    private final DistributionSegment[] segments;
    /** Probability segments for inverse cumulative probability calculations. */
    private final ProbabilitySegment[] pSegments;

    /**
     * Creates a triangular real distribution using the given lower limit,
     * upper limit, and mode.
     *
     * @param a Lower limit of this distribution (inclusive).
     * @param b Upper limit of this distribution (inclusive).
     * @param c Mode of this distribution.
     * @throws NumberIsTooLargeException if {@code a >= b} or if {@code c > b}.
     * @throws NumberIsTooSmallException if {@code c < a}.
     */
    public TriangularDistribution(double a, double c, double b)
        throws NumberIsTooLargeException, NumberIsTooSmallException {
        this(new Well19937c(), a, c, b);
    }

    /**
     * Creates a triangular distribution.
     *
     * @param rng Random number generator.
     * @param a Lower limit of this distribution (inclusive).
     * @param b Upper limit of this distribution (inclusive).
     * @param c Mode of this distribution.
     * @throws NumberIsTooLargeException if {@code a >= b} or if {@code c > b}.
     * @throws NumberIsTooSmallException if {@code c < a}.
     * @since 3.1
     */
    public TriangularDistribution(RandomGenerator rng,
                                  double a,
                                  double c,
                                  double b)
        throws NumberIsTooLargeException, NumberIsTooSmallException {
        super(rng);

        if (a >= b) {
            throw new NumberIsTooLargeException(
                            LocalizedFormats.LOWER_BOUND_NOT_BELOW_UPPER_BOUND,
                            a, b, false);
        }
        if (c < a) {
            throw new NumberIsTooSmallException(
                    LocalizedFormats.NUMBER_TOO_SMALL, c, a, true);
        }
        if (c > b) {
            throw new NumberIsTooLargeException(
                    LocalizedFormats.NUMBER_TOO_LARGE, c, b, true);
        }

        this.a = a;
        this.c = c;
        this.b = b;
        solverAbsoluteAccuracy = FastMath.max(FastMath.ulp(a), FastMath.ulp(b));

        this.segments = new DistributionSegment[] {
            new BelowASegment(),
            new BetweenAAndCSegment(),
            new AtCSegment(),
            new BetweenCAndBSegment(),
            new AboveBSegment()
        };

        this.pSegments = new ProbabilitySegment[] {
            new PZeroSegment(),
            new POneSegment(),
            new PLowerSegment(),
            new PUpperSegment()
        };
    }

    /**
     * Returns the mode {@code c} of this distribution.
     *
     * @return the mode {@code c} of this distribution
     */
    public double getMode() {
        return c;
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * For this distribution, the returned value is not really meaningful,
     * since exact formulas are implemented for the computation of the
     * {@link #inverseCumulativeProbability(double)} (no solver is invoked).
     * </p>
     * <p>
     * For lower limit {@code a} and upper limit {@code b}, the current
     * implementation returns {@code max(ulp(a), ulp(b)}.
     * </p>
     */
    @Override
    protected double getSolverAbsoluteAccuracy() {
        return solverAbsoluteAccuracy;
    }

    /**
     * {@inheritDoc}
     *
     * For lower limit {@code a}, upper limit {@code b} and mode {@code c}, the
     * PDF is given by
     * <ul>
     * <li>{@code 2 * (x - a) / [(b - a) * (c - a)]} if {@code a <= x < c},</li>
     * <li>{@code 2 / (b - a)} if {@code x = c},</li>
     * <li>{@code 2 * (b - x) / [(b - a) * (b - c)]} if {@code c < x <= b},</li>
     * <li>{@code 0} otherwise.
     * </ul>
     */
    public double density(double x) {
        for (DistributionSegment segment : segments) {
            if (segment.contains(x)) {
                return segment.density(x);
            }
        }
        return 0;
    }

    /**
     * {@inheritDoc}
     *
     * For lower limit {@code a}, upper limit {@code b} and mode {@code c}, the
     * CDF is given by
     * <ul>
     * <li>{@code 0} if {@code x < a},</li>
     * <li>{@code (x - a)^2 / [(b - a) * (c - a)]} if {@code a <= x < c},</li>
     * <li>{@code (c - a) / (b - a)} if {@code x = c},</li>
     * <li>{@code 1 - (b - x)^2 / [(b - a) * (b - c)]} if {@code c < x <= b},</li>
     * <li>{@code 1} if {@code x > b}.</li>
     * </ul>
     */
    public double cumulativeProbability(double x)  {
        for (DistributionSegment segment : segments) {
            if (segment.contains(x)) {
                return segment.cumulativeProbability(x);
            }
        }
        return 1;
    }

    /**
     * {@inheritDoc}
     *
     * For lower limit {@code a}, upper limit {@code b}, and mode {@code c},
     * the mean is {@code (a + b + c) / 3}.
     */
    public double getNumericalMean() {
        return (a + b + c) / 3;
    }

    /**
     * {@inheritDoc}
     *
     * For lower limit {@code a}, upper limit {@code b}, and mode {@code c},
     * the variance is {@code (a^2 + b^2 + c^2 - a * b - a * c - b * c) / 18}.
     */
    public double getNumericalVariance() {
        return (a * a + b * b + c * c - a * b - a * c - b * c) / 18;
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is equal to the lower limit parameter
     * {@code a} of the distribution.
     *
     * @return lower bound of the support
     */
    public double getSupportLowerBound() {
        return a;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is equal to the upper limit parameter
     * {@code b} of the distribution.
     *
     * @return upper bound of the support
     */
    public double getSupportUpperBound() {
        return b;
    }

    /** {@inheritDoc} */
    public boolean isSupportLowerBoundInclusive() {
        return true;
    }

    /** {@inheritDoc} */
    public boolean isSupportUpperBoundInclusive() {
        return true;
    }

    /**
     * {@inheritDoc}
     *
     * The support of this distribution is connected.
     *
     * @return {@code true}
     */
    public boolean isSupportConnected() {
        return true;
    }

    @Override
    public double inverseCumulativeProbability(double p)
        throws OutOfRangeException {
        if (p < 0 || p > 1) {
            throw new OutOfRangeException(p, 0, 1);
        }
        for (ProbabilitySegment segment : pSegments) {
            if (segment.contains(p)) {
                return segment.inverseCumulativeProbability(p);
            }
        }
        return a;
    }

    /**
     * Strategy interface representing a segment of the distribution's range.
     */
    private interface DistributionSegment {
        boolean contains(double x);
        double density(double x);
        double cumulativeProbability(double x);
    }

    /**
     * Segment for values of x strictly less than 'a'.
     */
    private class BelowASegment implements DistributionSegment {
        public boolean contains(double x) {
            return x < a;
        }
        public double density(double x) {
            return 0;
        }
        public double cumulativeProbability(double x) {
            return 0;
        }
    }

    /**
     * Segment for values of x between 'a' (inclusive) and 'c' (exclusive).
     */
    private class BetweenAAndCSegment implements DistributionSegment {
        public boolean contains(double x) {
            return a <= x && x < c;
        }
        public double density(double x) {
            double divident = 2 * (x - a);
            double divisor = (b - a) * (c - a);
            return divident / divisor;
        }
        public double cumulativeProbability(double x) {
            double divident = (x - a) * (x - a);
            double divisor = (b - a) * (c - a);
            return divident / divisor;
        }
    }

    /**
     * Segment for x exactly equal to the mode 'c'.
     */
    private class AtCSegment implements DistributionSegment { 
        public boolean contains(double x) {
            return x == c;
        }
        public double density(double x) {
            return 2 / (b - a);
        }
        public double cumulativeProbability(double x) {
            return (c - a) / (b - a);
        }
    }

    /**
     * Segment for values of x between 'c' (exclusive) and 'b' (inclusive).
     */
    private class BetweenCAndBSegment implements DistributionSegment {
        public boolean contains(double x) {
            return c < x && x <= b;
        }
        public double density(double x) {
            double divident = 2 * (b - x);
            double divisor = (b - a) * (b - c);
            return divident / divisor;
        }
        public double cumulativeProbability(double x) {
            double divident = (b - x) * (b - x);
            double divisor = (b - a) * (b - c);
            return 1 - (divident / divisor);
        }
    }

    /**
     * Segment for values of x strictly greater than 'b'.
     */
    private class AboveBSegment implements DistributionSegment {
        public boolean contains(double x) {
            return x > b;
        }
        public double density(double x) {
            return 0;
        }
        public double cumulativeProbability(double x) {
            return 1;
        }
    }

    /**
     * Strategy interface representing a segment of the probability range.
     */
    private interface ProbabilitySegment {
        boolean contains(double p);
        double inverseCumulativeProbability(double p);
    }

    /**
     * Segment for probability strictly equal to 0.
     */
    private class PZeroSegment implements ProbabilitySegment {
        public boolean contains(double p) {
            return p == 0;
        }
        public double inverseCumulativeProbability(double p) {
            return a;
        }
    }

    /**
     * Segment for probability strictly equal to 1.
     */
    private class POneSegment implements ProbabilitySegment {
        public boolean contains(double p) {
            return p == 1;
        }
        public double inverseCumulativeProbability(double p) {
            return b;
        }
    }

    /**
     * Segment for lower range probabilities.
     */
    private class PLowerSegment implements ProbabilitySegment {
        public boolean contains(double p) {
            return p > 0 && p < (c - a) / (b - a);
        }
        public double inverseCumulativeProbability(double p) {
            return a + FastMath.sqrt(p * (b - a) * (c - a));
        }
    }

    /**
     * Segment for upper range probabilities.
     */
    private class PUpperSegment implements ProbabilitySegment {
        public boolean contains(double p) {
            return p >= (c - a) / (b - a) && p < 1;
        }
        public double inverseCumulativeProbability(double p) {
            return b - FastMath.sqrt((1 - p) * (b - a) * (b - c));
        }
    }
}
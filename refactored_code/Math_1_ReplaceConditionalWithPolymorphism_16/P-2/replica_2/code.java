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

    /** Polymorphic strategy list for interval evaluations. */
    private transient Interval[] intervals;
    /** Polymorphic strategy list for probability range evaluations. */
    private transient ProbabilityRange[] probabilityRanges;

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
        initStrategies();
    }

    /**
     * Initializes the state-dependent interval and probability range strategy helpers.
     */
    private void initStrategies() {
        this.intervals = new Interval[] {
            new BelowLowerBound(),
            new LowerSubInterval(),
            new AtMode(),
            new UpperSubInterval(),
            new AboveUpperBound()
        };
        this.probabilityRanges = new ProbabilityRange[] {
            new AtZero(),
            new AtOne(),
            new BelowMode(),
            new AboveMode()
        };
    }

    /**
     * Re-initializes non-serializable strategies upon deserialization.
     */
    private void readObject(java.io.ObjectInputStream in)
        throws java.io.IOException, ClassNotFoundException {
        in.defaultReadObject();
        initStrategies();
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
     */
    @Override
    protected double getSolverAbsoluteAccuracy() {
        return solverAbsoluteAccuracy;
    }

    /**
     * {@inheritDoc}
     */
    public double density(double x) {
        for (Interval interval : intervals) {
            if (interval.contains(x)) {
                return interval.density(x);
            }
        }
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    public double cumulativeProbability(double x) {
        for (Interval interval : intervals) {
            if (interval.contains(x)) {
                return interval.cumulativeProbability(x);
            }
        }
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    public double getNumericalMean() {
        return (a + b + c) / 3;
    }

    /**
     * {@inheritDoc}
     */
    public double getNumericalVariance() {
        return (a * a + b * b + c * c - a * b - a * c - b * c) / 18;
    }

    /**
     * {@inheritDoc}
     */
    public double getSupportLowerBound() {
        return a;
    }

    /**
     * {@inheritDoc}
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
     */
    public boolean isSupportConnected() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double inverseCumulativeProbability(double p)
        throws OutOfRangeException {
        if (p < 0 || p > 1) {
            throw new OutOfRangeException(p, 0, 1);
        }
        for (ProbabilityRange range : probabilityRanges) {
            if (range.contains(p)) {
                return range.inverseCumulativeProbability(p);
            }
        }
        return b;
    }

    private interface Interval {
        boolean contains(double x);
        double density(double x);
        double cumulativeProbability(double x);
    }

    private class BelowLowerBound implements Interval {
        @Override
        public boolean contains(double x) {
            return x < a;
        }

        @Override
        public double density(double x) {
            return 0;
        }

        @Override
        public double cumulativeProbability(double x) {
            return 0;
        }
    }

    private class LowerSubInterval implements Interval {
        @Override
        public boolean contains(double x) {
            return a <= x && x < c;
        }

        @Override
        public double density(double x) {
            double divident = 2 * (x - a);
            double divisor = (b - a) * (c - a);
            return divident / divisor;
        }

        @Override
        public double cumulativeProbability(double x) {
            double divident = (x - a) * (x - a);
            double divisor = (b - a) * (c - a);
            return divident / divisor;
        }
    }

    private class AtMode implements Interval {
        @Override
        public boolean contains(double x) {
            return x == c;
        }

        @Override
        public double density(double x) {
            return 2 / (b - a);
        }

        @Override
        public double cumulativeProbability(double x) {
            return (c - a) / (b - a);
        }
    }

    private class UpperSubInterval implements Interval {
        @Override
        public boolean contains(double x) {
            return c < x && x <= b;
        }

        @Override
        public double density(double x) {
            double divident = 2 * (b - x);
            double divisor = (b - a) * (b - c);
            return divident / divisor;
        }

        @Override
        public double cumulativeProbability(double x) {
            double divident = (b - x) * (b - x);
            double divisor = (b - a) * (b - c);
            return 1 - (divident / divisor);
        }
    }

    private class AboveUpperBound implements Interval {
        @Override
        public boolean contains(double x) {
            return x > b;
        }

        @Override
        public double density(double x) {
            return 0;
        }

        @Override
        public double cumulativeProbability(double x) {
            return 1;
        }
    }

    private interface ProbabilityRange {
        boolean contains(double p);
        double inverseCumulativeProbability(double p);
    }

    private class AtZero implements ProbabilityRange {
        @Override
        public boolean contains(double p) {
            return p == 0;
        }

        @Override
        public double inverseCumulativeProbability(double p) {
            return a;
        }
    }

    private class AtOne implements ProbabilityRange {
        @Override
        public boolean contains(double p) {
            return p == 1;
        }

        @Override
        public double inverseCumulativeProbability(double p) {
            return b;
        }
    }

    private class BelowMode implements ProbabilityRange {
        @Override
        public boolean contains(double p) {
            return p < (c - a) / (b - a);
        }

        @Override
        public double inverseCumulativeProbability(double p) {
            return a + FastMath.sqrt(p * (b - a) * (c - a));
        }
    }

    private class AboveMode implements ProbabilityRange {
        @Override
        public boolean contains(double p) {
            return true;
        }

        @Override
        public double inverseCumulativeProbability(double p) {
            return b - FastMath.sqrt((1 - p) * (b - a) * (b - c));
        }
    }
}
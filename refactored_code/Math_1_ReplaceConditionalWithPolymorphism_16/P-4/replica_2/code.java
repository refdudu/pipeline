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

    /** Regions of the distribution. */
    private final Region[] regions;
    /** Probability segments for inverse cumulative probability calculation. */
    private final ProbabilitySegment[] segments;

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

        this.regions = new Region[] {
            new BelowSupport(),
            new LowerSlope(),
            new Peak(),
            new UpperSlope(),
            new AboveSupport()
        };

        this.segments = new ProbabilitySegment[] {
            new ZeroProbability(),
            new OneProbability(),
            new LeftSegment(),
            new RightSegment()
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
     */
    @Override
    protected double getSolverAbsoluteAccuracy() {
        return solverAbsoluteAccuracy;
    }

    /**
     * {@inheritDoc}
     */
    public double density(double x) {
        for (Region region : regions) {
            if (region.contains(x)) {
                return region.density(x);
            } // end if
        } // end for
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    public double cumulativeProbability(double x)  {
        for (Region region : regions) {
            if (region.contains(x)) {
                return region.cumulativeProbability(x);
            } // end if
        } // end for
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

    @Override
    public double inverseCumulativeProbability(double p)
        throws OutOfRangeException {
        if (p < 0 || p > 1) {
            throw new OutOfRangeException(p, 0, 1);
        }
        for (ProbabilitySegment segment : segments) {
            if (segment.contains(p)) {
                return segment.inverseCumulativeProbability(p);
            } // end if
        } // end for
        return b;
    }

    private abstract class Region {
        abstract boolean contains(double x);
        abstract double density(double x);
        abstract double cumulativeProbability(double x);
    }

    private final class BelowSupport extends Region {
        @Override
        boolean contains(double x) {
            return x < a;
        }

        @Override
        double density(double x) {
            return 0;
        }

        @Override
        double cumulativeProbability(double x) {
            return 0;
        }
    }

    private final class LowerSlope extends Region {
        @Override
        boolean contains(double x) {
            return a <= x && x < c;
        }

        @Override
        double density(double x) {
            return (2 * (x - a)) / ((b - a) * (c - a));
        }

        @Override
        double cumulativeProbability(double x) {
            return ((x - a) * (x - a)) / ((b - a) * (c - a));
        }
    }

    private final class Peak extends Region {
        @Override
        boolean contains(double x) {
            return x == c;
        }

        @Override
        double density(double x) {
            return 2 / (b - a);
        }

        @Override
        double cumulativeProbability(double x) {
            return (c - a) / (b - a);
        }
    }

    private final class UpperSlope extends Region {
        @Override
        boolean contains(double x) {
            return c < x && x <= b;
        }

        @Override
        double density(double x) {
            return (2 * (b - x)) / ((b - a) * (b - c));
        }

        @Override
        double cumulativeProbability(double x) {
            return 1 - ((b - x) * (b - x)) / ((b - a) * (b - c));
        }
    }

    private final class AboveSupport extends Region {
        @Override
        boolean contains(double x) {
            return x > b;
        }

        @Override
        double density(double x) {
            return 0;
        }

        @Override
        double cumulativeProbability(double x) {
            return 1;
        }
    }

    private abstract class ProbabilitySegment {
        abstract boolean contains(double p);
        abstract double inverseCumulativeProbability(double p);
    }

    private final class ZeroProbability extends ProbabilitySegment {
        @Override
        boolean contains(double p) {
            return p == 0;
        }

        @Override
        double inverseCumulativeProbability(double p) {
            return a;
        }
    }

    private final class OneProbability extends ProbabilitySegment {
        @Override
        boolean contains(double p) {
            return p == 1;
        }

        @Override
        double inverseCumulativeProbability(double p) {
            return b;
        }
    }

    private final class LeftSegment extends ProbabilitySegment {
        @Override
        boolean contains(double p) {
            return p < (c - a) / (b - a);
        }

        @Override
        double inverseCumulativeProbability(double p) {
            return a + FastMath.sqrt(p * (b - a) * (c - a));
        }
    }

    private final class RightSegment extends ProbabilitySegment {
        @Override
        boolean contains(double p) {
            return true;
        }

        @Override
        double inverseCumulativeProbability(double p) {
            return b - FastMath.sqrt((1 - p) * (b - a) * (b - c));
        }
    }
}
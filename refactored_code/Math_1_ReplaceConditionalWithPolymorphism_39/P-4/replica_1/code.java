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

import java.io.Serializable;

import org.apache.commons.math3.exception.MathInternalError;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomDataImpl;
import org.apache.commons.math3.util.FastMath;

/**
 * Base class for integer-valued discrete distributions.  Default
 * implementations are provided for some of the methods that do not vary
 * from distribution to distribution.
 *
 * @version $Id$
 */
public abstract class AbstractIntegerDistribution implements IntegerDistribution, Serializable {

    /** Serializable version identifier */
    private static final long serialVersionUID = -1146319659338487221L;

    /**
     * RandomData instance used to generate samples from the distribution.
     * @deprecated As of 3.1, to be removed in 4.0. Please use the
     * {@link #random} instance variable instead.
     */
    @Deprecated
    protected final RandomDataImpl randomData = new RandomDataImpl();

    /**
     * RNG instance used to generate samples from the distribution.
     * @since 3.1
     */
    protected final RandomGenerator random;

    /**
     * @deprecated As of 3.1, to be removed in 4.0. Please use
     * {@link #AbstractIntegerDistribution(RandomGenerator)} instead.
     */
    @Deprecated
    protected AbstractIntegerDistribution() {
        // Legacy users are only allowed to access the deprecated "randomData".
        // New users are forbidden to use this constructor.
        random = null;
    }

    /**
     * @param rng Random number generator.
     * @since 3.1
     */
    protected AbstractIntegerDistribution(RandomGenerator rng) {
        random = rng;
    }

    /**
     * {@inheritDoc}
     *
     * The default implementation uses the identity
     * <p>{@code P(x0 < X <= x1) = P(X <= x1) - P(X <= x0)}</p>
     */
    public double cumulativeProbability(int x0, int x1) throws NumberIsTooLargeException {
        if (x1 < x0) {
            throw new NumberIsTooLargeException(LocalizedFormats.LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT,
                    x0, x1, true);
        }
        return cumulativeProbability(x1) - cumulativeProbability(x0);
    }

    /**
     * {@inheritDoc}
     *
     * The default implementation returns
     * <ul>
     * <li>{@link #getSupportLowerBound()} for {@code p = 0},</li>
     * <li>{@link #getSupportUpperBound()} for {@code p = 1}, and</li>
     * <li>{@link #solveInverseCumulativeProbability(double, int, int)} for
     *     {@code 0 < p < 1}.</li>
     * </ul>
     */
    public int inverseCumulativeProbability(final double p) throws OutOfRangeException {
        if (p < 0.0 || p > 1.0) {
            throw new OutOfRangeException(p, 0, 1);
        }

        int lower = getSupportLowerBound();
        if (p == 0.0) {
            return lower;
        }

        final LowerBoundStrategy lbStrategy = getLowerBoundStrategy(lower);
        if (lbStrategy.isCompleted(p, lower)) {
            return lower;
        }
        lower = lbStrategy.getAdjustedLower(lower);

        int upper = getSupportUpperBound();
        if (p == 1.0) {
            return upper;
        }

        final double mu = getNumericalMean();
        final double sigma = FastMath.sqrt(getNumericalVariance());
        final ChebyshevStrategy chebyshevStrategy = getChebyshevStrategy(mu, sigma);
        final Bracket bracket = chebyshevStrategy.narrow(p, mu, sigma, lower, upper);
        lower = bracket.getLower();
        upper = bracket.getUpper();

        return solveInverseCumulativeProbability(p, lower, upper);
    }

    /**
     * Strategy interface for handling lower bound check and adjustment based on its initial value.
     */
    private interface LowerBoundStrategy {
        int getAdjustedLower(int lower);
        boolean isCompleted(double p, int lower);
    }

    /**
     * Strategy implementation for the case where the lower bound equals {@code Integer.MIN_VALUE}.
     */
    private class MinValueLowerBoundStrategy implements LowerBoundStrategy {
        public int getAdjustedLower(int lower) {
            return lower;
        }

        public boolean isCompleted(double p, int lower) {
            return checkedCumulativeProbability(lower) >= p;
        }
    }

    /**
     * Strategy implementation for cases where the lower bound is greater than {@code Integer.MIN_VALUE}.
     */
    private static class NormalLowerBoundStrategy implements LowerBoundStrategy {
        public int getAdjustedLower(int lower) {
            return lower - 1;
        }

        public boolean isCompleted(double p, int lower) {
            return false;
        }
    }

    /**
     * Factory method to select the appropriate {@link LowerBoundStrategy}.
     */
    private LowerBoundStrategy getLowerBoundStrategy(int lower) {
        if (lower == Integer.MIN_VALUE) {
            return new MinValueLowerBoundStrategy();
        }
        return new NormalLowerBoundStrategy();
    }

    /**
     * Represents a bracket range [lower, upper].
     */
    private static class Bracket {
        private final int lower;
        private final int upper;

        public Bracket(int lower, int upper) {
            this.lower = lower;
            this.upper = upper;
        }

        public int getLower() {
            return lower;
        }

        public int getUpper() {
            return upper;
        }
    }

    /**
     * Strategy interface for applying the Chebyshev inequality narrow down.
     */
    private interface ChebyshevStrategy {
        Bracket narrow(double p, double mu, double sigma, int lower, int upper);
    }

    /**
     * Active strategy where Chebyshev inequality is applicable.
     */
    private static class ActiveChebyshevStrategy implements ChebyshevStrategy {
        public Bracket narrow(double p, double mu, double sigma, int lower, int upper) {
            double k = FastMath.sqrt((1.0 - p) / p);
            double tmp = mu - k * sigma;
            int newLower = lower;
            if (tmp > lower) {
                newLower = ((int) Math.ceil(tmp)) - 1;
            }
            k = 1.0 / k;
            tmp = mu + k * sigma;
            int newUpper = upper;
            if (tmp < upper) {
                newUpper = ((int) Math.ceil(tmp)) - 1;
            }
            return new Bracket(newLower, newUpper);
        }
    }

    /**
     * Inactive strategy where Chebyshev inequality cannot be applied.
     */
    private static class InactiveChebyshevStrategy implements ChebyshevStrategy {
        public Bracket narrow(double p, double mu, double sigma, int lower, int upper) {
            return new Bracket(lower, upper);
        }
    }

    /**
     * Factory method to select the appropriate {@link ChebyshevStrategy}.
     */
    private ChebyshevStrategy getChebyshevStrategy(double mu, double sigma) {
        final boolean chebyshevApplies = !(Double.isInfinite(mu) || Double.isNaN(mu) ||
                Double.isInfinite(sigma) || Double.isNaN(sigma) || sigma == 0.0);
        if (chebyshevApplies) {
            return new ActiveChebyshevStrategy();
        }
        return new InactiveChebyshevStrategy();
    }

    /**
     * This is a utility function used by {@link
     * #inverseCumulativeProbability(double)}. It assumes {@code 0 < p < 1} and
     * that the inverse cumulative probability lies in the bracket {@code
     * (lower, upper]}. The implementation does simple bisection to find the
     * smallest {@code p}-quantile <code>inf{x in Z | P(X<=x) >= p}</code>.
     *
     * @param p the cumulative probability
     * @param lower a value satisfying {@code cumulativeProbability(lower) < p}
     * @param upper a value satisfying {@code p <= cumulativeProbability(upper)}
     * @return the smallest {@code p}-quantile of this distribution
     */
    protected int solveInverseCumulativeProbability(final double p, int lower, int upper) {
        while (lower + 1 < upper) {
            int xm = (lower + upper) / 2;
            if (xm < lower || xm > upper) {
                /*
                 * Overflow.
                 * There will never be an overflow in both calculation methods
                 * for xm at the same time
                 */
                xm = lower + (upper - lower) / 2;
            }

            double pm = checkedCumulativeProbability(xm);
            if (pm >= p) {
                upper = xm;
            } else {
                lower = xm;
            }
        }
        return upper;
    }

    /** {@inheritDoc} */
    public void reseedRandomGenerator(long seed) {
        random.setSeed(seed);
        randomData.reSeed(seed);
    }

    /**
     * {@inheritDoc}
     *
     * The default implementation uses the
     * <a href="http://en.wikipedia.org/wiki/Inverse_transform_sampling">
     * inversion method</a>.
     */
    public int sample() {
        return inverseCumulativeProbability(random.nextDouble());
    }

    /**
     * {@inheritDoc}
     *
     * The default implementation generates the sample by calling
     * {@link #sample()} in a loop.
     */
    public int[] sample(int sampleSize) {
        if (sampleSize <= 0) {
            throw new NotStrictlyPositiveException(
                    LocalizedFormats.NUMBER_OF_SAMPLES, sampleSize);
        }
        int[] out = new int[sampleSize];
        for (int i = 0; i < sampleSize; i++) {
            out[i] = sample();
        }
        return out;
    }

    /**
     * Computes the cumulative probability function and checks for {@code NaN}
     * values returned. Throws {@code MathInternalError} if the value is
     * {@code NaN}. Rethrows any exception encountered evaluating the cumulative
     * probability function. Throws {@code MathInternalError} if the cumulative
     * probability function returns {@code NaN}.
     *
     * @param argument input value
     * @return the cumulative probability
     * @throws MathInternalError if the cumulative probability is {@code NaN}
     */
    private double checkedCumulativeProbability(int argument)
        throws MathInternalError {
        double result = Double.NaN;
        result = cumulativeProbability(argument);
        if (Double.isNaN(result)) {
            throw new MathInternalError(LocalizedFormats
                    .DISCRETE_CUMULATIVE_PROBABILITY_RETURNED_NAN, argument);
        }
        return result;
    }
}
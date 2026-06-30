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
package org.apache.commons.math3.primes;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

/**
 * Implementation of the Pollard's rho factorization algorithm.
 * @version $Id$
 * @since 3.2
 */
class PollardRho {

    /**
     * Hide utility class.
     */
    private PollardRho() {
    }

    /**
     * Factorization using Pollard's rho algorithm.
     * @param n number to factors, must be &gt; 0
     * @return the list of prime factors of n.
     */
    public static List<Integer> primeFactors(int n) {
        final List<Integer> factors = new ArrayList<Integer>();

        n = SmallPrimes.smallTrialDivision(n, factors);
        if (isFullyFactorized(n)) {
            return factors;
        }

        if (SmallPrimes.millerRabinPrimeTest(n)) {
            factors.add(n);
            return factors;
        }

        addRhoFactors(n, factors);
        return factors;
    }

    private static boolean isFullyFactorized(int n) {
        return n == 1;
    }

    private static void addRhoFactors(int n, List<Integer> factors) {
        final int divisor = rhoBrent(n);
        factors.add(divisor);
        factors.add(n / divisor);
    }

    /**
     * Implementation of the Pollard's rho factorization algorithm.
     * <p>
     * This implementation follows the paper "An improved Monte Carlo factorization algorithm"
     * by Richard P. Brent. This avoids the triple computation of f(x) typically found in Pollard's
     * rho implementations. It also batches several gcd computation into 1.
     * <p>
     * The backtracking is not implemented as we deal only with semi-primes.
     * 
     * @param n number to factor, must be semi-prime.
     * @return a prime factor of n.
     */
    static int rhoBrent(final int n) {
        return new BrentSearch(n).findFactor();
    }

    /**
     * Gcd between two positive numbers.
     * <p>
     * Gets the greatest common divisor of two numbers, using the "binary gcd" method,
     * which avoids division and modulo operations. See Knuth 4.5.2 algorithm B.
     * This algorithm is due to Josef Stein (1961).
     * </p>
     * Special cases:
     * <ul>
     * <li>The result of {@code gcd(x, x)}, {@code gcd(0, x)} and {@code gcd(x, 0)} is the value of {@code x}.</li>
     * <li>The invocation {@code gcd(0, 0)} is the only one which returns {@code 0}.</li>
     * </ul>
     *
     * @param a first number, must be &ge; 0
     * @param b second number, must be &ge; 0
     * @return gcd(a,b)
     */
    static int gcdPositive(int a, int b) {
        if (a == 0) {
            return b;
        } else if (b == 0) {
            return a;
        }

        final int aTwos = Integer.numberOfTrailingZeros(a);
        final int bTwos = Integer.numberOfTrailingZeros(b);
        final int shift = FastMath.min(aTwos, bTwos);

        a >>= aTwos;
        b >>= bTwos;

        return performBinaryGcd(a, b) << shift;
    }

    private static int performBinaryGcd(int a, int b) {
        while (a != b) {
            final int delta = a - b;
            b = FastMath.min(a, b);
            a = FastMath.abs(delta);
            a >>= Integer.numberOfTrailingZeros(a);
        }
        return a;
    }

    /**
     * Helper class representing the state and steps of Brent's search algorithm.
     */
    private static class BrentSearch {
        private static final int X0 = 2;
        private static final int M = 25;

        private final int n;
        private int cst = SmallPrimes.PRIMES_LAST;
        private int y = X0;
        private int r = 1;
        private int x = X0;

        BrentSearch(int n) {
            this.n = n;
        }

        int findFactor() {
            do {
                x = y;
                y = hopY(r);
                int k = 0;
                do {
                    final int bound = FastMath.min(M, r - k);
                    final int factor = evaluateBound(bound);
                    if (factor != 0 && factor != -1) {
                        return factor;
                    }
                    if (factor == -1) {
                        k = -M;
                    }
                    k = k + M;
                } while (k < r);
                r = 2 * r;
            } while (true);
        }

        private int hopY(int steps) {
            int currentY = y;
            for (int i = 0; i < steps; i++) {
                currentY = nextY(currentY);
            }
            return currentY;
        }

        private int nextY(int currentY) {
            final long y2 = ((long) currentY) * currentY;
            return (int) ((y2 + cst) % n);
        }

        private int evaluateBound(int bound) {
            int q = 1;
            boolean resetTriggered = false;
            for (int i = -3; i < bound; i++) {
                y = nextY(y);
                final long divisor = FastMath.abs(x - y);
                if (0 == divisor) {
                    resetTriggered = true;
                    break;
                }
                final long prod = divisor * q;
                q = (int) (prod % n);
                if (0 == q) {
                    return gcdPositive(FastMath.abs((int) divisor), n);
                }
            }

            final int out = gcdPositive(FastMath.abs(q), n);
            if (1 != out) {
                return out;
            }

            if (resetTriggered) {
                cst += SmallPrimes.PRIMES_LAST;
                y = X0;
                r = 1;
                return -1;
            }

            return 0;
        }
    }
}
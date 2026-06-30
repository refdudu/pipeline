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
        if (1 == n) {
            return factors;
        }

        if (SmallPrimes.millerRabinPrimeTest(n)) {
            factors.add(n);
        } else {
            addFactorsOfSemiPrime(n, factors);
        }
        return factors;
    }

    /**
     * Extracts factors of a semi-prime and adds them to the list of factors.
     * 
     * @param n the semi-prime number to factor
     * @param factors the list of factors to add to
     */
    private static void addFactorsOfSemiPrime(int n, List<Integer> factors) {
        int divisor = rhoBrent(n);
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
        final int x0 = 2;
        final int m = 25;
        int cst = SmallPrimes.PRIMES_LAST;
        int y = x0;
        int r = 1;
        do {
            int x = y;
            y = updateY(y, r, cst, n);
            int k = 0;
            do {
                final int bound = FastMath.min(m, r - k);
                int q = 1;
                for (int i = -3; i < bound; i++) { //start at -3 to ensure we enter this loop at least 3 times
                    y = nextY(y, cst, n);
                    final long divisor = FastMath.abs(x - y);
                    if (0 == divisor) {
                        cst += SmallPrimes.PRIMES_LAST;
                        k = -m;
                        y = x0;
                        r = 1;
                        break;
                    }
                    q = updateProduct(divisor, q, n);
                    if (0 == q) {
                        return gcdPositive(FastMath.abs((int) divisor), n);
                    }
                }
                final int out = gcdPositive(FastMath.abs(q), n);
                if (1 != out) {
                    return out;
                }
                k = k + m;
            } while (k < r);
            r = 2 * r;
        } while (true);
    }

    /**
     * Computes the next y value in the pseudo-random sequence.
     */
    private static int nextY(int y, int cst, int n) {
        final long y2 = ((long) y) * y;
        return (int) ((y2 + cst) % n);
    }

    /**
     * Updates y by applying the pseudo-random step multiple times.
     */
    private static int updateY(int y, int steps, int cst, int n) {
        for (int i = 0; i < steps; i++) {
            y = nextY(y, cst, n);
        }
        return y;
    }

    /**
     * Updates the product term used in the Brent's Pollard-rho gcd batching.
     */
    private static int updateProduct(long divisor, int q, int n) {
        final long prod = divisor * q;
        return (int) (prod % n);
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
    static int gcdPositive(int a, int b){
        // both a and b must be positive, it is not checked here
        // gdc(a,0) = a
        if (a == 0) {
            return b;
        } else if (b == 0) {
            return a;
        }

        // make a and b odd, keep in mind the common power of twos
        final int shift = FastMath.min(Integer.numberOfTrailingZeros(a), Integer.numberOfTrailingZeros(b));
        a = removeTrailingZeros(a);
        b = removeTrailingZeros(b);

        // a and b >0
        // if a > b then gdc(a,b) = gcd(a-b,b)
        // if a < b then gcd(a,b) = gcd(b-a,a)
        // so next a is the absolute difference and next b is the minimum of current values
        while (a != b) {
            final int delta = a - b;
            b = FastMath.min(a, b);
            a = removeTrailingZeros(FastMath.abs(delta));
        }

        // gcd(a,a) = a, just "add" the common power of twos
        return a << shift;
    }

    /**
     * Helper to remove trailing zeros of a number by right-shifting.
     */
    private static int removeTrailingZeros(int value) {
        return value >> Integer.numberOfTrailingZeros(value);
    }
}
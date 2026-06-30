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
            return factors;
        }

        return findFactorsOfSemiPrime(n, factors);
    }

    private static List<Integer> findFactorsOfSemiPrime(int n, List<Integer> factors) {
        int divisor = rhoBrent(n);
        factors.add(divisor);
        factors.add(n / divisor);
        return factors;
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
        final BrentState state = new BrentState(n);
        do {
            state.updateX();
            state.advanceY(state.r);
            
            final int factor = checkLoopForFactor(state, n);
            if (factor != 0) {
                return factor;
            }
            state.r = 2 * state.r;
        } while (true);
    }

    private static int checkLoopForFactor(BrentState state, int n) {
        int k = 0;
        do {
            final int bound = FastMath.min(BrentState.M, state.r - k);
            int q = 1;
            for (int i = -3; i < bound; i++) {
                state.y = state.nextValue(state.y);
                final long divisor = FastMath.abs(state.x - state.y);
                if (0 == divisor) {
                    state.reset();
                    k = -BrentState.M;
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
            k = k + BrentState.M;
        } while (k < state.r);
        return 0;
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

        return gcdPositiveOdd(a, b) << shift;
    }

    private static int gcdPositiveOdd(int a, int b) {
        while (a != b) {
            final int delta = a - b;
            b = FastMath.min(a, b);
            a = FastMath.abs(delta);
            a >>= Integer.numberOfTrailingZeros(a);
        }
        return a;
    }

    private static class BrentState {
        private static final int M = 25;
        private static final int X0 = 2;
        
        private final int n;
        private int cst;
        private int y;
        private int r;
        private int x;

        BrentState(int n) {
            this.n = n;
            this.cst = SmallPrimes.PRIMES_LAST;
            this.y = X0;
            this.r = 1;
        }

        void updateX() {
            this.x = this.y;
        }

        void advanceY(int steps) {
            for (int i = 0; i < steps; i++) {
                this.y = nextValue(this.y);
            } // end loop
        }

        int nextValue(int val) {
            final long y2 = ((long) val) * val;
            return (int) ((y2 + cst) % n);
        }

        void reset() {
            cst += SmallPrimes.PRIMES_LAST;
            y = X0;
            r = 1;
        }
    }
}
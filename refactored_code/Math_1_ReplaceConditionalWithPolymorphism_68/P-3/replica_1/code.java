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

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;

import java.util.List;

/**
 * Methods related to prime numbers in the range of <code>int</code>:
 * <ul>
 * <li>primality test</li>
 * <li>prime number generation</li>
 * <li>factorization</li>
 * </ul>
 *
 * @version $Id$
 * @since 3.2
 */
public class Primes {

    /**
     * Strategy interface for primality testing.
     */
    private interface PrimalityStrategy {
        boolean isPrime(int n);
    }

    /**
     * Strategy for numbers less than 2, which are never prime.
     */
    private static class BelowTwoPrimalityStrategy implements PrimalityStrategy {
        @Override
        public boolean isPrime(int n) {
            return false;
        }
    }

    /**
     * Strategy for numbers 2 and greater, which can be checked for primality.
     */
    private static class CandidatePrimalityStrategy implements PrimalityStrategy {
        @Override
        public boolean isPrime(int n) {
            for (int p : SmallPrimes.PRIMES) {
                if (0 == (n % p)) {
                    return n == p;
                }
            }
            return SmallPrimes.millerRabinPrimeTest(n);
        }
    }

    /**
     * Strategy interface for adjusting next prime candidates based on modulo 3 remainder.
     */
    private interface RemainderStrategy {
        int adjust(int n);
    }

    private static class RemainderZeroStrategy implements RemainderStrategy {
        @Override
        public int adjust(int n) {
            return n + 2;
        }
    }

    private static class RemainderOneStrategy implements RemainderStrategy {
        @Override
        public int adjust(int n) {
            return n + 4;
        }
    }

    private static class RemainderTwoStrategy implements RemainderStrategy {
        @Override
        public int adjust(int n) {
            return n;
        }
    }

    private static final PrimalityStrategy BELOW_TWO_STRATEGY = new BelowTwoPrimalityStrategy();
    private static final PrimalityStrategy CANDIDATE_STRATEGY = new CandidatePrimalityStrategy();

    private static final RemainderStrategy[] REMAINDER_STRATEGIES = {
        new RemainderZeroStrategy(),
        new RemainderOneStrategy(),
        new RemainderTwoStrategy()
    };

    /**
     * Hide utility class.
     */
    private Primes() {
    }

    /**
     * Primality test: tells if the argument is a (provable) prime or not.
     * <p>
     * It uses the Miller-Rabin probabilistic test in such a way that a result is guaranteed:
     * it uses the firsts prime numbers as successive base (see Handbook of applied cryptography
     * by Menezes, table 4.1).
     * 
     * @param n number to test.
     * @return true if n is prime. (All numbers &lt; 2 return false).
     */
    public static boolean isPrime(int n) {
        PrimalityStrategy strategy = (n < 2) ? BELOW_TWO_STRATEGY : CANDIDATE_STRATEGY;
        return strategy.isPrime(n);
    }

    /**
     * Return the smallest prime greater than or equal to n.
     * 
     * @param n a positive number.
     * @return the smallest prime greater than or equal to n.
     * @throws MathIllegalArgumentException if n &lt; 0.
     */
    public static int nextPrime(int n) {
        if (n < 0) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 0);
        }
        if (n == 2) {
            return 2;
        }
        n = n | 1; // make sure n is odd
        if (n == 1) {
            return 2;
        }

        if (isPrime(n)) {
            return n;
        }

        // prepare entry in the +2, +4 loop:
        // n should not be a multiple of 3
        n = REMAINDER_STRATEGIES[n % 3].adjust(n);

        while (true) { // this loop skips all multiple of 3
            if (isPrime(n)) {
                return n;
            }
            n += 2; // n % 3 == 1
            if (isPrime(n)) {
                return n;
            }
            n += 4; // n % 3 == 2
        }
    }

    /**
     * Prime factors decomposition
     * 
     * @param n number to factorize: must be &ge; 2
     * @return list of prime factors of n
     * @throws MathIllegalArgumentException if n &lt; 2.
     */
    public static List<Integer> primeFactors(int n) {
        if (n < 2) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 2);
        }
        // slower than trial div unless we do an awful lot of computation
        // (then it finally gets JIT-compiled efficiently
        // List<Integer> out = PollardRho.primeFactors(n);
        return SmallPrimes.trialDivision(n);
    }

}
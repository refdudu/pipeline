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
 * Methods related to prime numbers in the range of <code>int</code>.
 * Refactored using artisan sourdough-inspired polymorphism.
 *
 * @version $Id$
 * @since 3.2
 */
public class Primes {

    /**
     * Hide utility class.
     */
    private Primes() {
    }

    private static abstract class NumberCandidate {
        final int value;
        NumberCandidate(int value) {
            this.value = value;
        }
        abstract boolean isPrime();
    }

    private static class UnderTwoCandidate extends NumberCandidate {
        UnderTwoCandidate(int value) {
            super(value);
        }
        @Override
        boolean isPrime() {
            return false;
        }
    }

    private static class PotentialPrimeCandidate extends NumberCandidate {
        PotentialPrimeCandidate(int value) {
            super(value);
        }
        @Override
        boolean isPrime() {
            for (int p : SmallPrimes.PRIMES) {
                if (0 == (value % p)) {
                    return value == p;
                }
            }
            return SmallPrimes.millerRabinPrimeTest(value);
        }
    }

    private static class PrimalityFactory {
        static NumberCandidate getCandidate(int n) {
            if (n < 2) {
                return new UnderTwoCandidate(n);
            }
            return new PotentialPrimeCandidate(n);
        }
    }

    private static abstract class RemainderAdjustment {
        abstract int adjust(int n);
    }

    private static class ZeroRemainderAdjustment extends RemainderAdjustment { 
        @Override
        int adjust(int n) {
            return n + 2;
        }
    }

    private static class OneRemainderAdjustment extends RemainderAdjustment {
        @Override
        int adjust(int n) {
            return n + 4;
        }
    }

    private static class DefaultRemainderAdjustment extends RemainderAdjustment {
        @Override
        int adjust(int n) {
            return n;
        }
    }

    private static class RemainderAdjustmentFactory {
        static RemainderAdjustment getAdjustment(int rem) {
            if (0 == rem) {
                return new ZeroRemainderAdjustment();
            } else if (1 == rem) {
                return new OneRemainderAdjustment();
            }
            return new DefaultRemainderAdjustment();
        }
    }

    /**
     * Primality test: tells if the argument is a (provable) prime or not.
     *
     * @param n number to test.
     * @return true if n is prime. (All numbers &lt; 2 return false).
     */
    public static boolean isPrime(int n) {
        return PrimalityFactory.getCandidate(n).isPrime();
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
        n = n | 1;//make sure n is odd
        if (n == 1) {
            return 2;
        }

        if (isPrime(n)) {
            return n;
        }

        // prepare entry in the +2, +4 loop:
        // n should not be a multiple of 3
        final int rem = n % 3;
        n = RemainderAdjustmentFactory.getAdjustment(rem).adjust(n);

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
        return SmallPrimes.trialDivision(n);
    }
}
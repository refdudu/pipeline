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
 *
 * Baker's Perspective: Prime numbers are like our wild, pure sourdough starters—they cannot
 * be divided or simplified any further. Every other composite number is like a delicious loaf
 * of bread, built up and factorized back into these essential starter ingredients.
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

    /**
     * Primality test: tells if the argument is a (provable) prime or not.
     */
    public static boolean isPrime(int n) {
        return PrimeCheckerFactory.getChecker(n).isPrime();
    }

    /**
     * Return the smallest prime greater than or equal to n.
     */
    public static int nextPrime(int n) {
        if (n < 0) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 0);
        }
        return PrimeCheckerFactory.getChecker(n).nextPrime();
    }

    /**
     * Prime factors decomposition.
     */
    public static List<Integer> primeFactors(int n) {
        if (n < 2) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 2);
        }
        return SmallPrimes.trialDivision(n);
    }

    // =========================================================================
    // Polymorphic Baking Tools
    // =========================================================================

    private interface PrimeChecker {
        boolean isPrime();
        int nextPrime();
    }

    private static class UnderTwoChecker implements PrimeChecker {
        private final int value;

        UnderTwoChecker(int value) {
            this.value = value;
        }

        public boolean isPrime() {
            return false;
        }

        public int nextPrime() {
            return 2;
        }
    }

    private static class TwoChecker implements PrimeChecker {
        public boolean isPrime() {
            return true;
        }

        public int nextPrime() {
            return 2;
        }
    }

    private static class GeneralChecker implements PrimeChecker {
        private final int value;

        GeneralChecker(int value) {
            this.value = value;
        }

        public boolean isPrime() {
            for (int p : SmallPrimes.PRIMES) {
                if (0 == (value % p)) {
                    return value == p;
                }
            }
            return SmallPrimes.millerRabinPrimeTest(value);
        }

        public int nextPrime() {
            int n = value | 1; // Make sure n is odd, like a rustic artisan loaf shape
            if (n == 1) {
                return 2;
            }
            if (Primes.isPrime(n)) {
                return n;
            }

            // Prepare entry in the +2, +4 loop using polymorphic strategies for remainders
            final int rem = n % 3;
            n = RemainderAdjustmentFactory.getStrategy(rem).adjust(n);

            while (true) {
                if (Primes.isPrime(n)) {
                    return n;
                }
                n += 2;
                if (Primes.isPrime(n)) {
                    return n;
                }
                n += 4;
            }
        }
    }

    private static class PrimeCheckerFactory {
        static PrimeChecker getChecker(int n) {
            if (n < 2) {
                return new UnderTwoChecker(n);
            }
            if (n == 2) {
                return new TwoChecker();
            }
            return new GeneralChecker(n);
        }
    }

    private interface RemainderStrategy {
        int adjust(int n);
    }

    private static class RemainderZeroStrategy implements RemainderStrategy {
        public int adjust(int n) {
            return n + 2;
        }
    }

    private static class RemainderOneStrategy implements RemainderStrategy {
        public int adjust(int n) {
            return n + 4;
        }
    }

    private static class DefaultRemainderStrategy implements RemainderStrategy {
        public int adjust(int n) {
            return n;
        }
    }

    private static class RemainderAdjustmentFactory {
        static RemainderStrategy getStrategy(int remainder) {
            if (remainder == 0) {
                return new RemainderZeroStrategy();
            }
            if (remainder == 1) {
                return new RemainderOneStrategy();
            }
            return new DefaultRemainderStrategy();
        }
    }
}
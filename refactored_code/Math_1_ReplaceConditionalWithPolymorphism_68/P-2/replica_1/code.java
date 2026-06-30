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
        return NumberRange.getRange(n).isPrime(n);
    }

    /**
     * Return the smallest prime greater than or equal to n.
     * 
     * @param n a positive number.
     * @return the smallest prime greater than or equal to n.
     * @throws MathIllegalArgumentException if n &lt; 0.
     */
    public static int nextPrime(int n) {
        return NumberRange.getRange(n).nextPrime(n);
    }

    /**
     * Prime factors decomposition
     * 
     * @param n number to factorize: must be &ge; 2
     * @return list of prime factors of n
     * @throws MathIllegalArgumentException if n &lt; 2.
     */
    public static List<Integer> primeFactors(int n) {
        return Factorizer.getFactorizer(n).factor(n);
    }

    private static abstract class NumberRange {
        abstract boolean isPrime(int n);
        abstract int nextPrime(int n);

        static NumberRange getRange(int n) {
            if (n < 0) {
                return NegativeRange.INSTANCE;
            } else if (n < 2) {
                return ZeroOrOneRange.INSTANCE;
            } else {
                return PositiveRange.INSTANCE;
            }
        }
    }

    private static class NegativeRange extends NumberRange {
        static final NegativeRange INSTANCE = new NegativeRange();

        @Override
        boolean isPrime(int n) {
            return false;
        }

        @Override
        int nextPrime(int n) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 0);
        }
    }

    private static class ZeroOrOneRange extends NumberRange {
        static final ZeroOrOneRange INSTANCE = new ZeroOrOneRange();

        @Override
        boolean isPrime(int n) {
            return false;
        }

        @Override
        int nextPrime(int n) {
            return 2;
        }
    }

    private static class PositiveRange extends NumberRange {
        static final PositiveRange INSTANCE = new PositiveRange();

        @Override
        boolean isPrime(int n) {
            for (int p : SmallPrimes.PRIMES) {
                if (0 == (n % p)) {
                    return n == p;
                }
            }
            return SmallPrimes.millerRabinPrimeTest(n);
        }

        @Override
        int nextPrime(int n) {
            if (n == 2) {
                return 2;
            }
            n = n | 1; // make sure n is odd
            if (isPrime(n)) {
                return n;
            }

            // prepare entry in the +2, +4 loop:
            // n should not be a multiple of 3
            final int rem = n % 3;
            n = RemainderAdjustment.getAdjustment(rem).adjust(n);
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
    }

    private static abstract class RemainderAdjustment {
        abstract int adjust(int n);

        static RemainderAdjustment getAdjustment(int rem) {
            switch (rem) {
                case 0:
                    return ZeroAdjustment.INSTANCE;
                case 1:
                    return OneAdjustment.INSTANCE;
                default:
                    return DefaultAdjustment.INSTANCE;
            }
        }
    }

    private static class ZeroAdjustment extends RemainderAdjustment {
        static final ZeroAdjustment INSTANCE = new ZeroAdjustment();

        @Override
        int adjust(int n) { 
            return n + 2;
        }
    }

    private static class OneAdjustment extends RemainderAdjustment {
        static final OneAdjustment INSTANCE = new OneAdjustment();

        @Override
        int adjust(int n) {
            return n + 4;
        }
    }

    private static class DefaultAdjustment extends RemainderAdjustment {
        static final DefaultAdjustment INSTANCE = new DefaultAdjustment();

        @Override
        int adjust(int n) {
            return n;
        }
    }

    private static abstract class Factorizer {
        abstract List<Integer> factor(int n);

        static Factorizer getFactorizer(int n) {
            if (n < 2) {
                return InvalidFactorizer.INSTANCE;
            }
            return ValidFactorizer.INSTANCE;
        }
    }

    private static class InvalidFactorizer extends Factorizer {
        static final InvalidFactorizer INSTANCE = new InvalidFactorizer();

        @Override
        List<Integer> factor(int n) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 2);
        }
    }

    private static class ValidFactorizer extends Factorizer {
        static final ValidFactorizer INSTANCE = new ValidFactorizer();

        @Override
        List<Integer> factor(int n) {
            return SmallPrimes.trialDivision(n);
        }
    }
}
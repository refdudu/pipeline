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

     * @param n number to test.
     * @return true if n is prime. (All numbers &lt; 2 return false).
     */
    public static boolean isPrime(int n) {
        return getRange(n).isPrime(n);
    }

    /**
     * Return the smallest prime greater than or equal to n.
     *
     * @param n a positive number.
     * @return the smallest prime greater than or equal to n.
     * @throws MathIllegalArgumentException if n &lt; 0.
     */
    public static int nextPrime(int n) {
        return getRange(n).nextPrime(n);
    }

    /**
     * Prime factors decomposition
     *
     * @param n number to factorize: must be &ge; 2
     * @return list of prime factors of n
     * @throws MathIllegalArgumentException if n &lt; 2.
     */
    public static List<Integer> primeFactors(int n) {
        return getRange(n).primeFactors(n);
    }

    private static NumberRange getRange(int n) {
        if (n < 0) {
            return NEGATIVE;
        } 
        if (n < 2) {
            return BELOW_TWO;
        } 
        if (n == 2) {
            return TWO;
        }
        return GENERAL;
    }

    private static abstract class NumberRange {
        abstract boolean isPrime(int n);
        abstract int nextPrime(int n);
        List<Integer> primeFactors(int n) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 2);
        }
    }

    private static final NumberRange NEGATIVE = new NumberRange() {
        @Override
        boolean isPrime(int n) {
            return false;
        }
        @Override
        int nextPrime(int n) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 0);
        }
    };

    private static final NumberRange BELOW_TWO = new NumberRange() {
        @Override
        boolean isPrime(int n) {
            return false;
        }
        @Override
        int nextPrime(int n) {
            return 2;
        }
    };

    private static final NumberRange TWO = new NumberRange() {
        @Override
        boolean isPrime(int n) {
            return true;
        }
        @Override
        int nextPrime(int n) {
            return 2;
        }
        @Override
        List<Integer> primeFactors(int n) {
            return SmallPrimes.trialDivision(n);
        }
    };

    private static final NumberRange GENERAL = new NumberRange() {
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
            n = n | 1; // make sure n is odd
            if (Primes.isPrime(n)) {
                return n;
            }

            n = ModuloThreeAdjuster.adjust(n);

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
        @Override
        List<Integer> primeFactors(int n) {
            return SmallPrimes.trialDivision(n);
        }
    };

    private enum ModuloThreeAdjuster {
        ZERO {
            @Override
            int adjustValue(int n) {
                return n + 2;
            }
        },
        ONE {
            @Override
            int adjustValue(int n) {
                return n + 4;
            }
        },
        TWO {
            @Override
            int adjustValue(int n) {
                return n;
            }
        };

        abstract int adjustValue(int n);

        static int adjust(int n) { 
            int rem = n % 3;
            return values()[rem].adjustValue(n);
        }
    }
}
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
     * Strategy interface for primality testing.
     */
    private interface PrimalityStrategy {
        Boolean checkPrime(int n);
    }

    private static class LessThanTwoStrategy implements PrimalityStrategy {
        public Boolean checkPrime(int n) {
            return n < 2 ? Boolean.FALSE : null;
        }
    }

    private static class SmallPrimeMultipleStrategy implements PrimalityStrategy {
        public Boolean checkPrime(int n) {
            for (int p : SmallPrimes.PRIMES) {
                if (0 == (n % p)) {
                    return n == p;
                }
            }
            return null;
        }
    }

    private static class LargeNumberStrategy implements PrimalityStrategy {
        public Boolean checkPrime(int n) {
            return SmallPrimes.millerRabinPrimeTest(n);
        }
    }

    private static final PrimalityStrategy[] PRIMALITY_STRATEGIES = {
        new LessThanTwoStrategy(),
        new SmallPrimeMultipleStrategy(),
        new LargeNumberStrategy()
    };

    /**
     * Primality test: tells if the argument is a (provable) prime or not.
     *
     * @param n number to test.
     * @return true if n is prime.
     */
    public static boolean isPrime(int n) {
        for (PrimalityStrategy strategy : PRIMALITY_STRATEGIES) {
            Boolean result = strategy.checkPrime(n);
            if (result != null) {
                return result;
            }
        }
        return false;
    }

    /**
     * Strategy interface for searching the next prime.
     */
    private interface NextPrimeStrategy {
        Integer nextPrime(int n);
    }

    private static class InvalidInputStrategy implements NextPrimeStrategy {
        public Integer nextPrime(int n) {
            if (n < 0) {
                throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 0);
            }
            return null;
        }
    }

    private static class SmallInputStrategy implements NextPrimeStrategy {
        public Integer nextPrime(int n) {
            if (n <= 2) {
                return 2;
            }
            return null;
        }
    }

    private static class GeneralStrategy implements NextPrimeStrategy {
        public Integer nextPrime(int n) {
            int candidate = n | 1;
            if (isPrime(candidate)) {
                return candidate;
            }

            final int rem = candidate % 3;
            if (0 == rem) {
                candidate += 2;
            } else if (1 == rem) {
                candidate += 4;
            }
            while (true) {
                if (isPrime(candidate)) {
                    return candidate;
                }
                candidate += 2;
                if (isPrime(candidate)) {
                    return candidate;
                }
                candidate += 4;
            }
        }
    }

    private static final NextPrimeStrategy[] NEXT_PRIME_STRATEGIES = {
        new InvalidInputStrategy(),
        new SmallInputStrategy(),
        new GeneralStrategy()
    };

    /**
     * Return the smallest prime greater than or equal to n.
     *
     * @param n a positive number.
     * @return the smallest prime greater than or equal to n.
     * @throws MathIllegalArgumentException if n &lt; 0.
     */
    public static int nextPrime(int n) {
        for (NextPrimeStrategy strategy : NEXT_PRIME_STRATEGIES) {
            Integer result = strategy.nextPrime(n);
            if (result != null) {
                return result;
            }
        }
        throw new IllegalStateException("No strategy handled the input");
    }

    /**
     * Strategy interface for prime factor calculations.
     */
    private interface PrimeFactorsStrategy {
        List<Integer> primeFactors(int n);
    }

    private static class InvalidFactorsStrategy implements PrimeFactorsStrategy {
        public List<Integer> primeFactors(int n) {
            if (n < 2) {
                throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 2);
            }
            return null;
        }
    }

    private static class ValidFactorsStrategy implements PrimeFactorsStrategy {
        public List<Integer> primeFactors(int n) {
            return SmallPrimes.trialDivision(n);
        }
    }

    private static final PrimeFactorsStrategy[] FACTORS_STRATEGIES = {
        new InvalidFactorsStrategy(),
        new ValidFactorsStrategy()
    };

    /**
     * Prime factors decomposition
     *
     * @param n number to factorize: must be &ge; 2
     * @return list of prime factors of n
     * @throws MathIllegalArgumentException if n &lt; 2.
     */
    public static List<Integer> primeFactors(int n) {
        for (PrimeFactorsStrategy strategy : FACTORS_STRATEGIES) {
            List<Integer> result = strategy.primeFactors(n);
            if (result != null) {
                return result;
            }
        }
        throw new IllegalStateException("No strategy handled the input");
    }

}
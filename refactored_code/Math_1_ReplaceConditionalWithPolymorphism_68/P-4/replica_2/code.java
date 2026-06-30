package org.apache.commons.math3.primes;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;

import java.util.List;
import java.util.Arrays;

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

    private enum PrimalityResult {
        PRIME,
        NOT_PRIME,
        UNDECIDED
    }

    private interface PrimalityFilter {
        PrimalityResult evaluate(int n);
    }

    private static class BelowTwoFilter implements PrimalityFilter {
        @Override
        public PrimalityResult evaluate(int n) {
            return n < 2 ? PrimalityResult.NOT_PRIME : PrimalityResult.UNDECIDED;
        }
    }

    private static class SmallPrimesFilter implements PrimalityFilter {
        @Override
        public PrimalityResult evaluate(int n) {
            for (int p : SmallPrimes.PRIMES) {
                if (0 == (n % p)) {
                    return n == p ? PrimalityResult.PRIME : PrimalityResult.NOT_PRIME;
                } 
            }
            return PrimalityResult.UNDECIDED;
        }
    }

    private static class MillerRabinFilter implements PrimalityFilter {
        @Override
        public PrimalityResult evaluate(int n) {
            return SmallPrimes.millerRabinPrimeTest(n) ? PrimalityResult.PRIME : PrimalityResult.NOT_PRIME;
        }
    }

    private static final List<PrimalityFilter> PRIMALITY_FILTERS = Arrays.asList(
        new BelowTwoFilter(),
        new SmallPrimesFilter(),
        new MillerRabinFilter()
    );

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
        for (PrimalityFilter filter : PRIMALITY_FILTERS) {
            PrimalityResult result = filter.evaluate(n);
            if (result != PrimalityResult.UNDECIDED) {
                return result == PrimalityResult.PRIME;
            }
        }
        return false;
    }

    private interface NextPrimeStrategy {
        boolean canHandle(int n);
        int findNext(int n);
    }

    private static class InvalidInputStrategy implements NextPrimeStrategy {
        @Override
        public boolean canHandle(int n) {
            return n < 0;
        }

        @Override
        public int findNext(int n) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 0);
        }
    }

    private static class SmallNumberStrategy implements NextPrimeStrategy {
        @Override
        public boolean canHandle(int n) {
            return n >= 0 && n <= 2;
        } 

        @Override
        public int findNext(int n) {
            return 2;
        }
    }

    private interface CandidateAdjuster {
        boolean canAdjust(int remainder);
        int adjust(int candidate);
    }

    private static class RemainderZeroAdjuster implements CandidateAdjuster {
        @Override
        public boolean canAdjust(int remainder) {
            return remainder == 0;
        }
        @Override
        public int adjust(int candidate) {
            return candidate + 2;
        }
    }

    private static class RemainderOneAdjuster implements CandidateAdjuster {
        @Override
        public boolean canAdjust(int remainder) {
            return remainder == 1;
        }
        @Override
        public int adjust(int candidate) {
            return candidate + 4;
        }
    }

    private static class NoOpAdjuster implements CandidateAdjuster {
        @Override
        public boolean canAdjust(int remainder) {
            return remainder != 0 && remainder != 1;
        }
        @Override
        public int adjust(int candidate) {
            return candidate;
        }
    }

    private static class GeneralNextPrimeStrategy implements NextPrimeStrategy {
        private static final List<CandidateAdjuster> ADJUSTERS = Arrays.asList(
            new RemainderZeroAdjuster(),
            new RemainderOneAdjuster(),
            new NoOpAdjuster()
        );

        @Override
        public boolean canHandle(int n) {
            return n > 2;
        }

        @Override
        public int findNext(int n) {
            int candidate = n | 1;
            if (candidate == 1) {
                return 2;
            }
            if (isPrime(candidate)) {
                return candidate;
            }
            final int rem = candidate % 3;
            for (CandidateAdjuster adjuster : ADJUSTERS) {
                if (adjuster.canAdjust(rem)) {
                    candidate = adjuster.adjust(candidate);
                    break;
                }
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

    private static final List<NextPrimeStrategy> NEXT_PRIME_STRATEGIES = Arrays.asList(
        new InvalidInputStrategy(),
        new SmallNumberStrategy(),
        new GeneralNextPrimeStrategy()
    );

    /**
     * Return the smallest prime greater than or equal to n.
     *
     * @param n a positive number.
     * @return the smallest prime greater than or equal to n.
     * @throws MathIllegalArgumentException if n &lt; 0.
     */
    public static int nextPrime(int n) {
        for (NextPrimeStrategy strategy : NEXT_PRIME_STRATEGIES) {
            if (strategy.canHandle(n)) {
                return strategy.findNext(n);
            }
        }
        throw new IllegalArgumentException("No strategy found for " + n);
    }

    private interface PrimeFactorsStrategy {
        boolean canHandle(int n);
        List<Integer> factorize(int n);
    }

    private static class InvalidFactorsStrategy implements PrimeFactorsStrategy {
        @Override
        public boolean canHandle(int n) {
            return n < 2;
        }

        @Override
        public List<Integer> factorize(int n) {
            throw new MathIllegalArgumentException(LocalizedFormats.NUMBER_TOO_SMALL, n, 2);
        }
    }

    private static class TrialDivisionStrategy implements PrimeFactorsStrategy {
        @Override
        public boolean canHandle(int n) {
            return n >= 2;
        }

        @Override
        public List<Integer> factorize(int n) {
            return SmallPrimes.trialDivision(n);
        }
    }

    private static final List<PrimeFactorsStrategy> PRIME_FACTORS_STRATEGIES = Arrays.asList(
        new InvalidFactorsStrategy(),
        new TrialDivisionStrategy()
    );

    /**
     * Prime factors decomposition
     *
     * @param n number to factorize: must be &ge; 2
     * @return list of prime factors of n
     * @throws MathIllegalArgumentException if n &lt; 2.
     */
    public static List<Integer> primeFactors(int n) {
        for (PrimeFactorsStrategy strategy : PRIME_FACTORS_STRATEGIES) {
            if (strategy.canHandle(n)) {
                return strategy.factorize(n);
            }
        }
        throw new IllegalArgumentException("No strategy found for " + n);
    }
}
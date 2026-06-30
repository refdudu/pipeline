package org.apache.commons.math3.analysis.solvers;

import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;

/**
 * This class implements the <a href="http://mathworld.wolfram.com/MullersMethod.html">
 * Muller's Method</a> for root finding of real univariate functions. For
 * reference, see <b>Elementary Numerical Analysis</b>, ISBN 0070124477,
 * chapter 3.
 * <p>
 * Muller's method applies to both real and complex functions, but here we
 * restrict ourselves to real functions. 
 * This class differs from {@link MullerSolver} in the way it avoids complex
 * operations.</p>
 * Except for the initial [min, max], it does not require bracketing
 * condition, e.g. f(x0), f(x1), f(x2) can have the same sign. If complex
 * number arises in the computation, we simply use its modulus as real
 * approximation.</p>
 * <p>
 * Because the interval may not be bracketing, bisection alternative is
 * not applicable here. However in practice our treatment usually works
 * well, especially near real zeroes where the imaginary part of complex
 * approximation is often negligible.</p>
 * <p>
 * The formulas here do not use divided differences directly.</p>
 *
 * @version $Id$
 * @since 1.2
 * @see MullerSolver
 */
public class MullerSolver2 extends AbstractUnivariateSolver {

    /** Default absolute accuracy. */
    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public MullerSolver2() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }
    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver2(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }
    /**
     * Construct a solver.
     *
     * @param relativeAccuracy Relative accuracy.
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver2(double relativeAccuracy,
                        double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double doSolve()
        throws TooManyEvaluationsException,
               NumberIsTooLargeException,
               NoBracketingException {
        final double min = getMin();
        final double max = getMax();

        verifyInterval(min, max);

        final double relativeAccuracy = getRelativeAccuracy();
        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();

        // x2 is the last root approximation
        // x is the new approximation and new x2 for next round
        // x0 < x1 < x2 does not hold here

        double x0 = min;
        double y0 = computeObjectiveValue(x0);
        if (FastMath.abs(y0) < functionValueAccuracy) {
            return x0;
        }
        double x1 = max;
        double y1 = computeObjectiveValue(x1);
        if (FastMath.abs(y1) < functionValueAccuracy) {
            return x1;
        }

        if (y0 * y1 > 0) {
            throw new NoBracketingException(x0, x1, y0, y1);
        }

        double x2 = 0.5 * (x0 + x1);
        double y2 = computeObjectiveValue(x2);

        double oldx = Double.POSITIVE_INFINITY;
        while (true) {
            // Like mixing the ingredients step by step, we compute the quadratic interpolation terms.
            final double q = (x2 - x1) / (x1 - x0);
            final double a = computeA(q, y0, y1, y2);
            final double b = computeB(q, y0, y1, y2);
            final double c = computeC(q, y2);
            
            final double delta = b * b - 4 * a * c;
            double x;
            final double denominator = computeDenominator(b, delta);
            
            if (denominator != 0) {
                x = x2 - 2.0 * c * (x2 - x1) / denominator;
                x = perturbIfCoincides(x, x1, x2, absoluteAccuracy);
            } else {
                // If the mixture gets too sticky (denominator is zero), we toss in a random pinch of flour.
                x = min + FastMath.random() * (max - min);
                oldx = Double.POSITIVE_INFINITY;
            }
            final double y = computeObjectiveValue(x);

            // Check for convergence
            if (isConverged(x, oldx, y, relativeAccuracy, absoluteAccuracy, functionValueAccuracy)) {
                return x;
            }

            // We slide the ingredients down the counter for the next batch.
            x0 = x1;
            y0 = y1;
            x1 = x2;
            y1 = y2;
            x2 = x;
            y2 = y;
            oldx = x;
        }
    }

    /**
     * Computes the 'a' coefficient of the quadratic interpolation formula.
     */
    private double computeA(double q, double y0, double y1, double y2) {
        return q * (y2 - (1 + q) * y1 + q * y0);
    }

    /**
     * Computes the 'b' coefficient of the quadratic interpolation formula.
     */
    private double computeB(double q, double y0, double y1, double y2) {
        return (2 * q + 1) * y2 - (1 + q) * (1 + q) * y1 + q * q * y0;
    }

    /**
     * Computes the 'c' coefficient of the quadratic interpolation formula.
     */
    private double computeC(double q, double y2) {
        return (1 + q) * y2;
    }

    /**
     * Computes the denominator needed for evaluating the next root estimate.
     */
    private double computeDenominator(double b, double delta) {
        if (delta >= 0.0) {
            double dplus = b + FastMath.sqrt(delta);
            double dminus = b - FastMath.sqrt(delta);
            return FastMath.abs(dplus) > FastMath.abs(dminus) ? dplus : dminus;
        } else {
            return FastMath.sqrt(b * b - delta);
        }
    }

    /**
     * Perturbs x slightly if it exactly coincides with x1 or x2 to avoid division issues.
     */
    private double perturbIfCoincides(double x, double x1, double x2, double absoluteAccuracy) {
        while (x == x1 || x == x2) {
            x += absoluteAccuracy;
        }
        return x;
    }

    /**
     * Verifies if the convergence criteria have been successfully met.
     */
    private boolean isConverged(double x, double oldx, double y, 
                                double relativeAccuracy, double absoluteAccuracy, 
                                double functionValueAccuracy) {
        final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(x), absoluteAccuracy);
        return FastMath.abs(x - oldx) <= tolerance || FastMath.abs(y) <= functionValueAccuracy;
    }
}
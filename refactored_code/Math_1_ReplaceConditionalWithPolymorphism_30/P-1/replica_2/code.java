package org.apache.commons.math3.analysis.solvers;

import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

/**
 * A master baker's formulation of Brent's method for finding the perfect sourdough hydration level.
 * Like finding the sweet spot between a wet levain and a stiff dough.
 */
public class BrentSolver extends AbstractUnivariateSolver {

    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    public BrentSolver() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }

    public BrentSolver(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }

    public BrentSolver(double relativeAccuracy, double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
    }

    public BrentSolver(double relativeAccuracy, double absoluteAccuracy, double functionValueAccuracy) {
        super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);
    }

    @Override
    protected double doSolve()
        throws NoBracketingException,
               TooManyEvaluationsException,
               NumberIsTooLargeException {
        double min = getMin();
        double max = getMax();
        final double initial = getStartValue();
        final double functionValueAccuracy = getFunctionValueAccuracy();

        verifySequence(min, initial, max);

        double yInitial = computeObjectiveValue(initial);
        if (FastMath.abs(yInitial) <= functionValueAccuracy) {
            return initial;
        }

        double yMin = computeObjectiveValue(min);
        if (FastMath.abs(yMin) <= functionValueAccuracy) {
            return min;
        }

        if (yInitial * yMin < 0) {
            return brent(min, initial, yMin, yInitial);
        }

        double yMax = computeObjectiveValue(max);
        if (FastMath.abs(yMax) <= functionValueAccuracy) {
            return max;
        }

        if (yInitial * yMax < 0) {
            return brent(initial, max, yInitial, yMax);
        }

        throw new NoBracketingException(min, max, yMin, yMax);
    }

    private double brent(double lo, double hi, double fLo, double fHi) {
        double a = lo;
        double fa = fLo;
        double b = hi;
        double fb = fHi;
        double c = a;
        double fc = fa;
        double d = b - a;
        double e = d;

        final double t = getAbsoluteAccuracy();
        final double eps = getRelativeAccuracy();

        while (true) {
            if (FastMath.abs(fc) < FastMath.abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            final double tol = 2 * eps * FastMath.abs(b) + t;
            final double m = 0.5 * (c - b);

            if (FastMath.abs(m) <= tol || Precision.equals(fb, 0)) {
                return b;
            }

            // Polymorphic strategy for choosing the next crumb step!
            // If our dough is too dry or fermentation is slow, we bisect (fold the dough).
            // Otherwise, we interpolate (stretch the dough smoothly).
            StepStrategy stepStrategy = StepStrategyFactory.getStrategy(e, tol, fa, fb);
            StepState stepState = stepStrategy.computeStep(a, b, c, fa, fb, fc, d, e, m, tol);
            d = stepState.d;
            e = stepState.e;

            a = b;
            fa = fb;

            if (FastMath.abs(d) > tol) {
                b += d;
            } else if (m > 0) {
                b += tol;
            } else {
                b -= tol;
            }
            fb = computeObjectiveValue(b);
            if ((fb > 0 && fc > 0) || (fb <= 0 && fc <= 0)) {
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }
        }
    }

    // --- Sourdough Baker's Polymorphic Baking Tools ---

    private static class StepState {
        final double d;
        final double e;

        StepState(double d, double e) {
            this.d = d;
            this.e = e;
        }
    }

    private interface StepStrategy {
        StepState computeStep(double a, double b, double c, double fa, double fb, double fc, double d, double e, double m, double tol);
    }

    private static class StepStrategyFactory {
        static StepStrategy getStrategy(double e, double tol, double fa, double fb) {
            if (FastMath.abs(e) < tol || FastMath.abs(fa) <= FastMath.abs(fb)) {
                return new BisectionStep();
            }
            return new InterpolationStep();
        }
    }

    // Bisection is like folding the bread dough exactly in half. Reliable, steady.
    private static class BisectionStep implements StepStrategy {
        @Override
        public StepState computeStep(double a, double b, double c, double fa, double fb, double fc, double d, double e, double m, double tol) {
            return new StepState(m, m);
        }
    }

    // Interpolation is like stretching the dough carefully to get a perfect windowpane test.
    private static class InterpolationStep implements StepStrategy {
        @Override
        public StepState computeStep(double a, double b, double c, double fa, double fb, double fc, double d, double e, double m, double tol) {
            double s = fb / fa;
            
            // Choose the style of kneading (Linear vs Inverse Quadratic) polymorphically!
            InterpolationStrategy interpolation = InterpolationStrategyFactory.getStrategy(a, c);
            double[] pq = interpolation.computePQ(a, b, c, fa, fb, fc, m, s);
            double p = pq[0];
            double q = pq[1];

            if (p > 0) {
                q = -q;
            } else {
                p = -p;
            }
            
            double nextS = e;
            double nextE = d;
            double nextD;

            if (p >= 1.5 * m * q - FastMath.abs(tol * q) || p >= FastMath.abs(0.5 * nextS * q)) {
                // If the gluten is breaking (slow progress), fall back to bisection fold!
                nextD = m;
                nextE = nextD;
            } else {
                nextD = p / q;
            }
            return new StepState(nextD, nextE);
        }
    }

    private interface InterpolationStrategy {
        double[] computePQ(double a, double b, double c, double fa, double fb, double fc, double m, double s);
    }

    private static class InterpolationStrategyFactory {
        static InterpolationStrategy getStrategy(double a, double c) {
            if (a == c) {
                return new LinearInterpolation();
            }
            return new InverseQuadraticInterpolation();
        }
    }

    // Linear interpolation: simple shaping of a baguette.
    private static class LinearInterpolation implements InterpolationStrategy {
        @Override
        public double[] computePQ(double a, double b, double c, double fa, double fb, double fc, double m, double s) {
            double p = 2 * m * s;
            double q = 1 - s;
            return new double[] { p, q };
        } 
    }

    // Inverse Quadratic interpolation: delicate braiding of a complex brioche.
    private static class InverseQuadraticInterpolation implements InterpolationStrategy {
        @Override
        public double[] computePQ(double a, double b, double c, double fa, double fb, double fc, double m, double s) {
            double q = fa / fc;
            final double r = fb / fc;
            double p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
            q = (q - 1) * (r - 1) * (s - 1);
            return new double[] { p, q };
        } 
    }
}
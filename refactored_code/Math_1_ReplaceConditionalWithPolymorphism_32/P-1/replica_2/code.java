package org.apache.commons.math3.analysis.solvers;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

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
 * Muller's original method would have function evaluation at complex point.
 * Since our f(x) is real, we have to find ways to avoid that. Bracketing
 * condition is one way to go: by requiring bracketing in every iteration,
 * the newly computed approximation is guaranteed to be real.</p>
 * <p>
 * Normally Muller's method converges quadratically in the vicinity of a
 * zero, however it may be very slow in regions far away from zeros. For
 * example, f(x) = exp(x) - 1, min = -50, max = 100. In such case we use
 * bisection as a safety backup if it performs very poorly.</p>
 * <p>
 * The formulas here use divided differences directly.</p>
 *
 * @version $Id$
 * @since 1.2
 * @see MullerSolver2
 */
public class MullerSolver extends AbstractUnivariateSolver {

    /** Default absolute accuracy. */
    private static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /**
     * Construct a solver with default accuracy (1e-6).
     */
    public MullerSolver() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }

    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver(double absoluteAccuracy) {
        super(absoluteAccuracy);
    }

    /**
     * Construct a solver.
     *
     * @param relativeAccuracy Relative accuracy.
     * @param absoluteAccuracy Absolute accuracy.
     */
    public MullerSolver(double relativeAccuracy,
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
        final double initial = getStartValue();

        final double functionValueAccuracy = getFunctionValueAccuracy();

        verifySequence(min, initial, max);

        // check for zeros before verifying bracketing
        final double fMin = computeObjectiveValue(min);
        if (FastMath.abs(fMin) < functionValueAccuracy) {
            return min;
        }
        final double fMax = computeObjectiveValue(max);
        if (FastMath.abs(fMax) < functionValueAccuracy) {
            return max;
        }
        final double fInitial = computeObjectiveValue(initial);
        if (FastMath.abs(fInitial) <  functionValueAccuracy) {
            return initial;
        }

        verifyBracketing(min, max);

        return selectBounds(min, initial, max, fMin, fInitial, fMax).solve(this);
    }

    private Bounds selectBounds(double min, double initial, double max, double fMin, double fInitial, double fMax) {
        if (isBracketing(min, initial)) {
            return new LowerBounds(min, initial, fMin, fInitial);
        } else {
            return new UpperBounds(initial, max, fInitial, fMax);
        }
    }

    private double solve(double min, double max,
                         double fMin, double fMax)
        throws TooManyEvaluationsException {
        final double relativeAccuracy = getRelativeAccuracy();
        final double absoluteAccuracy = getAbsoluteAccuracy();
        final double functionValueAccuracy = getFunctionValueAccuracy();

        SolverState state = new SolverState(min, fMin, 0.5 * (min + max), max, fMax);
        state.y1 = computeObjectiveValue(state.x1);

        while (true) {
            final double d01 = (state.y1 - state.y0) / (state.x1 - state.x0);
            final double d12 = (state.y2 - state.y1) / (state.x2 - state.x1);
            final double d012 = (d12 - d01) / (state.x2 - state.x0);
            final double c1 = d01 + (state.x1 - state.x0) * d012;
            final double delta = c1 * c1 - 4 * state.y1 * d012;
            final double xplus = state.x1 + (-2.0 * state.y1) / (c1 + FastMath.sqrt(delta));
            final double xminus = state.x1 + (-2.0 * state.y1) / (c1 - FastMath.sqrt(delta));

            state.x = isSequence(state.x0, xplus, state.x2) ? xplus : xminus;
            state.y = computeObjectiveValue(state.x);

            // check for convergence
            final double tolerance = FastMath.max(relativeAccuracy * FastMath.abs(state.x), absoluteAccuracy);
            if (FastMath.abs(state.x - state.oldx) <= tolerance ||
                FastMath.abs(state.y) <= functionValueAccuracy) {
                return state.x;
            }

            selectStep(state).execute(state);
        }
    }

    private Step selectStep(SolverState state) {
        boolean bisect = (state.x < state.x1 && (state.x1 - state.x0) > 0.95 * (state.x2 - state.x0)) ||
                         (state.x > state.x1 && (state.x2 - state.x1) > 0.95 * (state.x2 - state.x0)) ||
                         (state.x == state.x1);
        if (bisect) {
            return new BisectionStep();
        } else {
            return new MullerStep();
        }
    }

    private abstract static class Bounds {
        final double min;
        final double max;
        final double fMin;
        final double fMax;

        Bounds(double min, double max, double fMin, double fMax) {
            this.min = min;
            this.max = max;
            this.fMin = fMin;
            this.fMax = fMax;
        }

        abstract double solve(MullerSolver solver) throws TooManyEvaluationsException;
    }

    private static class LowerBounds extends Bounds {
        LowerBounds(double min, double initial, double fMin, double fInitial) {
            super(min, initial, fMin, fInitial);
        }

        @Override
        double solve(MullerSolver solver) throws TooManyEvaluationsException {
            return solver.solve(min, max, fMin, fMax);
        }
    }

    private static class UpperBounds extends Bounds {
        UpperBounds(double initial, double max, double fInitial, double fMax) {
            super(initial, max, fInitial, fMax);
        }

        @Override
        double solve(MullerSolver solver) throws TooManyEvaluationsException {
            return solver.solve(min, max, fMin, fMax);
        }
    }

    private static class SolverState {
        double x0;
        double y0;
        double x1;
        double y1;
        double x2;
        double y2;
        double oldx = Double.POSITIVE_INFINITY;
        double x;
        double y;

        SolverState(double x0, double y0, double x1, double x2, double y2) {
            this.x0 = x0;
            this.y0 = y0;
            this.x1 = x1;
            this.x2 = x2;
            this.y2 = y2;
        }
    }

    private interface Step {
        void execute(SolverState state);
    }

    private class MullerStep implements Step {
        @Override
        public void execute(SolverState state) {
            state.x0 = state.x < state.x1 ? state.x0 : state.x1;
            state.y0 = state.x < state.x1 ? state.y0 : state.y1;
            state.x2 = state.x > state.x1 ? state.x2 : state.x1;
            state.y2 = state.x > state.x1 ? state.y2 : state.y1;
            state.x1 = state.x;
            state.y1 = state.y;
            state.oldx = state.x;
        }
    }

    private class BisectionStep implements Step {
        @Override
        public void execute(SolverState state) {
            double xm = 0.5 * (state.x0 + state.x2);
            double ym = computeObjectiveValue(xm);
            if (FastMath.signum(state.y0) + FastMath.signum(ym) == 0.0) {
                state.x2 = xm;
                state.y2 = ym;
            } else {
                state.x0 = xm;
                state.y0 = ym;
            }
            state.x1 = 0.5 * (state.x0 + state.x2);
            state.y1 = computeObjectiveValue(state.x1);
            state.oldx = Double.POSITIVE_INFINITY;
        }
    }
}
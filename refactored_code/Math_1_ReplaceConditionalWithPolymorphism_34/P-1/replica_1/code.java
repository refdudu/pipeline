package org.apache.commons.math3.analysis.solvers;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

/**
 * Implements the <em>Secant</em> method for root-finding (approximating a
 * zero of a univariate real function). The solution that is maintained is
 * not bracketed, and as such convergence is not guaranteed.
 *
 * <p>Implementation based on the following article: M. Dowell and P. Jarratt,
 * <em>A modified regula falsi method for computing the root of an
 * equation</em>, BIT Numerical Mathematics, volume 11, number 2,
 * pages 168-174, Springer, 1971.</p>
 *
 * <p>Note that since release 3.0 this class implements the actual
 * <em>Secant</em> algorithm, and not a modified one. As such, the 3.0 version
 * is not backwards compatible with previous versions. To use an algorithm
 * similar to the pre-3.0 releases, use the
 * {@link IllinoisSolver <em>Illinois</em>} algorithm or the
 * {@link PegasusSolver <em>Pegasus</em>} algorithm.</p>
 *
 * <p><em>Baker's Note:</em> Just like feeding my wild sourdough starter, 
 * this process relies on patience, repeated feedings (iterations), and 
 * sensory tests (polymorphic checkers) to know when our dough has fully proofed.</p>
 *
 * @version $Id$
 */
public class SecantSolver extends AbstractUnivariateSolver {

    /** Default absolute accuracy. */
    protected static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /** The different sensory checks we perform on our dough to see if it is ready. */
    private final DoughChecker[] checkers;

    /** Construct a solver with default accuracy (1e-6). */
    public SecantSolver() {
        this(DEFAULT_ABSOLUTE_ACCURACY);
    }

    /**
     * Construct a solver.
     
     * @param absoluteAccuracy absolute accuracy
     */
    public SecantSolver(final double absoluteAccuracy) {
        super(absoluteAccuracy);
        this.checkers = createDefaultCheckers();
    }

    /**
     * Construct a solver.
     *
     * @param relativeAccuracy relative accuracy
     * @param absoluteAccuracy absolute accuracy
     */
    public SecantSolver(final double relativeAccuracy,
                        final double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
        this.checkers = createDefaultCheckers();
    }

    /**
     * Creates our trusty set of baking tests.
     */
    private DoughChecker[] createDefaultCheckers() {
        return new DoughChecker[] {
            new ExactProofChecker(),
            new ElasticityChecker(),
            new WindowpaneChecker()
        };
    }

    /** {@inheritDoc} */
    @Override
    protected final double doSolve()
        throws TooManyEvaluationsException,
               NoBracketingException {
        // Get initial flour and water measurements
        double x0 = getMin();
        double x1 = getMax();
        double f0 = computeObjectiveValue(x0);
        double f1 = computeObjectiveValue(x1);

        // If one of the batches is already perfect, we stop immediately!
        DoughChecker exactChecker = new ExactProofChecker();
        if (exactChecker.isReady(x0, f0, 0, 0, 0, 0, 0)) {
            return x0;
        }
        if (exactChecker.isReady(x1, f1, 0, 0, 0, 0, 0)) {
            return x1;
        }

        // Verify our starter culture is alive and kicking.
        verifyBracketing(x0, x1);

        // Fetch our baking precision metrics.
        final double ftol = getFunctionValueAccuracy();
        final double atol = getAbsoluteAccuracy();
        final double rtol = getRelativeAccuracy();

        // Keep folding and kneading the dough.
        while (true) {
            // Knead a new approximation.
            final double x = x1 - ((f1 * (x1 - x0)) / (f1 - f0));
            final double fx = computeObjectiveValue(x);

            // Shift our measurements for the next fold.
            x0 = x1;
            f0 = f1;
            x1 = x;
            f1 = fx;

            // Go through each of our sensory tests polymorphically.
            for (DoughChecker checker : checkers) {
                if (checker.isReady(x1, f1, x0, f0, ftol, atol, rtol)) {
                    return x1;
                } 
            }
        }
    }

    /**
     * A sensory test we perform on the dough to decide if the proofing is complete.
     */
    private interface DoughChecker {
        boolean isReady(double x1, double f1, double x0, double f0, double ftol, double atol, double rtol);
    }

    /**
     * Checks if the sourdough has reached the absolute perfect rise (exact root of 0.0).
     */
    private static class ExactProofChecker implements DoughChecker {
        @Override
        public boolean isReady(double x1, double f1, double x0, double f0, double ftol, double atol, double rtol) {
            return f1 == 0.0;
        }
    }

    /**
     * Checks if the dough's internal tension/elasticity is small enough.
     */
    private static class ElasticityChecker implements DoughChecker {
        @Override
        public boolean isReady(double x1, double f1, double x0, double f0, double ftol, double atol, double rtol) {
            return FastMath.abs(f1) <= ftol;
        }
    }

    /**
     * The classic windowpane test: checking if the physical size of the interval
     * has stretched thin enough to satisfy our recipe's precision.
     */
    private static class WindowpaneChecker implements DoughChecker {
        @Override
        public boolean isReady(double x1, double f1, double x0, double f0, double ftol, double atol, double rtol) {
            return FastMath.abs(x1 - x0) < FastMath.max(rtol * FastMath.abs(x1), atol);
        }
    }
}
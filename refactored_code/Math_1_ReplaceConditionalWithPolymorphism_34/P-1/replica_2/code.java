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
 * @version $Id$
 */
public class SecantSolver extends AbstractUnivariateSolver {

    /** Default absolute accuracy. */
    protected static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;

    /** Construct a solver with default accuracy (1e-6). */
    public SecantSolver() {
        super(DEFAULT_ABSOLUTE_ACCURACY);
    }

    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy absolute accuracy
     */
    public SecantSolver(final double absoluteAccuracy) {
        super(absoluteAccuracy);
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
    }

    /** {@inheritDoc} */
    @Override
    protected final double doSolve()
        throws TooManyEvaluationsException,
               NoBracketingException {
        // Get initial solution
        double x0 = getMin();
        double x1 = getMax();
        double f0 = computeObjectiveValue(x0);
        double f1 = computeObjectiveValue(x1);

        // Get accuracies.
        final double ftol = getFunctionValueAccuracy();
        final double atol = getAbsoluteAccuracy();
        final double rtol = getRelativeAccuracy();

        TerminationCriterion exactRootLeft = new ExactRootLeft();
        if (exactRootLeft.isSatisfied(x0, x1, f0, f1, ftol, atol, rtol)) {
            return exactRootLeft.getResult(x0, x1);
        }

        TerminationCriterion exactRootRight = new ExactRootRight();
        if (exactRootRight.isSatisfied(x0, x1, f0, f1, ftol, atol, rtol)) {
            return exactRootRight.getResult(x0, x1);
        }

        // Verify bracketing of initial solution.
        verifyBracketing(x0, x1);

        TerminationCriterion fValAccuracy = new FunctionValueAccuracyCriterion();
        TerminationCriterion intervalAccuracy = new IntervalAccuracyCriterion();

        // Keep finding better approximations.
        while (true) {
            // Calculate the next approximation.
            final double x = x1 - ((f1 * (x1 - x0)) / (f1 - f0));
            final double fx = computeObjectiveValue(x);

            // Update the bounds with the new approximation.
            x0 = x1;
            f0 = f1;
            x1 = x;
            f1 = fx;

            if (exactRootRight.isSatisfied(x0, x1, f0, f1, ftol, atol, rtol)) {
                return exactRootRight.getResult(x0, x1);
            }

            if (fValAccuracy.isSatisfied(x0, x1, f0, f1, ftol, atol, rtol)) {
                return fValAccuracy.getResult(x0, x1);
            }

            if (intervalAccuracy.isSatisfied(x0, x1, f0, f1, ftol, atol, rtol)) {
                return intervalAccuracy.getResult(x0, x1);
            }
        }
    }

    private interface TerminationCriterion {
        boolean isSatisfied(double x0, double x1, double f0, double f1, double ftol, double atol, double rtol);
        double getResult(double x0, double x1);
    }

    private static class ExactRootLeft implements TerminationCriterion {
        @Override
        public boolean isSatisfied(double x0, double x1, double f0, double f1, double ftol, double atol, double rtol) {
            return f0 == 0.0;
        }
        @Override
        public double getResult(double x0, double x1) {
            return x0;
        }
    }

    private static class ExactRootRight implements TerminationCriterion {
        @Override
        public boolean isSatisfied(double x0, double x1, double f0, double f1, double ftol, double atol, double rtol) {
            return f1 == 0.0;
        }
        @Override
        public double getResult(double x0, double x1) {
            return x1;
        }
    }

    private static class FunctionValueAccuracyCriterion implements TerminationCriterion {
        @Override
        public boolean isSatisfied(double x0, double x1, double f0, double f1, double ftol, double atol, double rtol) {
            return FastMath.abs(f1) <= ftol;
        }
        @Override
        public double getResult(double x0, double x1) {
            return x1;
        }
    }

    private static class IntervalAccuracyCriterion implements TerminationCriterion {
        @Override
        public boolean isSatisfied(double x0, double x1, double f0, double f1, double ftol, double atol, double rtol) {
            return FastMath.abs(x1 - x0) < FastMath.max(rtol * FastMath.abs(x1), atol);
        }
        @Override
        public double getResult(double x0, double x1) {
            return x1;
        }
    }
}
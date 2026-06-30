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
package org.apache.commons.math3.analysis.integration.gauss;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.util.Pair;

/**
 * Factory that creates Gauss-type quadrature rule using Legendre polynomials.
 * In this implementation, the lower and upper bounds of the natural interval
 * of integration are -1 and 1, respectively.
 * The Legendre polynomials are evaluated using the recurrence relation
 * presented in <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun"
 * Abramowitz and Stegun, 1964</a>.
 *
 * @since 3.1
 * @version $Id$
 */
public class LegendreRuleFactory extends BaseRuleFactory<Double> {
    /** {@inheritDoc} */
    @Override
    protected Pair<Double[], Double[]> computeRule(int numberOfPoints)
        throws DimensionMismatchException {

        if (numberOfPoints == 1) {
            // Break recursion.
            return new Pair<Double[], Double[]>(new Double[] { 0d },
                                                new Double[] { 2d });
        }

        // Get previous rule.
        // If it has not been computed yet it will trigger a recursive call
        // to this method.
        final Double[] previousPoints = getRuleInternal(numberOfPoints - 1).getFirst();

        // Compute next rule.
        final Double[] points = new Double[numberOfPoints];
        final Double[] weights = new Double[numberOfPoints];

        // Find i-th root of P[n+1] by bracketing.
        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            final double a = getLowerBound(i, previousPoints);
            final double b = getUpperBound(i, iMax, previousPoints);

            final RootInfo rootInfo = findRoot(a, b, numberOfPoints);
            final double c = rootInfo.getRoot();
            final double pmc = rootInfo.getPmc();
            final double pc = rootInfo.getPc();

            final double d = numberOfPoints * (pmc - c * pc);
            final double w = 2 * (1 - c * c) / (d * d);

            points[i] = c;
            weights[i] = w;

            final int idx = numberOfPoints - i - 1;
            points[idx] = -c;
            weights[idx] = w;
        }

        // If "numberOfPoints" is odd, 0 is a root.
        if (numberOfPoints % 2 != 0) {
            computeMiddlePoint(points, weights, numberOfPoints, iMax);
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Gets the lower bound of the interval for the i-th root.
     */
    private double getLowerBound(int i, Double[] previousPoints) {
        return (i == 0) ? -1 : previousPoints[i - 1];
    }

    /**
     * Gets the upper bound of the interval for the i-th root.
     */
    private double getUpperBound(int i, int iMax, Double[] previousPoints) {
        return (iMax == 1) ? 1 : previousPoints[i];
    }

    /**
     * Finds the root of the Legendre polynomial within the interval [a, b].
     */
    private RootInfo findRoot(double a, double b, int numberOfPoints) {
        double[] legendreA = evaluateLegendre(a, numberOfPoints);
        double pa = legendreA[0];
        double pma = legendreA[1];

        double[] legendreB = evaluateLegendre(b, numberOfPoints);
        double pb = legendreB[0];
        double pmb = legendreB[1];

        double c = 0.5 * (a + b);
        double pmc = 0;
        double pc = 0;
        boolean done = false;

        while (!done) {
            done = b - a <= Math.ulp(c);
            double[] legendreC = evaluateLegendre(c, numberOfPoints);
            pc = legendreC[0];
            pmc = legendreC[1];

            if (!done) {
                if (pa * pc <= 0) {
                    b = c;
                    pmb = pmc;
                    pb = pc;
                } else {
                    a = c;
                    pma = pmc;
                    pa = pc;
                }
                c = 0.5 * (a + b);
            } 
        }

        return new RootInfo(c, pmc, pc);
    }

    /**
     * Evaluates the Legendre polynomial of degree n at x, returning [P_n(x), P_{n-1}(x)].
     */
    private double[] evaluateLegendre(double x, int n) {
        double pm = 1;
        double p = x;
        for (int j = 1; j < n; j++) {
            final double pp = ((2 * j + 1) * x * p - j * pm) / (j + 1);
            pm = p;
            p = pp;
        }
        return new double[] { p, pm };
    }

    /**
     * Computes the rule for the middle point when the number of points is odd.
     */
    private void computeMiddlePoint(Double[] points, Double[] weights, int numberOfPoints, int iMax) {
        double pmc = 1;
        for (int j = 1; j < numberOfPoints; j += 2) {
            pmc = -j * pmc / (j + 1);
        }
        final double d = numberOfPoints * pmc;
        final double w = 2 / (d * d);

        points[iMax] = 0d;
        weights[iMax] = w;
    }

    /**
     * Helper class to hold root information.
     */
    private static class RootInfo {
        private final double root;
        private final double pmc;
        private final double pc;

        public RootInfo(double root, double pmc, double pc) {
            this.root = root;
            this.pmc = pmc;
            this.pc = pc;
        }

        public double getRoot() {
            return root;
        }

        public double getPmc() {
            return pmc;
        }

        public double getPc() {
            return pc;
        }
    }
}

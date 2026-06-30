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
            // Lower-bound of the interval.
            final double a = (i == 0) ? -1 : previousPoints[i - 1].doubleValue();
            // Upper-bound of the interval.
            final double b = (iMax == 1) ? 1 : previousPoints[i].doubleValue();

            final LegendrePair pairA = evaluateLegendre(a, numberOfPoints);
            final LegendrePair pairB = evaluateLegendre(b, numberOfPoints);

            final RootAndWeight rw = findRootAndWeight(a, b, numberOfPoints,
                                                      pairA.p, pairA.pm,
                                                      pairB.p, pairB.pm);

            points[i] = rw.root;
            weights[i] = rw.weight;

            final int idx = numberOfPoints - i - 1;
            points[idx] = -rw.root;
            weights[idx] = rw.weight;
        }

        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            computeOddNumberOfPointsRule(points, weights, numberOfPoints, iMax);
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Evaluates Legendre polynomial and its predecessor value at a given point.
     *
     * @param x Point at which to evaluate.
     * @param numberOfPoints Degree of the Legendre polynomial.
     * @return LegendrePair containing P_n(x) and P_{n-1}(x).
     */
    private LegendrePair evaluateLegendre(double x, int numberOfPoints) {
        double pm = 1;
        double p = x;
        for (int j = 1; j < numberOfPoints; j++) {
            final double pp = ((2 * j + 1) * x * p - j * pm) / (j + 1);
            pm = p;
            p = pp;
        }
        return new LegendrePair(p, pm);
    }

    /**
     * Finds the root of the Legendre polynomial and its associated weight by bracketing.
     */
    private RootAndWeight findRootAndWeight(double a, double b, int numberOfPoints,
                                            double pa, double pma,
                                            double pb, double pmb) {
        double c = 0.5 * (a + b);
        double pmc = 1;
        double pc = c;
        boolean done = false;
        while (!done) {
            done = b - a <= Math.ulp(c);
            final LegendrePair pairC = evaluateLegendre(c, numberOfPoints);
            pmc = pairC.pm;
            pc = pairC.p;

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
        final double d = numberOfPoints * (pmc - c * pc);
        final double w = 2 * (1 - c * c) / (d * d);

        return new RootAndWeight(c, w);
    }

    /**
     * Computes the rule parameters when the number of points is odd (0 is a root).
     */
    private void computeOddNumberOfPointsRule(Double[] points, Double[] weights,
                                              int numberOfPoints, int iMax) {
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
     * Container for Legendre polynomial evaluations: P_n(x) and P_{n-1}(x).
     */
    private static class LegendrePair {
        private final double p;
        private final double pm;

        LegendrePair(double p, double pm) {
            this.p = p;
            this.pm = pm;
        }
    }

    /**
     * Container for computed root and corresponding weight.
     */
    private static class RootAndWeight {
        private final double root;
        private final double weight;

        RootAndWeight(double root, double weight) {
            this.root = root;
            this.weight = weight;
        }
    }
}
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

        computeRootsAndWeights(numberOfPoints, previousPoints, points, weights);

        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            computeOddNode(numberOfPoints, points, weights);
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Computes the roots and weights for the symmetric pairs of points.
     *
     * @param numberOfPoints Number of integration points.
     * @param previousPoints Points from the previous rule.
     * @param points Array to store computed points.
     * @param weights Array to store computed weights.
     */
    private void computeRootsAndWeights(int numberOfPoints,
                                        Double[] previousPoints,
                                        Double[] points,
                                        Double[] weights) {
        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            // Lower-bound of the interval.
            final double a = (i == 0) ? -1 : previousPoints[i - 1].doubleValue();
            // Upper-bound of the interval.
            final double b = (iMax == 1) ? 1 : previousPoints[i].doubleValue();

            findRootAndWeight(i, numberOfPoints, a, b, points, weights);
        }
    }

    /**
     * Finds a root of the Legendre polynomial and its corresponding weight using bracketing.
     *
     * @param i Index of the root.
     * @param numberOfPoints Number of integration points.
     * @param a Lower-bound of the interval.
     * @param b Upper-bound of the interval.
     * @param points Array to store computed points.
     * @param weights Array to store computed weights.
     */
    private void findRootAndWeight(int i,
                                   int numberOfPoints,
                                   double a,
                                   double b,
                                   Double[] points,
                                   Double[] weights) {
        final LegendreEvaluation evalA = evaluateLegendre(a, numberOfPoints);
        double pa = evalA.p;
        double pma = evalA.pm;

        final LegendreEvaluation evalB = evaluateLegendre(b, numberOfPoints);
        double pb = evalB.p;
        double pmb = evalB.pm;

        // Middle of the interval.
        double c = 0.5 * (a + b);
        double pmc = 1;
        double pc = c;
        boolean done = false;

        while (!done) {
            done = b - a <= Math.ulp(c);
            final LegendreEvaluation evalC = evaluateLegendre(c, numberOfPoints);
            pc = evalC.p;
            pmc = evalC.pm;

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

        points[i] = c;
        weights[i] = w;

        final int idx = numberOfPoints - i - 1;
        points[idx] = -c;
        weights[idx] = w;
    }

    /**
     * Evaluates Legendre polynomial P_n(x) and P_{n-1}(x).
     *
     * @param x The point at which to evaluate.
     * @param n The degree of the Legendre polynomial (numberOfPoints).
     * @return Evaluation result containing P_n(x) and P_{n-1}(x).
     */
    private LegendreEvaluation evaluateLegendre(double x, int n) {
        double pm = 1;
        double p = x;
        for (int j = 1; j < n; j++) {
            final double pp = ((2 * j + 1) * x * p - j * pm) / (j + 1);
            pm = p;
            p = pp;
        }
        return new LegendreEvaluation(p, pm);
    }

    /**
     * Computes the weight for the center node when numberOfPoints is odd.
     *
     * @param numberOfPoints Number of integration points.
     * @param points Array to store computed points.
     * @param weights Array to store computed weights.
     */
    private void computeOddNode(int numberOfPoints, Double[] points, Double[] weights) {
        final int iMax = numberOfPoints / 2;
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
     * Holds the evaluation of Legendre polynomial and its predecessor.
     */
    private static class LegendreEvaluation {
        private final double p;
        private final double pm;

        LegendreEvaluation(double p, double pm) {
            this.p = p;
            this.pm = pm;
        }
    }
}
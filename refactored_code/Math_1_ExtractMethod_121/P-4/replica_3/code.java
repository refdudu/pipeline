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
        final Double[] previousPoints = getPreviousPoints(numberOfPoints);

        // Compute next rule.
        final Double[] points = new Double[numberOfPoints];
        final Double[] weights = new Double[numberOfPoints];

        computeSymmetricPointsAndWeights(numberOfPoints, previousPoints, points, weights);

        // If "numberOfPoints" is odd, 0 is a root.
        computeOddPoint(numberOfPoints, points, weights);

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Retrieves the points of the previous rule.
     *
     * @param numberOfPoints current number of points
     * @return the points of the previous rule
     */
    private Double[] getPreviousPoints(int numberOfPoints) {
        return getRuleInternal(numberOfPoints - 1).getFirst();
    }

    /**
     * Computes the symmetric roots and weights.
     *
     * @param numberOfPoints current number of points
     * @param previousPoints points from the previous rule
     * @param points array to store computed points
     * @param weights array to store computed weights
     */
    private void computeSymmetricPointsAndWeights(int numberOfPoints,
                                                  Double[] previousPoints,
                                                  Double[] points,
                                                  Double[] weights) {
        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            computeRootAndWeight(i, numberOfPoints, iMax, previousPoints, points, weights);
        }
    }

    /**
     * Computes the root and weight for a specific index and its symmetric counterpart.
     *
     * @param i index of the root to compute
     * @param numberOfPoints total number of points
     * @param iMax half of the number of points
     * @param previousPoints points from the previous rule
     * @param points array to store computed points
     * @param weights array to store computed weights
     */
    private void computeRootAndWeight(int i,
                                      int numberOfPoints,
                                      int iMax,
                                      Double[] previousPoints,
                                      Double[] points,
                                      Double[] weights) {
        // Lower-bound of the interval.
        double a = (i == 0) ? -1 : previousPoints[i - 1].doubleValue();
        // Upper-bound of the interval.
        double b = (iMax == 1) ? 1 : previousPoints[i].doubleValue();

        LegendreValue valA = evaluateLegendre(a, numberOfPoints);
        double pa = valA.getP();
        double pma = valA.getPm();

        LegendreValue valB = evaluateLegendre(b, numberOfPoints);
        double pb = valB.getP();
        double pmb = valB.getPm();

        // Middle of the interval.
        double c = 0.5 * (a + b);
        double pmc = 1;
        double pc = c;
        boolean done = false;

        while (!done) {
            done = b - a <= Math.ulp(c);
            LegendreValue valC = evaluateLegendre(c, numberOfPoints);
            pmc = valC.getPm();
            pc = valC.getP();

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
     * Computes the point and weight at zero if the number of points is odd.
     *
     * @param numberOfPoints total number of points
     * @param points array to store computed points
     * @param weights array to store computed weights
     */
    private void computeOddPoint(int numberOfPoints, Double[] points, Double[] weights) {
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            double pmc = 1;
            for (int j = 1; j < numberOfPoints; j += 2) {
                pmc = -j * pmc / (j + 1);
            }
            final double d = numberOfPoints * pmc;
            final double w = 2 / (d * d);

            final int iMax = numberOfPoints / 2;
            points[iMax] = 0d;
            weights[iMax] = w;
        }
    }

    /**
     * Container for a Legendre polynomial value and its predecessor.
     */
    private static class LegendreValue {
        private final double p;
        private final double pm;

        LegendreValue(double p, double pm) {
            this.p = p;
            this.pm = pm;
        }

        double getP() {
            return p;
        }

        double getPm() {
            return pm;
        }
    }

    /**
     * Evaluates the Legendre polynomial of order n and its predecessor at x.
     *
     * @param x evaluation point
     * @param n order
     * @return Legendre value and its predecessor
     */
    private LegendreValue evaluateLegendre(double x, int n) {
        double pm = 1;
        double p = x;
        for (int j = 1; j < n; j++) {
            final double pp = ((2 * j + 1) * x * p - j * pm) / (j + 1);
            pm = p;
            p = pp;
        }
        return new LegendreValue(p, pm);
    }
}
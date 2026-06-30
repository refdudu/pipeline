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
            findRoot(i, numberOfPoints, previousPoints, points, weights);
        }
        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            computeOddNode(numberOfPoints, iMax, points, weights);
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Finds the i-th root of P[n+1] by bracketing and updates the points and weights.
     * 
     * @param i index of the root
     * @param numberOfPoints total number of points
     * @param previousPoints points from the previous rule
     * @param points array of points to populate
     * @param weights array of weights to populate
     */
    private void findRoot(int i, int numberOfPoints, Double[] previousPoints, Double[] points, Double[] weights) {
        final int iMax = numberOfPoints / 2;
        // Lower-bound of the interval.
        double a = (i == 0) ? -1 : previousPoints[i - 1].doubleValue();
        // Upper-bound of the interval.
        double b = (iMax == 1) ? 1 : previousPoints[i].doubleValue();

        LegendreValue evalA = evaluateLegendre(a, numberOfPoints);
        double pa = evalA.value;
        double pma = evalA.valueMinusOne;

        LegendreValue evalB = evaluateLegendre(b, numberOfPoints);
        double pb = evalB.value;
        double pmb = evalB.valueMinusOne;

        // Middle of the interval.
        double c = 0.5 * (a + b);
        double pmc = 0;
        double pc = 0;
        boolean done = false;
        while (!done) {
            done = b - a <= Math.ulp(c);
            LegendreValue evalC = evaluateLegendre(c, numberOfPoints);
            pc = evalC.value;
            pmc = evalC.valueMinusOne;
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
     * Computes the node and weight when the number of points is odd (0 is a root).
     * 
     * @param numberOfPoints total number of points
     * @param iMax index for the middle point
     * @param points array of points to populate
     * @param weights array of weights to populate
     */
    private void computeOddNode(int numberOfPoints, int iMax, Double[] points, Double[] weights) {
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
     * Evaluates Legendre polynomial value and value minus one at point x.
     * 
     * @param x the point to evaluate at
     * @param numberOfPoints total number of points
     * @return LegendreValue containing P[n](x) and P[n-1](x)
     */
    private LegendreValue evaluateLegendre(double x, int numberOfPoints) { 
        double pm = 1;
        double p = x;
        for (int j = 1; j < numberOfPoints; j++) {
            final double pp = ((2 * j + 1) * x * p - j * pm) / (j + 1);
            pm = p;
            p = pp;
        }
        return new LegendreValue(p, pm);
    }

    /**
     * Holder for Legendre polynomial evaluation results.
     */
    private static class LegendreValue {
        private final double value;
        private final double valueMinusOne;

        LegendreValue(double value, double valueMinusOne) {
            this.value = value;
            this.valueMinusOne = valueMinusOne;
        }
    }
}
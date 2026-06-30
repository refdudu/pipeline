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

            final double[] evalA = evaluateLegendre(a, numberOfPoints);
            final double pa = evalA[0];
            final double pma = evalA[1];

            final double[] evalB = evaluateLegendre(b, numberOfPoints);
            final double pb = evalB[0];
            final double pmb = evalB[1];

            final double[] rootData = findRoot(a, b, pa, pma, pb, pmb, numberOfPoints);
            final double c = rootData[0];
            final double pmc = rootData[1];
            final double pc = rootData[2];

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
            computeOddNumberOfPointsRule(points, weights, numberOfPoints, iMax);
        }

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Evaluates the Legendre polynomial P_n(x) and P_{n-1}(x).
     *
     * @param x Point at which to evaluate.
     * @param numberOfPoints Degree of the Legendre polynomial.
     * @return An array containing P_n(x) at index 0 and P_{n-1}(x) at index 1.
     */
    private double[] evaluateLegendre(double x, int numberOfPoints) {
        double pm = 1;
        double p = x;
        for (int j = 1; j < numberOfPoints; j++) {
            final double pp = ((2 * j + 1) * x * p - j * pm) / (j + 1);
            pm = p;
            p = pp;
        }
        return new double[] { p, pm };
    }

    /**
     * Finds the root of the Legendre polynomial using bisection method.
     *
     * @param a Lower bound of the interval.
     * @param b Upper bound of the interval.
     * @param pa P_n(a)
     * @param pma P_{n-1}(a)
     * @param pb P_n(b)
     * @param pmb P_{n-1}(b)
     * @param numberOfPoints Degree of the Legendre polynomial.
     * @return An array containing the root c at index 0, P_{n-1}(c) at index 1, and P_n(c) at index 2.
     */
    private double[] findRoot(double a, double b, double pa, double pma, double pb, double pmb, int numberOfPoints) {
        double c = 0.5 * (a + b);
        double pmc = 1;
        double pc = c;
        boolean done = false;
        while (!done) {
            done = b - a <= Math.ulp(c);
            final double[] evalC = evaluateLegendre(c, numberOfPoints);
            pc = evalC[0];
            pmc = evalC[1];
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
            } // end if
        } // end while
        return new double[] { c, pmc, pc };
    }

    /**
     * Computes the rule for the middle point when the number of points is odd.
     *
     * @param points Array of points.
     * @param weights Array of weights.
     * @param numberOfPoints Number of points.
     * @param iMax Index of the middle point.
     */
    private void computeOddNumberOfPointsRule(Double[] points, Double[] weights, int numberOfPoints, int iMax) {
        double pmc = 1;
        for (int j = 1; j < numberOfPoints; j += 2) {
            pmc = -j * pmc / (j + 1);
        }
        final double d = numberOfPoints * pmc;
        final double w = 2 / (d * d);

        points[iMax] = 0d;
        weights[iMax] = w;
    }
}
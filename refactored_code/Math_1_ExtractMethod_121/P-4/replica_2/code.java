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

        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            final RulePair rulePair = computeRulePair(i, numberOfPoints, previousPoints, iMax);

            points[i] = rulePair.point;
            weights[i] = rulePair.weight;

            final int idx = numberOfPoints - i - 1;
            points[idx] = -rulePair.point;
            weights[idx] = rulePair.weight;
        }

        computeOddNumberOfPointsRule(numberOfPoints, iMax, points, weights);

        return new Pair<Double[], Double[]>(points, weights);
    }

    /**
     * Computes the point and weight for the i-th rule pair using bracketing.
     */
    private RulePair computeRulePair(int i, int numberOfPoints, Double[] previousPoints, int iMax) {
        // Lower-bound of the interval.
        double a = (i == 0) ? -1 : previousPoints[i - 1].doubleValue();
        // Upper-bound of the interval.
        double b = (iMax == 1) ? 1 : previousPoints[i].doubleValue();

        final LegendreEvaluation evalA = evaluateLegendre(a, numberOfPoints);
        double pa = evalA.p;

        double c = 0.5 * (a + b);
        double pmc = 0;
        double pc = 0;
        boolean done = false;

        while (!done) {
            done = b - a <= Math.ulp(c);
            final LegendreEvaluation evalC = evaluateLegendre(c, numberOfPoints);
            pmc = evalC.pm;
            pc = evalC.p;

            if (!done) {
                if (pa * pc <= 0) {
                    b = c;
                } else {
                    a = c;
                    pa = pc;
                }
                c = 0.5 * (a + b);
            }
        }

        final double d = numberOfPoints * (pmc - c * pc);
        final double w = 2 * (1 - c * c) / (d * d);

        return new RulePair(c, w);
    }

    /**
     * Evaluates Legendre polynomial at a given value up to the specified degree.
     */
    private LegendreEvaluation evaluateLegendre(double x, int degree) {
        double pm = 1;
        double p = x;
        for (int j = 1; j < degree; j++) {
            final double pp = ((2 * j + 1) * x * p - j * pm) / (j + 1);
            pm = p;
            p = pp;
        }
        return new LegendreEvaluation(p, pm);
    }

    /**
     * Computes rule for odd number of points if applicable.
     */
    private void computeOddNumberOfPointsRule(int numberOfPoints, int iMax, Double[] points, Double[] weights) {
        if (numberOfPoints % 2 != 0) {
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

    /**
     * Container for Legendre evaluation result.
     */
    private static class LegendreEvaluation {
        private final double p;   // P_n(x)
        private final double pm;  // P_{n-1}(x)

        public LegendreEvaluation(double p, double pm) {
            this.p = p;
            this.pm = pm;
        } 
    }

    /**
     * Container for computed point and weight rule pair.
     */
    private static class RulePair {
        private final double point;
        private final double weight;

        public RulePair(double point, double weight) {
            this.point = point;
            this.weight = weight;
        } 
    }
}
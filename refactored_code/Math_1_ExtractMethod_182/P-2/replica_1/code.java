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
package org.apache.commons.math3.genetics;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;

/**
 * Cycle Crossover [CX] builds offspring from <b>ordered</b> chromosomes by identifying cycles
 * between two parent chromosomes. To form the children, the cycles are copied from the
 * respective parents.
 * <p>
 * To form a cycle the following procedure is applied:
 * <ol>
 *   <li>start with the first gene of parent 1</li>
 *   <li>look at the gene at the same position of parent 2</li>
 *   <li>go to the position with the same gene in parent 1</li>
 *   <li>add this gene index to the cycle</li>
 *   <li>repeat the steps 2-5 until we arrive at the starting gene of this cycle</li>
 * </ol>
 * The indices that form a cycle are then used to form the children in alternating order, i.e.
 * in cycle 1, the genes of parent 1 are copied to child 1, while in cycle 2 the genes of parent 1
 * are copied to child 2, and so forth ...
 * </p>
 *
 * Example (zero-start cycle):
 * <pre>
 * p1 = (8 4 7 3 6 2 5 1 9 0)    X   c1 = (8 1 2 3 4 5 6 7 9 0)
 * p2 = (0 1 2 3 4 5 6 7 8 9)    X   c2 = (0 4 7 3 6 2 5 1 8 9)
 *
 * cycle 1: 8 0 9
 * cycle 2: 4 1 7 2 5 6
 * cycle 3: 3
 * </pre>
 *
 * This policy works only on {@link AbstractListChromosome}, and therefore it
 * is parameterized by T. Moreover, the chromosomes must have same lengths.
 *
 * @see <a href="http://www.rubicite.com/Tutorials/GeneticAlgorithms/CrossoverOperators/CycleCrossoverOperator.aspx">
 * Cycle Crossover Operator</a>
 *
 * @param <T> generic type of the {@link AbstractListChromosome}s for crossover
 * @since 3.1
 * @version $Id$
 */
public class CycleCrossover<T> implements CrossoverPolicy {

    /** If the start index shall be chosen randomly. */
    private final boolean randomStart;

    /**
     * Creates a new {@link CycleCrossover} policy.
     */
    public CycleCrossover() {
        this(false);
    }

    /**
     * Creates a new {@link CycleCrossover} policy using the given {@code randomStart} behavior.
     *
     * @param randomStart whether the start index shall be chosen randomly or be set to 0
     */
    public CycleCrossover(final boolean randomStart) {
        this.randomStart = randomStart;
    }

    /**
     * Returns whether the starting index is chosen randomly or set to zero.
     *
     * @return {@code true} if the starting index is chosen randomly, {@code false} otherwise
     */
    public boolean isRandomStart() {
        return randomStart;
    }

    /**
     * {@inheritDoc}
     *
     * @throws MathIllegalArgumentException if the chromosomes are not an instance of {@link AbstractListChromosome}
     * @throws DimensionMismatchException if the length of the two chromosomes is different
     */
    @SuppressWarnings("unchecked")
    public ChromosomePair crossover(final Chromosome first, final Chromosome second)
        throws DimensionMismatchException, MathIllegalArgumentException {

        if (!(first instanceof AbstractListChromosome<?> && second instanceof AbstractListChromosome<?>)) {
            throw new MathIllegalArgumentException(LocalizedFormats.INVALID_FIXED_LENGTH_CHROMOSOME);
        }
        return mate((AbstractListChromosome<T>) first, (AbstractListChromosome<T>) second);
    }

    /**
     * Helper for {@link #crossover(Chromosome, Chromosome)}. Performs the actual crossover.
     *
     * @param first the first chromosome
     * @param second the second chromosome
     * @return the pair of new chromosomes that resulted from the crossover
     * @throws DimensionMismatchException if the length of the two chromosomes is different
     */
    protected ChromosomePair mate(final AbstractListChromosome<T> first, final AbstractListChromosome<T> second)
        throws DimensionMismatchException {

        final int length = first.getLength();
        checkLength(length, second.getLength());

        // array representations of the parents
        final List<T> parent1Rep = first.getRepresentation();
        final List<T> parent2Rep = second.getRepresentation();
        // and of the children: do a crossover copy to simplify the later processing
        final List<T> child1Rep = new ArrayList<T>(second.getRepresentation());
        final List<T> child2Rep = new ArrayList<T>(first.getRepresentation());

        // the set of all visited indices so far
        final Set<Integer> visitedIndices = new HashSet<Integer>(length);

        // determine the starting index
        int idx = calculateStartIdx(length);
        int cycle = 1;

        while (visitedIndices.size() < length) {
            final List<Integer> indices = findCycleIndices(idx, parent1Rep, parent2Rep);

            // for even cycles: swap the child elements on the indices found in this cycle
            if (cycle % 2 != 0) {
                swapElements(child1Rep, child2Rep, indices);
            }
            cycle++;

            visitedIndices.addAll(indices);
            // find next starting index: last one + 1 until we find an unvisited index
            idx = findNextStartIdx(indices.get(0), length, visitedIndices);
            indices.clear();
        }

        return new ChromosomePair(first.newFixedLengthChromosome(child1Rep),
                                  second.newFixedLengthChromosome(child2Rep));
    }

    /**
     * Helper method to verify both chromosomes have the same length.
     */
    private void checkLength(final int length1, final int length2) throws DimensionMismatchException {
        if (length1 != length2) {
            throw new DimensionMismatchException(length2, length1);
        }
    }

    /**
     * Determines the starting index for the cycle crossover.
     */
    private int calculateStartIdx(final int length) {
        return randomStart ? GeneticAlgorithm.getRandomGenerator().nextInt(length) : 0;
    }

    /**
     * Finds the indices belonging to a single cycle starting at the given index.
     */
    private List<Integer> findCycleIndices(final int startIdx, final List<T> parent1Rep, final List<T> parent2Rep) {
        final List<Integer> indices = new ArrayList<Integer>();
        indices.add(startIdx);

        T item = parent2Rep.get(startIdx);
        int idx = parent1Rep.indexOf(item);

        while (idx != startIdx) {
            indices.add(idx);
            item = parent2Rep.get(idx);
            idx = parent1Rep.indexOf(item);
        }
        return indices;
    }

    /**
     * Swaps elements of the two child representation lists at the specified indices.
     */
    private void swapElements(final List<T> child1Rep, final List<T> child2Rep, final List<Integer> indices) {
        for (int i : indices) {
            final T tmp = child1Rep.get(i);
            child1Rep.set(i, child2Rep.get(i));
            child2Rep.set(i, tmp);
        }
    }

    /**
     * Finds the next starting index that has not been visited yet.
     */
    private int findNextStartIdx(final int currentStartIdx, final int length, final Set<Integer> visitedIndices) {
        int idx = (currentStartIdx + 1) % length;
        while (visitedIndices.contains(idx) && visitedIndices.size() < length) {
            idx++;
            if (idx >= length) {
                idx = 0;
            }
        }
        return idx;
    }
}
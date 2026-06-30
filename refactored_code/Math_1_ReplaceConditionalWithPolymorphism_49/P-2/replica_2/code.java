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
package org.apache.commons.math3.fraction;

import java.math.BigInteger;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.exception.NullArgumentException;

/**
 * Formats a BigFraction number in proper format.  The number format for each of
 * the whole number, numerator and, denominator can be configured.
 * <p>
 * Minus signs are only allowed in the whole number part - i.e.,
 * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
 * will result in a <code>ParseException</code>.</p>
 * 
 * @since 1.1
 * @version $Id$
 */
public class ProperBigFractionFormat extends BigFractionFormat {

    /** Serializable version identifier */
    private static final long serialVersionUID = -6337346779577272307L;

    /** Character parsers map to replace conditional with polymorphism */
    private static final Map<Character, CharacterParser> PARSERS = new HashMap<Character, CharacterParser>();

    static {
        PARSERS.put((char) 0, new ZeroCharacterParser());
        PARSERS.put('/', new SlashCharacterParser());
    }

    /** The format used for the whole number. */
    private NumberFormat wholeFormat;

    /**
     * Create a proper formatting instance with the default number format for
     * the whole, numerator, and denominator.
     */
    public ProperBigFractionFormat() {
        this(getDefaultNumberFormat());
    }

    /**
     * Create a proper formatting instance with a custom number format for the
     * whole, numerator, and denominator.
     * @param format the custom format for the whole, numerator, and
     *        denominator.
     */
    public ProperBigFractionFormat(final NumberFormat format) {
        this(format, (NumberFormat)format.clone(), (NumberFormat)format.clone());
    }

    /**
     * Create a proper formatting instance with a custom number format for each
     * of the whole, numerator, and denominator.
     * @param wholeFormat the custom format for the whole.
     * @param numeratorFormat the custom format for the numerator.
     * @param denominatorFormat the custom format for the denominator.
     */
    public ProperBigFractionFormat(final NumberFormat wholeFormat,
                                   final NumberFormat numeratorFormat,
                                   final NumberFormat denominatorFormat) {
        super(numeratorFormat, denominatorFormat);
        setWholeFormat(wholeFormat);
    }

    /**
     * Formats a {@link BigFraction} object to produce a string.  The BigFraction
     * is output in proper format.
     *
     * @param fraction the object to format.
     * @param toAppendTo where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     *            offsets of the alignment field
     * @return the value passed in as toAppendTo.
     */
    @Override
    public StringBuffer format(final BigFraction fraction,
                               final StringBuffer toAppendTo, final FieldPosition pos) {

        pos.setBeginIndex(0);
        pos.setEndIndex(0);

        BigInteger num = fraction.getNumerator();
        BigInteger den = fraction.getDenominator();
        BigInteger whole = num.divide(den);
        num = num.remainder(den);

        if (!BigInteger.ZERO.equals(whole)) {
            getWholeFormat().format(whole, toAppendTo, pos);
            toAppendTo.append(' ');
            if (num.compareTo(BigInteger.ZERO) < 0) {
                num = num.negate();
            }
        }
        getNumeratorFormat().format(num, toAppendTo, pos);
        toAppendTo.append(" / ");
        getDenominatorFormat().format(den, toAppendTo, pos);

        return toAppendTo;
    }

    /**
     * Access the whole format.
     * @return the whole format.
     */
    public NumberFormat getWholeFormat() { 
        return wholeFormat; 
    }

    /**
     * Parses a string to produce a {@link BigFraction} object.  This method
     * expects the string to be formatted as a proper BigFraction.
     * <p>
     * Minus signs are only allowed in the whole number part - i.e.,
     * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
     * will result in a <code>ParseException</code>.</p>
     *
     * @param source the string to parse
     * @param pos input/ouput parsing parameter.
     * @return the parsed {@link BigFraction} object.
     */
    @Override
    public BigFraction parse(final String source, final ParsePosition pos) {
        // try to parse improper BigFraction
        BigFraction ret = super.parse(source, pos);
        if (ret != null) {
            return ret;
        }

        final int initialIndex = pos.getIndex();

        // parse whitespace
        parseAndIgnoreWhitespace(source, pos);

        // parse whole
        BigInteger whole = parseNextBigInteger(source, pos);
        if (whole == null) {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.setIndex(initialIndex);
            return null;
        }

        // parse whitespace
        parseAndIgnoreWhitespace(source, pos);

        // parse numerator
        BigInteger num = parseNextBigInteger(source, pos);
        if (num == null) {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.setIndex(initialIndex);
            return null;
        }

        if (num.compareTo(BigInteger.ZERO) < 0) {
            // minus signs should be leading, invalid expression
            pos.setIndex(initialIndex);
            return null;
        }

        // parse '/'
        final int startIndex = pos.getIndex();
        final char c = parseNextCharacter(source, pos);

        CharacterParser parser = PARSERS.get(c);
        if (parser == null) {
            parser = new DefaultCharacterParser();
        }

        return parser.handle(this, source, pos, initialIndex, startIndex, whole, num);
    }

    /**
     * Modify the whole format.
     * @param format The new whole format value.
     * @throws NullArgumentException if {@code format} is {@code null}.
     */
    public void setWholeFormat(final NumberFormat format) {
        if (format == null) {
            throw new NullArgumentException(LocalizedFormats.WHOLE_FORMAT);
        }
        this.wholeFormat = format;
    }

    /**
     * CharacterParser interface for polymorphism.
     */
    private interface CharacterParser {
        BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos,
                           int initialIndex, int startIndex, BigInteger whole, BigInteger num);
    }

    private static class ZeroCharacterParser implements CharacterParser {
        @Override
        public BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos,
                                   int initialIndex, int startIndex, BigInteger whole, BigInteger num) {
            // no '/'
            // return num as a BigFraction
            return new BigFraction(num);
        }
    }

    private static class SlashCharacterParser implements CharacterParser {
        @Override
        public BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos,
                                   int initialIndex, int startIndex, BigInteger whole, BigInteger num) {
            // found '/', continue parsing denominator
            // parse whitespace
            formatter.parseAndIgnoreWhitespace(source, pos);

            // parse denominator
            final BigInteger den = formatter.parseNextBigInteger(source, pos);
            if (den == null) {
                // invalid integer number
                // set index back to initial, error index should already be set
                // character examined.
                pos.setIndex(initialIndex);
                return null;
            }

            if (den.compareTo(BigInteger.ZERO) < 0) {
                // minus signs must be leading, invalid
                pos.setIndex(initialIndex);
                return null;
            }

            boolean wholeIsNeg = whole.compareTo(BigInteger.ZERO) < 0;
            BigInteger processedWhole = whole;
            if (wholeIsNeg) {
                processedWhole = whole.negate();
            }
            BigInteger processedNum = processedWhole.multiply(den).add(num);
            if (wholeIsNeg) {
                processedNum = processedNum.negate();
            }

            return new BigFraction(processedNum, den);
        }
    }

    private static class DefaultCharacterParser implements CharacterParser {
        @Override
        public BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos,
                                   int initialIndex, int startIndex, BigInteger whole, BigInteger num) {
            // invalid '/'
            // set index back to initial, error index should be the last
            // character examined.
            pos.setIndex(initialIndex);
            pos.setErrorIndex(startIndex);
            return null;
        }
    }
}
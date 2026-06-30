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

    /** The format used for the whole number. */
    private NumberFormat wholeFormat;

    /** Polymorphic formatters for handling zero and non-zero whole parts */
    private static final Map<Boolean, WholeNumberFormatter> FORMATTERS = new HashMap<Boolean, WholeNumberFormatter>();

    /** Polymorphic handlers for characters parsed after the numerator */
    private static final Map<Character, SlashHandler> HANDLERS = new HashMap<Character, SlashHandler>();

    static {
        FORMATTERS.put(Boolean.TRUE, new ZeroWholeFormatter());
        FORMATTERS.put(Boolean.FALSE, new NonZeroWholeFormatter());

        HANDLERS.put((char) 0, new ZeroHandler());
        HANDLERS.put('/', new SlashHandlerImpl());
    }

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

        WholeNumberFormatter wholeFormatter = FORMATTERS.get(BigInteger.ZERO.equals(whole));
        num = wholeFormatter.formatWhole(this, whole, num, toAppendTo, pos);

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

        SlashHandler handler = HANDLERS.get(c);
        if (handler == null) {
            handler = new DefaultHandler();
        }
        return handler.handle(this, source, pos, initialIndex, startIndex, whole, num);
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
     * Strategy interface for formatting whole numbers.
     */
    private interface WholeNumberFormatter {
        BigInteger formatWhole(ProperBigFractionFormat formatter, BigInteger whole, BigInteger num, StringBuffer toAppendTo, FieldPosition pos);
    }

    private static class NonZeroWholeFormatter implements WholeNumberFormatter {
        public BigInteger formatWhole(ProperBigFractionFormat formatter, BigInteger whole, BigInteger num, StringBuffer toAppendTo, FieldPosition pos) {
            formatter.getWholeFormat().format(whole, toAppendTo, pos);
            toAppendTo.append(' ');
            if (num.compareTo(BigInteger.ZERO) < 0) {
                return num.negate();
            }
            return num;
        }
    }

    private static class ZeroWholeFormatter implements WholeNumberFormatter {
        public BigInteger formatWhole(ProperBigFractionFormat formatter, BigInteger whole, BigInteger num, StringBuffer toAppendTo, FieldPosition pos) {
            return num;
        }
    }

    /**
     * Strategy interface for handling the character parsed after the numerator.
     */
    private interface SlashHandler {
        BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos, int initialIndex, int startIndex, BigInteger whole, BigInteger num);
    }

    private static class ZeroHandler implements SlashHandler {
        public BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos, int initialIndex, int startIndex, BigInteger whole, BigInteger num) {
            return new BigFraction(num);
        } 
    }

    private static class SlashHandlerImpl implements SlashHandler {
        public BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos, int initialIndex, int startIndex, BigInteger whole, BigInteger num) {
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

            BigInteger localWhole = whole;
            BigInteger localNum = num;
            boolean wholeIsNeg = localWhole.compareTo(BigInteger.ZERO) < 0;
            if (wholeIsNeg) {
                localWhole = localWhole.negate();
            }
            localNum = localWhole.multiply(den).add(localNum);
            if (wholeIsNeg) {
                localNum = localNum.negate();
            }

            return new BigFraction(localNum, den);
        }
    }

    private static class DefaultHandler implements SlashHandler {
        public BigFraction handle(ProperBigFractionFormat formatter, String source, ParsePosition pos, int initialIndex, int startIndex, BigInteger whole, BigInteger num) {
            // invalid '/'
            // set index back to initial, error index should be the last
            // character examined.
            pos.setIndex(initialIndex);
            pos.setErrorIndex(startIndex);
            return null;
        }
    }
}
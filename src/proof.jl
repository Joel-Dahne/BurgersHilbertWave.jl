"""
    round_for_publishing(n₀, δ₀, D₀; sigdigits = 10)

Convert `n₀, δ₀, D₀` to `Float64`, rounding up to the prescribed
number of significant digits, and check that the inequality `δ₀ <= (1
- D₀)^2 / 4n₀` holds for the rounded values as well.

This is used to get upper bounds of the values in a simpler format
than the `Arb` type.
"""
function round_for_publishing(n₀::Arb, δ₀::Arb, D₀::Arb; sigdigits = 10)
    n₀_float = Arblib.get_d(ubound(n₀), RoundUp)
    δ₀_float = Arblib.get_d(ubound(δ₀), RoundUp)
    D₀_float = Arblib.get_d(ubound(D₀), RoundUp)

    # Check that the inequality holds before rounding. Conversion to
    # Float64 loses precision so this is not guaranteed.
    inequality_holds =
        D₀_float < 1 && Arb(δ₀_float) < (1 - Arb(D₀_float))^2 / 4Arb(n₀_float)

    if !inequality_holds
        @warn "Inequality doesn't hold after conversion" n₀_float, δ₀_float, D₀_float
        return false, n₀_float, δ₀_float, D₀_float
    end

    n₀_float_rounded = round(n₀_float, RoundUp; sigdigits)
    δ₀_float_rounded = round(δ₀_float, RoundUp; sigdigits)
    D₀_float_rounded = round(D₀_float, RoundUp; sigdigits)

    @assert n₀ <= n₀_float_rounded
    @assert δ₀ <= δ₀_float_rounded
    @assert D₀ <= D₀_float_rounded

    # Check that the inequality holds after rounding
    inequality_holds =
        D₀_float_rounded < 1 &&
        Arb(δ₀_float_rounded) < (1 - Arb(D₀_float_rounded))^2 / 4Arb(n₀_float_rounded)

    return inequality_holds, n₀_float_rounded, δ₀_float_rounded, D₀_float_rounded
end

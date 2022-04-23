export iswide

"""
    mince(x::Arb, n::Integer)

Return a vector with `n` balls covering the ball `x`.
"""
function mince(x::Arb, n::Integer)
    balls = Vector{Arb}(undef, n)
    xₗ, xᵤ = Arblib.getinterval(Arb, x)
    dx = (xᵤ - xₗ) / n
    for i in eachindex(balls)
        yₗ = xₗ + (i - 1) * dx
        yᵤ = xₗ + i * dx
        balls[i] = Arb((yₗ, yᵤ))
    end

    return balls
end

"""
    iswide(x; cutoff = 10)

Return true if `x` is wide in the meaning that the effective relative
accuracy of `x` measured in bits is more than `cutoff` lower than it's
precision. For `x` not of type `Arb` or `Acb` this always return
`false`.
"""
iswide(x::Union{Arb,Acb}; cutoff = 10) = Arblib.rel_accuracy_bits(x) < precision(x) - cutoff
iswide(::Number; cutoff = 10) = false

"""
    stieltjes(T, n::Integer)

Compute the Stieltjes constant `γₙ` in type `T`.
"""
function stieltjes(::Type{Arb}, n::Integer)
    res = zero(Acb)
    a = one(Acb)

    ccall(
        (:acb_dirichlet_stieltjes, Arblib.libarb),
        Cvoid,
        (Ref{Arblib.acb_struct}, Ref{fmpz_struct}, Ref{Arblib.acb_struct}, Clong),
        res,
        fmpz_struct(convert(Int, n)),
        a,
        precision(res),
    )

    return real(res)
end
stieltjes(T, n::Integer) = convert(T, stieltjes(Arb, n))

"""
    unique_integer(x::Arb)

If `x` contains a unique integer return `true, n` where `n` is the
integer. Otherwise return `false, 0`
"""
function unique_integer(x::Arb)
    res = fmpz_struct()
    unique = ccall(
        Arblib.@libarb(arb_get_unique_fmpz),
        Int,
        (Ref{fmpz_struct}, Ref{Arblib.arb_struct}),
        res,
        x,
    )

    return !iszero(unique), Int(res)
end

"""
    abs(x::ArbSeries)

Compute the absolute value of `x`.

If `x[0]` contains zero then all non-constant terms are set to `NaN`,
otherwise either `-x` or `x` is returned depending on the sign of
`x[0]`.
"""
function Base.abs(x::ArbSeries)
    if Arblib.contains_zero(Arblib.ref(x, 0))
        res = zero(x)
        Arblib.abs!(Arblib.ref(res, 0), Arblib.ref(x, 0))
        for i = 1:Arblib.degree(x)
            Arblib.indeterminate!(Arblib.ref(res, i))
        end
        # Since we manually set the coefficients of the polynomial we
        # need to also manually set the degree.
        res.arb_poly.length = Arblib.degree(x) + 1
        return res
    elseif Arblib.isnegative(Arblib.ref(x, 0))
        return -x
    else
        return copy(x)
    end
end

"""
    abspow(x, y)

Compute `abs(x)^y `in a way that works if `x` overlaps with zero.
"""
function abspow(x::Arb, y::Arb)
    iszero(y) && return one(x)

    if iszero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(zero(x))
        Arblib.ispositive(y) && return zero(x)
        return Arblib.unit_interval!(zero(x))
    end

    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(zero(x))
        x_upp = Arblib.abs_ubound(Arb, x)
        return Arb((zero(x), x_upp^y))
    end

    res = abs(x)
    return Arblib.pow!(res, res, y)
end

function abspow(x::ArbSeries, y::Arb)
    if Arblib.contains_zero(Arblib.ref(x, 0))
        # All non-constant terms are indeterminate, the constant term
        # is given by abs(x[0])^y
        res = ArbSeries(abspow(x[0], y), degree = Arblib.degree(x))
        for i = 1:Arblib.degree(x)
            Arblib.indeterminate!(Arblib.ref(res, i))
        end
        # Since we manually set the coefficients of the polynomial we
        # need to also manually set the degree.
        res.arb_poly.length = length(x)
        return res
    end

    res = abs(x)
    return Arblib.pow_arb_series!(res, res, y, length(res))
end

abspow(x, y) = abs(x)^y

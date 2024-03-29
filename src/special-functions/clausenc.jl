export clausenc, clausencmzeta

"""
    _reduce_argument_clausen(x::Arb)

Perform argument reduction for the Clausen functions. It returns `(y,
haszero, haspi, has2pi)` where `y` is the reduced argument, `haszero`,
and `has2pi` are true if `y` overlaps with zero, and `2π`
respectively, `haspi` is true if `x` overlaps with `π` **or** `-π`.

The possible results are divided into three cases
1. `x` covers a full period, i.e. the radius is at least `π`. In this
  case it returns the full interval ``[0, 2π]`` and `haszero, haspi,
  has2pi` are all true.
2. `x` doesn't contain a multiple of `2π`. In this case the result
  satisfies `0 < y < 2π`, `haszero` and `has2pi` are false and `haspi`
  is true if `y` contains `π`.
3. `x` contains precisely one multiple of `2π`. In this case `-2π < y
  < 2π` and `haszero` is true while `has2pi` is false, `haspi` is true
  if `y` contains `-π` or `π`.

Note that if `haszero` is false then `y` is guaranteed to satisfy `0 <
y < 2π`. This is used by many of the functions.

It has issues with very small but negative values of `x`. In this case
`x + 2π < 2π` in theory but due to overestimations we might get that
they overlap. For this reason it can be beneficial to pass the
absolute value of `x` to this function and handle the sign outside of
it.
"""
function _reduce_argument_clausen(x::Arb)
    pi = Arb(π)
    twopi = let tmp = Arb(π)
        Arblib.mul_2exp!(tmp, tmp, 1)
    end

    if Arblib.ispositive(x) && x < twopi
        # Happy path, when 0 < x < 2π
        y = copy(x)
        haszero = has2pi = false
        haspi = Arblib.overlaps(x, pi)
    else
        xdiv2pi = x / twopi

        if Arblib.cmp_2exp(Arblib.radref(xdiv2pi), -1) >= 0
            # If the radius of x / 2π is not less than 1 / 2 then y
            # spans the full interval [0, 2π]
            y = Arb((0, twopi))
            haszero = haspi = has2pi = true
        else
            # Add/subtract a multiple of 2π so that the midpoint is on
            # the interval [0, 2π]
            k = Arblib.floor!(zero(Arf), Arblib.midref(xdiv2pi))
            y = x - k * twopi

            # We want to avoid the case when y overlaps 2π
            if Arblib.overlaps(y, twopi)
                y -= twopi
            end

            haszero = Arblib.contains_zero(y)
            haspi = Arblib.overlaps(y, pi) || Arblib.overlaps(y, -pi)
            has2pi = Arblib.overlaps(y, twopi)

            # Neither of these checks should be required
            # mathematically, though overestimations might lead to
            # issues here?

            # We guarantee this
            @assert haszero || 0 < y < twopi "$x $y"
            # It should never contain -2π or 2π in this case
            @assert !(Arblib.overlaps(y, -twopi) || has2pi)
        end
    end

    return y, haszero, haspi, has2pi
end

"""
    _clausenc_polylog(x::Arb, s::Union{Arb,Integer})
    _clausenc_polylog(x::Arb, s::ArbSeries)
    _clausenc_polylog(x::Arb, s::Arb, β::Integer)

Evaluation of the `clausenc` function through the polylog function.

It uses the formula
```
clausenc(x, s) = real(polylog(s, exp(im * x)))
```
"""
function _clausenc_polylog(x::Arb, s::Union{Arb,Integer})
    z = exp(Acb(0, x, prec = precision(x)))
    s = s isa Integer ? s : Acb(s, prec = precision(x))
    return real(polylog(s, z))
end

function _clausenc_polylog(x::Arb, s::ArbSeries)
    z = exp(Acb(0, x, prec = precision(x)))
    return ArbSeries(real.(Arblib.coeffs(polylog(AcbSeries(s), z))))
end

function _clausenc_polylog(x::Arb, s::Arb, β::Integer)
    z = exp(Acb(0, x, prec = precision(x)))
    return real(polylog(AcbSeries([s, 1], degree = β), z)[β]) * factorial(β)
end

"""
    _clausenc_zeta(x::Arb, s::Arb)

Evaluation of the `clausenc` function through the zeta function.

It is based on [`equation_clausenc_zeta`](@ref), which we can also
write as
```
clausenc(x, s) = let v = 1 - s
    gamma(v) * inv(2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
end
```
It comes from the formula [25.13.2](https://dlmf.nist.gov/25.13) for
the periodic zeta function and then taking the real part.

The formula is well defined as long as `s` doesn't overlap with any
non-negative integer. See further down for how to handle the case when
`s` overlaps a unique non-negative integer. If `s` overlaps with more
than one non-negative integer we just get an indeterminate result.

It currently only handles `0 < x < 2π`.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using [`enclosure_series`](@ref).

# Handling `s` overlapping zero
If `s` is zero then `zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` diverges
but `cospi(v / 2)` goes to zero. To compute their product we use the
following approach. Let `zeta_deflated(s, a) = zeta(s, a) + 1 / (1 -
s)` be the deflated zeta function. We can then write the product as
```
cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) =
    cospi(v / 2) * (zeta_deflated(v, x / 2π) + zeta_deflated(v, 1 - x / 2π))
    + 2cospi(v / 2) / (v - 1)
```
The first term is well defined and can be evaluated directly. The
second term has a removable singularity and is handled by rewriting it
as `2cospi((t + 1) / 2) / t` to place the removable singularity at `t
= 0` and then use [`fx_div_x`](@ref).

# Handling `s` overlapping a positive integer
If `s` is a positive integer then `gamma(v)` diverges, if the
integer is even then `cospi(v / 2)` is zero and if the integer is odd
then `zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` is zero. To see that
`zeta(v, x / 2π) + zeta(v, 1 - x / 2π)` is zero when `s` is an odd
non-negative integer, i.e. when `v` is an even non-positive integer,
we can use formula [25.11.14](https://dlmf.nist.gov/25.11.E14)
together with [24.4.1](https://dlmf.nist.gov/24.4.E1).

For even `s` we thus want to compute an enclosure of
```
gamma(v) * cospi(v / 2)
```
and for odd `s` we want to compute an enclosure of
```
gamma(v) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
```
In both cases rewrite them using the reciprocal gamma function
[`rgamma`](@ref) as
```
cospi(v / 2) / rgamma(v)

(zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / rgamma(v)
```
To be able to use [`fx_div_x`](@ref) we let `t = v - v_integer`, where
`v_integer` is the integer that `v` overlaps with, and write them as
```
(cospi((t + v_integer) / 2) / t) / (rgamma(t + v_integer) / t)

((zeta(t + v_integer, x / 2π) + zeta(t + v_integer, 1 - x / 2π)) / t) / (rgamma(t + v_integer) / t)
```
"""
function _clausenc_zeta(x::Arb, s::Arb)
    # Check that x > 0
    Arblib.ispositive(x) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    inv2pi = let tmp = Arb(π, prec = precision(x))
        Arblib.mul_2exp!(tmp, tmp, 1)
        Arblib.inv!(tmp, tmp)
    end
    xinv2pi = x * inv2pi
    onemxinv2pi = let onemxinv2pi = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(onemxinv2pi, xinv2pi)
        Arblib.add!(onemxinv2pi, onemxinv2pi, 1)
    end

    # Check that 1 - x / 2π > 0, i.e. x < 2π
    Arblib.ispositive(onemxinv2pi) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    v = let v = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(v, s)
        Arblib.add!(v, v, 1)
    end

    unique, s_integer = unique_integer(s)

    if unique && s_integer == 0
        # Enclosure of cospi(v / 2) / (v - 1)
        cos_div_v = fx_div_x(t -> cospi((t + 1) / 2), v - 1, extra_degree = 2)

        # Enclosure of zeta_deflated(v, x / 2π) + zeta_deflated(v, 1 - x / 2π)
        z = zeta_deflated(v, xinv2pi) + zeta_deflated(v, onemxinv2pi)

        rest = ArbExtras.enclosure_series(v -> gamma(v) * inv2pi^v, v, degree = 2)

        return rest * (cospi(v / 2) * z + 2cos_div_v)
    elseif unique && s_integer > 0
        v_integer = 1 - s_integer

        # Enclosure of rgamma(v) / (v - v_integer)
        rgamma_div_v = fx_div_x(t -> rgamma(t + v_integer), v - v_integer, extra_degree = 2)

        if iseven(s_integer)
            # Enclosure of cospi(v / 2) / (v - v_integer)
            cos_div_v =
                fx_div_x(t -> cospi((t + v_integer) / 2), v - v_integer, extra_degree = 2)

            # Enclosure of gamma(v) * cospi(v / 2)
            gammacos = cos_div_v / rgamma_div_v

            if iswide(s)
                rest = ArbExtras.enclosure_series(
                    v -> inv2pi^v * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi)),
                    v,
                    degree = 2,
                )
            else
                rest = inv2pi^v * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))
            end

            return gammacos * rest
        else
            # Enclosure of (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / (v - v_integer)
            zeta_div_v = fx_div_x(
                t -> zeta(t + v_integer, xinv2pi) + zeta(t + v_integer, onemxinv2pi),
                v - v_integer,
                extra_degree = 2,
                force = true,
            )

            # Enclosure of gamma(v) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
            gammazeta = zeta_div_v / rgamma_div_v

            if iswide(s)
                rest =
                    ArbExtras.enclosure_series(v -> inv2pi^v * cospi(v / 2), v, degree = 2)
            else
                rest = inv2pi^v * cospi(v / 2)
            end

            return gammazeta * rest
        end
    else
        f(v) =
            gamma(v) * inv2pi^v * cospi(v / 2) * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))

        if iswide(s)
            return ArbExtras.enclosure_series(f, v, degree = 2)
        else
            return f(v)
        end
    end
end

"""
    _clausenc_zeta(x::Acb, s::Arb)

Evaluation of the `clausenc` function through the zeta function.

This uses the same formula as the method with `x::Arb`, valid for `0 <
real(x) < 2π`. If `s` overlaps with any non-negative integer the
result will be indeterminate.

This method is only used in the integration for bounding the error
term. It is therefore not as optimized as many of the other methods.
"""
function _clausenc_zeta(x::Acb, s::Arb)
    twopi = 2Arb(π)

    0 < real(x) < twopi || throw(
        DomainError(x, "method only supports x with real part on the interval (0, 2π)"),
    )

    inv2pi = inv(twopi)
    xdiv2pi = x / twopi
    v = 1 - s

    return ArbExtras.enclosure_series(
        v -> gamma(v) * inv2pi^v * cospi(v / 2),
        v,
        degree = 2,
    ) * (zeta(Acb(v), xdiv2pi) + zeta(Acb(v), 1 - xdiv2pi))
end

"""
    _clausenc_zeta(x::Arb, s::ArbSeries)

Evaluation of the `clausenc` function through the zeta function as a
power series in `s`.

It supports `s` overlapping non-negative integers in the same way
`_clausenc_zeta(x::Arb, s::Arb)` does, using [`fx_div_x`](@ref) to
handle the removable singularities.
"""
function _clausenc_zeta(x::Arb, s::ArbSeries)
    # Handle the case when s has degree 0, so that we can assume that
    # the degree is positive from here on.
    iszero(Arblib.degree(s)) && return ArbSeries(_clausenc_zeta(x, s[0]))

    # Check that x > 0
    Arblib.ispositive(x) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    inv2pi = inv(2Arb(π, prec = precision(x)))
    xinv2pi = x * inv2pi
    onemxinv2pi = let onemxinv2pi = zero(x) # We do it like this to preserve the precision
        Arblib.neg!(onemxinv2pi, xinv2pi)
        Arblib.add!(onemxinv2pi, onemxinv2pi, 1)
    end

    # Check that 1 - x / 2π > 0, i.e. x < 2π
    Arblib.ispositive(onemxinv2pi) ||
        throw(DomainError(x, "method only supports x on the interval (0, 2π)"))

    v = 1 - s

    s0 = s[0]

    unique, s0_integer = unique_integer(s0)

    if unique && s0_integer == 0
        cos_div_v = fx_div_x(t -> cospi((t + 1) / 2), v - 1, extra_degree = 1)

        return gamma(v) *
               inv2pi^v *
               (
                   cospi(v / 2) *
                   (zeta_deflated(v, xinv2pi) + zeta_deflated(v, onemxinv2pi)) + 2cos_div_v
               )
    elseif unique && s0_integer > 0
        v0_integer = 1 - s0_integer

        # Enclosure of rgamma(v) / (v - v_integer)
        rgamma_div_v =
            fx_div_x(t -> rgamma(t + v0_integer), v - v0_integer, extra_degree = 2)

        if iseven(s0_integer)
            # Enclosure of cospi(v / 2) / (v - v_integer)
            cos_div_v =
                fx_div_x(t -> cospi((t + v0_integer) / 2), v - v0_integer, extra_degree = 2)

            # Enclosure of gamma(v) * cospi(v / 2)
            gammacos = cos_div_v / rgamma_div_v

            rest = inv2pi^v * (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))

            return gammacos * rest
        else
            # Enclosure of (zeta(v, x / 2π) + zeta(v, 1 - x / 2π)) / (v - v_integer)
            zeta_div_v = fx_div_x(
                t -> zeta(t + v0_integer, xinv2pi) + zeta(t + v0_integer, onemxinv2pi),
                v - v0_integer,
                extra_degree = 2,
                force = true,
            )

            # Enclosure of gamma(v) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
            gammazeta = zeta_div_v / rgamma_div_v

            rest = inv2pi^v * cospi(v / 2)

            return gammazeta * rest
        end
    else
        return gamma(v) *
               inv2pi^v *
               cospi(v / 2) *
               (zeta(v, xinv2pi) + zeta(v, onemxinv2pi))
    end
end

"""
    _clausenc_zeta(x::Float64, s::Float64)

Evaluation of the `clausenc` function through the zeta function.

This uses the same formula as the method with `x::Arb`. But it doesn't
have any special handling when `s` is an integer, it will return
`NaN`.
"""
function _clausenc_zeta(x::Float64, s::Float64)
    x = mod2pi(x)
    v = 1 - s
    return gamma(v) / (2π)^v * cospi(v / 2) * (zeta(v, x / 2π) + zeta(v, 1 - x / 2π))
end

"""
    _clausenc_zeta(x::Arb, s::Arb, β::Integer)

Evaluation of the `clausenc(x, s, β)` function through the zeta
function.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure using a Taylor expansion in `s`.
**IMPROVE:** The degree used in this expansion could be tuned more. At
the moment it is set to a quite high value, but that seems to be
beneficial for `KdVZeroansatz` with a wide interval for `α`. Possibly
we want to use a higher value for higher values of `s` and also wider
values of `s`.
"""
function _clausenc_zeta(x::Arb, s::Arb, β::Integer)
    if iswide(s)
        f(s::Arb) = _clausenc_zeta(x, ArbSeries((s, 1), degree = β))[β] * factorial(β)
        f(s::ArbSeries) = begin
            @assert !iszero(Arblib.degree(s)) # Doesn't work in this case
            Arblib.derivative(
                _clausenc_zeta(x, ArbSeries(s, degree = Arblib.degree(s) + β)),
                β,
            )
        end

        res = ArbExtras.enclosure_series(f, s, degree = 10)
    else
        res = _clausenc_zeta(x, ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end

    return res
end

"""
    clausenc(x, s)

Compute the Clausen function ``C_s(x)``.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

To compute an good enclosure when `x` is a wide ball it uses that the
function is monotone for `0 < x < π` and even around `x = π`. The
monotonicity is a consequence of [`lemma_clausenc_monotone`](@ref).

If `x` contains zero and we don't have `s > 1` it returns an
indeterminate result.
- **IMPROVE:** Use that for even negative integers it is exactly zero.

If `x` contains zero and `s > 1` it computes an enclosure using that
the extrema is attained at `x = 0`, `abs_ubound(x)` or `π`, where we
have `clausenc(0, s) = zeta(s)` and `clausenc(π, s) = -eta(s)`.

If `x` doesn't contain zero we are assured by
[`_reduce_argument_clausen`](@ref) that `0 < x < 2π`. If `x` is a wide
ball, as determined by [`iswide`](@ref), it uses that the extrema is
attained at `lbound(x)`, `ubound(x)` or `π`. In this case it computes
the endpoints at a reduced precision given by
```
prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))
```
where `min_prec` is `32` in general but `64` if `s` is close to an
integer, determined by checking if the midpoint is withing `1e-2` of
an integer, in which case higher precision is typically needed.
- **IMPROVE:** This could be tuned more, but is probably not needed.
"""
function clausenc(x::Arb, s::Arb)
    # The function is even and this avoids issues with very small
    # negative x
    x = abs(x)

    x, haszero, haspi, has2pi = _reduce_argument_clausen(x)

    @assert !(has2pi && !haszero)

    # Handle the special case when x contains zero
    if haszero
        if !(s > 1)
            return Arblib.indeterminate!(zero(x))
        elseif haspi
            # Extrema at x = 0 and x = π
            return union(zeta(s), -Arblib.realref(eta(Acb(s))))
        elseif iszero(x)
            return zeta(s)
        else
            # Extrema at x = 0 and upper bound of abs(x)
            return union(zeta(s), clausenc(ArbExtras.enclosure_ubound(abs(x)), s))
        end
    end

    # We can now assume that 0 < x < 2π

    if iswide(x)
        s_f64 = Float64(s)
        min_prec = abs(s_f64 - round(s_f64)) < 1e-2 ? 64 : 32

        prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))

        xₗ, xᵤ = ArbExtras.enclosure_getinterval(setprecision(x, prec))
        res = union(clausenc(xₗ, s), clausenc(xᵤ, s))
        if haspi
            res = union(res, clausenc(Arb(π; prec), s))
        end
        return setprecision(res, precision(x))
    end

    return _clausenc_zeta(x, s)
end

function clausenc(x::Acb, s::Arb)
    x_real, haszero, _, _ = _reduce_argument_clausen(real(x))
    x = Acb(x_real, imag(x))

    # If x overlaps zero or s s is a non-negative use polylog formulation
    if haszero || (isinteger(s) && Arblib.isnonnegative(s))
        s = Acb(s)
        return (polylog(s, exp(im * x)) + polylog(s, exp(-im * x))) / 2
    else
        return _clausenc_zeta(x, s)
    end
end

function clausenc(x::Float64, s::Float64)
    if isinteger(s) && s >= 0
        return Float64(clausenc(Arb(mod2pi(x)), Arb(s)))
    else
        return _clausenc_zeta(x, s)
    end
end

clausenc(x::S, s::T) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausenc(convert(Arb, x), convert(Arb, s)))

"""
    clausenc(x::ArbSeries, s)

Compute the Taylor series of the Clausen function ``C_s(x)``.

It's computed by directly computing the Taylor coefficients by
differentiating `clausenc` and then composing with `x`.
"""
function clausenc(x::ArbSeries, s)
    x₀ = _reduce_argument_clausen(x[0])[1]

    res = zero(x)
    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausenc(x₀, s - i) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * clausens(x₀, s - i) / factorial(i)
        end
    end

    # Compose the Taylor series for the result with that of the input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

"""
    clausenc(x::Arb, s::ArbSeries)

Compute the Taylor series of the Clausen function ``C_s(x)`` in the
parameter `s`.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero it returns an indeterminate result. For `s > 1`
it would be possible to compute a finite result, but this has not been
needed.
"""
function clausenc(x::Arb, s::ArbSeries)
    # The function is even and this avoids issues with very small
    # negative x
    x = abs(x)

    x, haszero, _, _ = _reduce_argument_clausen(x)

    if haszero
        res = zero(s)
        for i = 0:Arblib.degree(res)
            res[0] = Arblib.indeterminate!(zero(x))
        end
        return res
    else
        return _clausenc_zeta(x, s)
    end
end

"""
    clausenc(x, s, β)

Compute ``C_s^{(β)}(x)``, that is `clausenc(x, s)` differentiated `β`
times w.r.t. `s`.

It first performs an argument reduction of `x`, using that the
function is `2π` periodic, with [`_reduce_argument_clausen`](@ref).

If `x` contains zero and we don't have `s > 1` it returns an
indeterminate result. If `s > 1` and `x`is not exactly zero it
computes an extremely naive enclosure using that `abs(clausenc(x, s,
β)) < abs(clausenc(0, s, β))`, that this is the case is easily seen
from the Fourier series.
- **IMPROVE:** Compute a tighter enclosure in this case. Either
  checking the derivative at the endpoint, similarly to
  [`clausens`](@ref), or expanding at `x = 0`.

If `x` is a wide ball (not containing zero), as determined by
`iswide(x)`, it computes a tighter enclosure by first checking if the
derivative doesn't contains zero, if not it uses monotonicity to only
evaluate at endpoints. If the derivative does contain zero it uses a
zeroth order approximation instead. In the wide case it computes the
endpoints at a reduced precision given by
```
prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))
```
where `min_prec` is `32` in general but `64` if `s` is close to an
integer, determined by checking if the midpoint withing `1e-2` of an
integer, in which case higher precision is typically needed.
"""
function clausenc(x::Arb, s::Arb, β::Integer)
    # The function is even and this avoids issues with very small
    # negative x
    x = abs(x)

    x, haszero, haspi, has2pi = _reduce_argument_clausen(x)

    @assert !(has2pi && !haszero)

    # Handle the special case when x contains zero
    if haszero
        if !(s > 1)
            # Only finite for s > 1
            return Arblib.indeterminate!(zero(x))
        else
            # Value at x = 0
            z = isone(β) ? dzeta(s) : zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
            if haspi
                # Absolute value upper bounded by value at x = 0
                return union(-z, z)
            elseif iszero(x)
                return z
            else
                # IMPROVE: Expand at x = 0? Check derivative at abs_ubound(x)?
                # Absolute value upper bounded by value at x = 0
                return union(-z, z)
            end
        end
    end

    # We can now assume that 0 < x < 2π

    if iswide(x)
        orig_prec = precision(x)
        s_f64 = Float64(s)
        min_prec = abs(s_f64 - round(s_f64)) < 1e-2 ? 64 : 32
        prec = min(max(Arblib.rel_accuracy_bits(x) + min_prec, min_prec), precision(x))
        x = setprecision(x, prec)

        # Compute derivative. We here call _clausens_zeta directly to
        # avoid potentially infinite recursion. This is okay since 0 <
        # x < 2π.
        dclausenc = -_clausens_zeta(x, s - 1, β)
        if Arblib.contains_zero(dclausenc)
            # Use a zero order approximation
            mid = Arblib.midpoint(Arb, x)
            res = Arblib.add_error!(clausenc(mid, s, β), (x - mid) * dclausenc)
        else
            # Use that it's monotone
            xₗ, xᵤ = ArbExtras.enclosure_getinterval(x)
            res = union(clausenc(xₗ, s, β), clausenc(xᵤ, s, β))
        end
        return setprecision(res, orig_prec)
    end

    return _clausenc_zeta(x, s, β)
end

clausenc(x::S, s::T, β::Integer) where {S<:Real,T<:Real} =
    convert(float(promote_type(S, T)), clausenc(convert(Arb, x), convert(Arb, s), β))

"""
    clausenc(x::ArbSeries, s, β)

Compute the Taylor series with respect to `x` of ``C_s^{(β)}(x)``,
that is `clausenc(x, s)` differentiated `β` times w.r.t. `s`.

It's computed by directly computing the Taylor coefficients by
differentiating ``C_s^{(\beta)}`` and then composing with `x`.
"""
function clausenc(x::ArbSeries, s, β::Integer)
    res = zero(x)
    x₀ = x[0]

    for i = 0:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausenc(x₀, s - i, β) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * clausens(x₀, s - i, β) / factorial(i)
        end
    end

    # Compose the Taylor series for the Clausian with that of the
    # input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

"""
    clausenc_expansion(x, s, M::Integer)

Compute the asymptotic expansion of `clausenc(x, s)` at `x = 0` up to
order `2M - 2`, meaning that the error term is of order `2M`.

This is the expansion given in [`equation_clausenc_expansion`](@ref).
The remainder term is bounded using
[`clausenc_expansion_remainder`](@ref).

It returns four things, the coefficient `C` and exponent `e` for the
non-analytic term, the analytic terms as a `ArbSeries` `P` and the
remainder term `E`.

It satisfies that
```
clausenc(y, s) ∈ C * abs(y)^e + P(y) + E * y^(2M)
```
for all `abs(y) <= abs(x)`.

If `skip_constant = true` it doesn't compute the constant term in the
expansion. This is useful if you want to compute the expansion for
`clausenc(x, s) - clausenc(0, s)`.

If `s` is wide, as determined by `iswide(s)` it computes a tighter
enclosure of the coefficients using a Taylor expansion in `s`.

The coefficient `C` is given by `gamma(1 - s) * sinpi(s / 2)` has a
singularity for positive integer values of `s` where the `gamma`
function blows up.

For even values this singularity is removable and in this case we can
still compute an enclosure of `C` by expanding at the integer,
including a remainder term, and explicitly cancel the removable
singularity. To do this we rewrite it using the reciprocal gamma
function [`rgamma`](@ref) as `sin(s / 2) / rgamma(1 - s)`

For odd values of `s` the singularity is not removable, instead the
non-analytic term coincides with one of the analytic terms, with their
coefficients blowing up in different directions. This case is
currently not handled in any special way and it just returns an
indeterminate result. Note that it is not clear how to handle this
case if we want to allow `s` overlapping odd integers.
"""
function clausenc_expansion(x::Arb, s::Arb, M::Integer; skip_constant = false)
    unique, s_integer = unique_integer(s)

    # Non-analytic term
    if unique && s_integer > 0
        if iseven(s_integer)
            rgamma_expansion =
                taylor_with_remainder(s -> rgamma(1 - s), Arb(s_integer), s, degree = 2)

            sin_expansion =
                taylor_with_remainder(s -> sinpi(s / 2), Arb(s_integer), s, degree = 2)

            # Expansion of gamma(1 - s) * sin(s / 2) with remainder term
            expansion =
                div_with_remainder(sin_expansion << 1, rgamma_expansion << 1, s - s_integer)

            C = expansion(s - s_integer)
        else
            C = Arblib.indeterminate!(zero(s))
        end
    else
        if iswide(s)
            C = ArbExtras.enclosure_series(s -> gamma(1 - s) * sinpi(s / 2), s, degree = 2)
        else
            C = gamma(1 - s) * sinpi(s / 2)
        end
    end
    e = s - 1

    # Analytic terms
    P = ArbSeries(degree = 2M - 2, prec = precision(x))
    start = skip_constant ? 1 : 0
    for m = start:M-1
        if iswide(s)
            z = ArbExtras.enclosure_series(s -> zeta(s - 2m), s, degree = 2)

            if !isfinite(z)
                # In some cases, when s overlaps zero (but not
                # always), the above returns NaN but the evaluation
                # below works.
                z = zeta(s - 2m)
            end
        else
            z = zeta(s - 2m)
        end
        P[2m] = (-1)^m * z / factorial(2m)
    end

    # Error term
    E = clausenc_expansion_remainder(x, s, M)

    return (C, e, P, E)
end

"""
    _clausen_expansion_remainder_ps(β::Integer)

Helper function for computing a representation of `p_k` for `k = 0:β`

Here `p_k` is defined by the recurrence ``p_{k + 1}(s) =
ψ^{(0)}(s)p_k(s) + p_k'(s)`` and ``p_0 = 1``.

It returns an `OffsetVector` with indices `0:β`. Each element is a
dictionary which maps a list of exponents to a coefficient. More
precisely the element `[q0, q1, ..., qk] => c` in the dictionary
represents the term
```
c * ψ(s, 0)^q0 * ψ(s, 1)^q1 * ... *ψ(s, k)^qk
```
where `ψ(s, i)` here means the function ``ψ^{(i)}(s)``. Note that this
doesn't actually evaluate the `ψ` functions in any way, it only stores
the representation of the terms.
"""
function _clausen_expansion_remainder_ps(β::Integer)
    p_0 = OrderedDict(Int[] => 1)
    ps = OffsetVector([p_0], 0:0)

    for k = 1:β
        p_k = empty(ps[k-1])

        # The term ψ^(0) * p_k
        for (exponents, coefficient) in ps[k-1]
            if isempty(exponents)
                new_exponents = [1] # ψ^(0) * 1 = ψ^(0)
            else
                new_exponents = copy(exponents)
                new_exponents[1] += 1 # Multiply by ψ^(0)
            end
            p_k[new_exponents] = coefficient
        end

        # The term p_k'
        for (exponents, coefficient) in ps[k-1]
            for l = 1:length(exponents)
                if !iszero(exponents[l])
                    # Differentiate ψ^(l+1)^exponents[l]
                    new_exponents = copy(exponents)

                    new_exponents[l] -= 1
                    if l == length(exponents)
                        push!(new_exponents, 1)
                    else
                        new_exponents[l+1] += 1
                    end

                    p_k[new_exponents] =
                        get(p_k, new_exponents, 0) + coefficient * exponents[l]
                end
            end
        end

        push!(ps, p_k)
    end

    return ps
end


"""
    clausenc_expansion_remainder(x::Arb, s::Arb, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausenc(x, s)` at zero up to order `2M - 2`, meaning that the
remainder is of order `2M`.

It is based on the bound given in [`lemma_clausen_remainder`](@ref).
It returns a ball centered at zero with this bound as the radius. Note
that since the bound is increasing in `x` it is valid for all `y`
satisfying `abs(y) <= abs(x)`.
"""
function clausenc_expansion_remainder(x::Arb, s::Arb, M::Integer)
    pi = Arb(π)

    abs(x) < 2pi || throw(DomainError(x, "x must be less than 2π"))
    s >= 0 || throw(DomainError(M, "s must be positive, got s = $s"))
    2M >= s + 1 || throw(DomainError(M, "must have 2M >= s + 1, got s = $s"))

    return Arblib.add_error!(
        zero(x),
        2(2pi)^(1 + s - 2M) * abs(sinpi(s / 2)) * zeta(2M + 1 - s) / (4pi^2 - x^2),
    )
end

"""
    clausenc_expansion_remainder(x::Arb, s::Arb, β::Integer, M::Integer)

Compute an enclosure of the remainder term in the asymptotic expansion
of `clausenc(x, s, β)` at zero up to order `2M - 2`, meaning that the
remainder is of order `2M`.

It is based on the bound given in
[`lemma_clausen_derivative_remainder`](@ref). It returns a ball
centered at zero with this bound as the radius. Note that since the
bound is increasing in `x` it is valid for all `y` satisfying `abs(y)
<= abs(x)`.

It uses [`_clausen_expansion_remainder_ps`](@ref) for computing an
expression for `p[0], ..., p[β]`.
"""
function clausenc_expansion_remainder(x::Arb, s::Arb, β::Integer, M::Integer)
    β == 0 && return clausenc_expansion_remainder(x, s, M)

    abs(x) < 2pi || throw(DomainError(x, "x must be less than 2π"))
    s >= 0 || throw(DomainError(M, "s must be positive, got s = $s"))
    2M >= s + 1 || throw(DomainError(M, "must have 2M >= s + 1, got s = $s"))

    twopi = 2Arb(π)

    # Function for computing multinomial
    # Denominator is smaller than numerator, so overflow
    # would always be caught by factorial(β)
    multinomial(β, j1, j2, j3) =
        Arb(factorial(β) // (factorial(j1) * factorial(j2) * factorial(j3)))

    # Upper bound of absolute value of polygamma functions for
    # 1 <= k <= β
    polygamma_bounds = [abs(real(polygamma(Acb(k), Acb(1 + 2M - s)))) for k = 1:β]

    # Upper bounds of zeta function and derivatives up to β
    zeta_expansion = zeta(ArbSeries((1 + 2M - s, 1), degree = β))
    zeta_bounds(j3) = abs(zeta_expansion[j3]) * factorial(j3)

    # p[k] for 0 <= k <= β
    ps = _clausen_expansion_remainder_ps(β)

    res = zero(x)

    for j1 = 0:β
        for j2 = 0:β-j1
            j3 = β - j1 - j2

            # Compute upper bound of
            # sum(abs(p_j2(1 + 2m - s)) * (x / 2π)^2m for m = M:Inf)
            term = zero(x)
            for (exponents, coefficient) in ps[j2]
                q0 = Arb(isempty(exponents) ? 0 : exponents[1])

                # Upper bound of polygamma-factors
                polygamma_factor = one(x)
                for i = 2:length(exponents)
                    polygamma_factor *= polygamma_bounds[i-1]^exponents[i]
                end

                term +=
                    coefficient * polygamma_factor * twopi^(-2M) / 2^(q0 / 2) *
                    lerch_phi((x / twopi)^2, -q0 / 2, M + Arb(1 // 2))
            end

            # Multiply with bounds for remaining factors

            # Upper bound of part of f factor
            F = 2(log(twopi) + Arb(π) / 2)^j1 * twopi^(s - 1)

            term *= multinomial(β, j1, j2, j3) * F * zeta_bounds(j3)

            res += term
        end
    end

    return Arblib.add_error!(zero(res), res)
end

"""
    clausencmzeta(x, s)

Compute `clausenc(x, s) - zeta(s)` for `s > 1`.

Note that for `s > 1` we have `clausenc(0, s) = zeta(s)` so this is
`clausenc` normalized to be zero at `x = 0`. When `s <= 1` we no
longer have `clausenc(0, s) = zeta(s)` and instead `clausenc(0, s)`
has a singularity, the expression `clausenc(x, s) - zeta(s)` is still
well defined however.

Typically this method just calls [`clausenc`](@ref) and [`zeta`](@ref)
directly and gives no direct benefit. However, if `s > 1` and is a
wide ball it can better handle the cancellations between the two terms
by computing them together.

For `s > 1` it uses [`lemma_clausencmzeta_monotone`](@ref), stating
that the function is non-decreasing in `s`.

We could compute a better enclosure for `s <= 1` by treating the terms
together and expanding in `s`. However this is not relevant for our
use case.
"""
function clausencmzeta(x::Arb, s::Arb)
    if s > 1 && iswide(s)
        sₗ, sᵤ = getinterval(Arb, s)
        sₗ > 1 || return Arblib.indeterminate!(zero(x))
        res_lower = clausenc(x, sₗ) - zeta(sₗ)
        res_upper = clausenc(x, sᵤ) - zeta(sᵤ)
        return Arb((res_lower, res_upper))
    else
        return clausenc(x, s) - zeta(s)
    end
end

function clausencmzeta(x::ArbSeries, s::Arb)
    res = zero(x)
    x₀ = x[0]

    res[0] = clausencmzeta(x₀, s)
    for i = 1:Arblib.degree(x)
        if i % 2 == 0
            res[i] = (-1)^(i ÷ 2) * clausenc(x₀, s - i) / factorial(i)
        else
            res[i] = -(-1)^(i ÷ 2) * clausens(x₀, s - i) / factorial(i)
        end
    end

    # Compose the Taylor series for the result with that of the input
    x_tmp = copy(x)
    x_tmp[0] = 0

    return Arblib.compose(res, x_tmp)
end

clausencmzeta(x::ArbSeries, s) = clausenc(x, s) - zeta(Arb(s, prec = precision(x)))

clausencmzeta(x::Arb, s::ArbSeries) = clausenc(x, s) - zeta(s)

function clausencmzeta(x, s)
    x, s = promote(x, s)
    return clausenc(x, s) - zeta(s)
end

"""
    clausencmzeta(x, s, β::Integer)

Compute `clausenc(x, s, β) - clausenc(0, s, β)`. Notice that
`clausenc(0, s, β)` is given by `zeta(s)` differentiated `β` times
with respect to `s`.

This method just calls [`clausenc`](@ref) and [`zeta`](@ref) directly
and gives no direct benefit other than converting `x` and `s` to the
same type to begin with.
"""
function clausencmzeta(x::Union{Arb,ArbSeries}, s, β::Integer)
    if isone(β)
        return clausenc(x, s, β) - dzeta(Arb(s))
    else
        return clausenc(x, s, β) - zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end
end

function clausencmzeta(x, s, β::Integer)
    x, s = promote(x, s)
    if isone(β)
        return clausenc(x, s, β) - dzeta(s)
    else
        return clausenc(x, s, β) - zeta(ArbSeries((s, 1), degree = β))[β] * factorial(β)
    end
end

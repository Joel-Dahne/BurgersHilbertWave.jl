"""
    dzeta(s)

Compute the Zeta function differentiated once with respect to `s`.
"""
dzeta(s::Arb) = zeta(ArbSeries((s, 1)))[1]
dzeta(s) = convert(float(typeof(s)), dzeta(Arb(s)))

"""
    rgamma(x)

Compute the reciprocal gamma function, defined by `rgamma(x) = 1 /
gamma(x)`.
"""
rgamma(x::Arb) = Arblib.rgamma!(zero(x), x)
rgamma(x::ArbSeries) = Arblib.rgamma_series!(zero(x), x, length(x))

"""
    zeta_deflated(s, a)

Compute the deflated zeta function.

The deflated zeta function is defined by
```
zeta_deflated(s, a) = zeta(s, a) + 1 / (1 - s)
```

**IMPROVE:** This current doesn't handle `s` overlapping one but not
exactly one.
"""
function zeta_deflated(s::Arb, a::Arb)
    res = ArbSeries(s)
    Arblib.zeta_series!(res, res, a, 1, length(res))
    return res[0]
end

function zeta_deflated(s::ArbSeries, a::Arb)
    return Arblib.zeta_series!(zero(s), s, a, 1, length(s))
end

"""
    lerch_phi(z, s, a)

Compute a naive enclosure of the Lerch transendent
```
lerch_phi(z, s, a) = sum(z^n / (n + a)^s for n = 0:Inf)
```

It only supports `abs(z) < 1`, `s < 0` and `a > 0`. We treat the whole
series as a tail, for bounding it see
[https://fredrikj.net/blog/2022/02/computing-the-lerch-transcendent/](Fredrik
Johanssons blog).

Note that this is a very naive implementation and only intended for
use in [`clausenc_expansion_remainder`](@ref). Once the version of Arb
that implements this function is released we should switch to using
that one.
"""
function lerch_phi(z::Arb, s::Arb, a::Arb)
    @assert abs(z) < 1
    @assert Arblib.isnonpositive(s)
    @assert Arblib.ispositive(a)

    # Find C
    C = exp(-s / a)

    abs(C * z) < 1 || throw(ArgumentError("C to large to satisfy C * z < 1"))

    return a^-s / (1 - C * abs(z))
end

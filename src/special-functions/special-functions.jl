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
    lerch_phi(z::Arb, s::Arb, a::Arb)

Compute the Lerch transcendent
```
lerch_phi(z, s, a) = sum(z^n / (n + a)^s for n = 0:Inf)
```
Returns an indeterminate value if the result is not real.
"""
function lerch_phi(z::Arb, s::Arb, a::Arb)
    res = Arblib.dirichlet_lerch_phi!(Acb(prec = precision(z)), Acb(z), Acb(s), Acb(a))

    if isreal(res)
        return real(res)
    else
        return indeterminate(z)
    end
end

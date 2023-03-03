"""
    lemma_asymptotic_expansion

Corresponds to Lemma 4.3 in the paper.

# Statement

The approximation `u0` has the asymptotic expansions
```
u0(x) = -1 / π * abs(x) * log(abs(x)) - (γ - 1) / π * abs(x) +
    sum(a⁰[j] * abspow(x, -α + j * p0) for j = 1:N0) +
    sum(m = 1:Inf) do m
        (-1)^m / factorial(2m) * (
            2 / π^2 * dzeta(2 - 2m)
            + sum(a[j] * zeta(1 - α + j * p0 - 2m) for j = 1:N0)
            + sum(b[n] * n^2m for n = 1:N1)
        ) * x^2m
    end
```
and
```
H(u0)(x) = -1 / 2π^2 * x^2 * log(abs(x))^2 + (3 - 2γ) / 2π^2 * x^2 * log(abs(x)) -
    sum(A⁰[j] * abspow(x, 1 - α + j * p0) for j = 1:N0) +
    1 / 2 *(
        (36γ - 12γ^2 - 24γ₁ - 42 + π^2) / 12π^2
        + sum(a[j] * zeta(-α + j * p0) for j = 1:N0)
        + sum(b[n] * n for n = 1:N1)
    ) * x^2 -
    sum(m = 2:Inf) do m
        (-1)^m / factorial(2m) * (
            2 / π^2 * dzeta(3 - 2m)
            + sum(a[j] * zeta(2 - α + j * p0 - 2m) for j = 1:N0)
            + sum(b[n] * n^(2m - 1) for n = 1:N1)
        ) * x^2m
    end
```
valid for `abs(x) < 2π`. The following notation is used
```
p0 = u0.p0
N0 = u0.N0
N1 = u0.N1
a[j] = u0.a[j]
b[n] = u0.b[n]

a⁰[j] = gamma(α - j * p0) * cospi((α - j * p0) / 2) * a[j]
A⁰[j] = gamma(α - 1 - j * p0) * cospi((α - 1 - j * p0) / 2) * a[j]
```
`γ₁` is the Stieltjes constant and `γ = γ₀` is the Euler's constant.
"""
function lemma_asymptotic_expansion end

"""
    lemma_inv_u0_normalised

Corresponds to Lemma 5.1 in the paper.

# Statement
The function
```
-abs(x) * log(abs(x)) / u0(x)
```
is positive and bounded at `x = 0` and for `abs(x) < 1` it has the
expansion
```
abs(x) * log(abs(x)) / u0(x) = -inv(
    -1 / π
    -(γ - 1) / π / log(abs(x))
    + sum(a⁰[j] * abspow(x, -1 - α + j * p0) / log(abs(x)) for j = 1:N0)
    + sum(m = 1:Inf) do m
        (-1)^m / factorial(2m) * (
            2 / π^2 * dzeta(2 - 2m)
            + sum(a[j] * zeta(1 - α + j * p0 - 2m) for j = 1:N0)
            + sum(b[n] * n^2m for n = 1:N1)
        ) * abs(x)^(2m - 1) / log(abs(x))
    end
)
```
"""
function lemma_inv_u0_normalised end

"""
    lemma_defect_normalised

Corresponds to Lemma 5.2 in the paper.

# Statement
The function
```
(H(u0)(x) + 1 / 2 * u0(x)^2) / (x^2 * log(abs(x)))
```
is non-zero and bounded at `x = 0` and for `abs(x) < 1`.

The lemma also gives the expansion at `x = 0`. We do however not use
this expression of the expansion directly. See the asymptotic version
of [`defect`](@ref) for [`BHAnsatz`](@ref) for how the expansion is
computed.
"""
lemma_defect_normalised

"""
    lemma_I_hat

Corresponds to Lemma 6.1 in the paper.

# Statement
For all `0 < x < π` the function
```
I_hat(x, t) = log(sin(abs(x * (1 - t) / 2)) * sin(x * (1 + t) / 2) / sin(x * t / 2)^2)
```
is decreasing and continuous in `t` for `0 < t < 1` and has the limit
`Inf` for `t` tending to `0` and `-Inf` for `t` tending to `1`.
Moreover, the unique root, `r(x)`, is decreasing in `x` and satisfies
the inequality
```
1 / 2 < r(x) < 1 / sqrt(2)
```
"""
function lemma_I_hat end

"""
    lemma_I

Corresponds to Lemma 6.2 in the paper.

# Statement
For `0 < x < π` we have
```
I(x, y) = log(sin(abs(x - y) / 2) * sin((x + y) / 2) / sin(y / 2)^2) < 0
```
for all `x < y < π`.
"""

"""
    lemma_U1_asymptotic

Corresponds to Lemma 6.3 in the paper.

# Statement
Let
```
U1(x) = x^2 * ∫ abs(I_hat(x, t)) * t * sqrt(log(1 + inv(x * t))) dt
```
where the integration is from `0` to `1` and `I_hat(x, t)` is as in
[`lemma_I_hat`](@ref).

For `0 <= x <= ϵ` with `ϵ < 1` we have
```
U1(x) / (-x^2 * log(x) * sqrt(log(1 + inv(x)))) <= 1 / sqrt(log(1 + inv(x))) * (
    log(2) / sqrt(log(inv(x)))
    + (c1 + log(2) * sqrt(log(1 + x))) / log(inv(x))
    + 3R1 / 8 * x^2 * (2 / sqrt(log(inv(x))) + (sqrt(π / 2) + 2sqrt(log(1 + x))) / log(inv(x)))
)
```
where
```
c1 = ∫ abs(log(1 / t^2 - 1)) * t * sqrt(log(inv(t))) dt
```
integrated from `0` to `1` and `R1` is the supremum for `0 <= y <= ϵ`
of
```
1 / 2 * abs(d²/dy²(log(sin(y) / y)))
```
where we use `d²/dy²` to denote differentiating twice w.r.t. `y`.
"""
function lemma_U1_asymptotic end

"""
    lemma_U2_asymptotic

Corresponds to Lemma 6.4 in the paper.

# Statement
Let
```
U2(x) = ∫ abs(I(x, y)) * y * sqrt(log(1 + inv(y))) dy
```
where the integration is from `x` to `π` and `I(x, t)` is as in
[`lemma_I`](@ref).

For `0 <= x <= ϵ` with `ϵ < 1 / 2` we have
```
U2(x) / (-x^2 * log(x) * sqrt(log(1 + inv(x)))) <= log(16 / 3sqrt(3)) / log(inv(x)) +
    (
        2 / 3 * log(inv(2x))^(3 / 2) / (log(inv(x)) * sqrt(log(1 + inv(x))))
        + R2 * sqrt(log(inv(2x))) / (8log(inv(x)) * sqrt(log(1 + inv(x))))
        + sqrt(log(2)) / sqrt(log(inv(1 + inv(x))))
        - log(2)^(3 / 2) / (log(inv(x)) * sqrt(log(1 + inv(x))))
        + R2 * sqrt(log(2)) * (1 - 4x^2) / (8log(inv(x)) * sqrt(log(1 + inv(x))))
    ) + sqrt(log(2)) / (log(inv(x)) * sqrt(log(1 + inv(x)))) * (
        1 / 2 * log((π^2 - x^2) / (1 - x^2))
        + log(1 - x^2) / 2x^2
        - π^2 * log(1 - x^2 / π^2) / 2x^2
    ) + D1 * c2 / (log(inv(x)) * sqrt(log(1 + inv(x))))
```
where
```
c2 = ∫ y * sqrt(log(1 + inv(y))) dt
```
integrated from `0` to `π`, `R2` is the supremum for `0 <= y <= 1 / 4`
of
```
1 / 2 * abs(d²/dy²(log(1 - y)))
```
where we use `d²/dy²` to denote differentiating twice w.r.t. `y` and
`D1` is the supremum for `0 < x < ϵ` of
```
-(log(sinc((1 - x / π) / 2)) + log(sinc((1 + x / π) / 2)) - 2log(sinc(1 / 2))) / x^2
```
Note that we are here using the Julia convention that `sinc(x) = sin(π
* x) / (π * x)`, which is different from the one used in the paper.
"""
function lemma_U2_asymptotic end

"""
    lemma_U2_term_bound

Corresponds to Lemma 6.5 in the paper.

# Statement
The function
```
f(x) = 2 / 3 * log(inv(2x))^(3 / 2) / (log(inv(x)) * sqrt(log(1 + inv(x))))
```
is bounded from above by `2 / 3` and is decreasing in `x` for `0 < x <
1 / 2`.
"""
function lemma_U2_term_bound end

"""
    lemma_bound_delta0

Corresponds to Lemma 7.2 in the paper.

We don't use the statement of the lemma in the code but some of the
details of the asymptotic version of [`F0`](@ref) are discussed in the
proof of the lemma.
"""
function lemma_bound_delta0 end

"""
    lemma_removable_singularities

Corresponds to Lemma A.1 in the paper.

# Statement
Let `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.

Let `m` be a non-negative integer and let `I` be an interval
containing zero. Consider a function `f(x)` with a zero of order `n`
at `x = 0` and such that `d(f, m + n)` is absolutely continuous on
`I`. Then for all `x ∈ I` we have
```
f(x) / x^n = sum(f[k + n](0) * x^k for k = 0:m) + f[m + n + 1](ξ) * x^(m + 1)
```
for some `ξ` between `0` and `x`. Here we use `f[k](x) = d(f, k)(x) /
factorial(k)`.

Furthermore, if `d(f, m + n + p)` is absolutely continuous for some
non-negative integer `p` we have
```
d(x -> f(x) / x^n, p)(x) =
    sum(factorial(k + p) / factorial(k) * f[k + n + p](0) * x^k for k = 0:m) +
    factorial(m + p + 1) / factorial(m + 1) * f[m + n + p + 1](ξ) * x^(m + 1)
```
for some `ξ` between `0` and `x`.
"""
function lemma_removable_singularities end

"""
    lemma_clausenc_monotone

Corresponds to Lemma B.1 in the paper.

# Statement
For all real `s` the Clausen function `clausenc(x, s)` is monotone in
`x` for `0 < x < π`. For `s > 0` it is non-increasing.
"""
function lemma_clausenc_monotone end

"""
    lemma_clausencmzeta_monotone

Corresponds to Lemma B.2 in the paper.

# Statement
For `s > 1` and `x` real the function `clausencmzeta(x, s)` is
non-decreasing in `s`.
"""
function lemma_clausencmzeta_monotone end

"""
    lemma_clausen_remainder

Corresponds to Lemma B.3 in the paper.

# Statement
Let `s >= 0, 2M >= s + 1` and `abs(x) < 2π`, we then have the
following bounds for the tails in
[`equation_clausenc_expansion`](@ref) and
[`equation_clausens_expansion`](@ref)
```
abs(sum((-1)^m * zeta(s - 2m) * x^2m / factorial(2m) for m = M:Inf)) <=
    2(2π)^(1 + s - 2M)abs(sinpi(s / 2)) * zeta(2M + 1 - s) * x^2M / (4π^2 - x^2)

abs(sum((-1)^m * zeta(s - 2m - 1) * x^(2m + 1) / factorial(2m + 1) for m = M:Inf)) <=
    2(2π)^(s - 2M)abs(cospi(s / 2)) * zeta(2M + 2 - s) * x^(2M + 1) / (4π^2 - x^2)
```
"""
function lemma_clausen_remainder end

"""
    lemma_clausen_derivative_remainder

Corresponds to Lemma B.4 in the paper.

# Statement
Let `d(f, k)(x)` denote the `k`th derivative of `f` evaluated at `x`.

Let `β >= 1, s >= 0, 2M >= s + 1` and `abs(x) < 2π`, we then have the
following bounds for the tails in ... and ...
```
abs(sum((-1)^m * d(zeta, β)(s - 2m) * x^2m / factorial(2m) for m = M:Inf)) <=
    sum(j1 + j2 + j3 = β) do j1, j2, j3
        binomial(β, (j1, j2, j3)) *
        2(log(2π) + π / 2)^j1 *
        (2π)^(s - 1) *
        abs(d(zeta, j3)(1 + 2M - s)) *
        sum(p[j2](1 + 2m - s) * (x / 2π)^2m for m = M:Inf)
    end

abs(sum((-1)^m * d(zeta, β)(s - 2m - 1) * x^(2m + 1) / factorial(2m + 1) for m = M:Inf)) <=
    sum(j1 + j2 + j3 = β) do j1, j2, j3
        binomial(β, (j1, j2, j3)) *
        2(log(2π) + π / 2)^j1 *
        (2π)^(s - 1) *
        abs(d(zeta, j3)(2 + 2M - s)) *
        sum(p[j2](2 + 2m - s) * (x / 2π)^(2m + 1) for m = M:Inf)
    end
```

Here `p[j2]` is recursively given by
```
p[k + 1](s) = polygamma(0, s) * p[k](s) + d(p[k], 1)(s)
p[0] = 1
```
It is given by a linear combination of terms of the form
```
polygamma(0, s)^q[0] * polygamma(1, s)^q[1] * ⋯ * polygamma(j2 - 1, s)^q[j2 - 1]
```
We have the following bounds
```
sum(
    polygamma(0, 1 + 2m - s)^q[0] * ⋯ * polygamma(j2 - 1, 1 + 2m - s)^q[j2 - 1] (x / 2π)^2m
    for m = M:Inf
) <= abs(polygamma(1, 1 + 2M - s)^q[1] * ⋯ * polygamma(j2 - 1, 1 + 2M - s)^q[j2 - 1]) *
    1 / 2^(q[0] / 2) * (2π)^(-2M) * lerch_phi(x^2 / 4π^2, -q[0] / 2, M + 1) * x^2M

sum(
    polygamma(0, 2 + 2m - s)^q[0] * ⋯ * polygamma(j2 - 1, 2 + 2m - s)^q[j2 - 1] (x / 2π)^(2m + 1)
    for m = M:Inf
) <= abs(polygamma(1, 2 + 2M - s)^q[1] * ⋯ * polygamma(j2 - 1, 2 + 2M - s)^q[j2 - 1]) *
    1 / 2^(q[0] / 2) * (2π)^(-2M - 1) * lerch_phi(x^2 / 4π^2, -q[0] / 2, M + 1) * x^(2M + 1)
```
"""
function lemma_clausen_derivative_remainder end

"""
    lemma_U_singular_integrals

Corresponds to Lemma C.1 in the paper.

# Statement
Let
```
U122(x) = -∫ I_hat(x, t) * t * sqrt(log(1 + inv(t))) dt
```
integrated from `1 - δ1` to `1` and
```
U21(x) = -∫ I(x, y) * y * sqrt(log(1 + inv(y))) dy
```
integrated from `x` to `x + δ2`.

We have
```
U122(x) = ξ1 / x * (
    (-clausens(0, 2) + clausens(2x, 2) - 2clausens(x, 2)) -
    (-clausens(x * δ1, s) + clausens(x * (2 - δ1), s) - 2clausens(x * (1 - δ1), s))
)
```
for some `ξ1` in the image of `t * sqrt(1 + inv(x * t))` for `1 - δ1
<= t <= 1` Furthermore
```
U21(x) = ξ2 * (
    (clausens(δ2, 2) + clausens(2x + δ2, 2) - 2clausens(x + δ2, 2)) -
    (-clausens(0, s) + clausens(2x, s) - 2clausens(x, s))
)
```
for some `ξ2` in the image of `y * sqrt(1 + inv(y))` for `x <= y <= x
+ δ2`.
"""
function lemma_U_singular_integrals end

"""
    lemma_U11_bound_zero

Corresponds to Lemma C.2 in the paper.

# Statement
For `0 < x < π` and `0 < t < 1 / 10x` the function
```
log(sin(x * t / 2)) * t * sqrt(log(1 + inv(x * t)))
```
is decreasing in `t`.
"""
function lemma_U11_bound_zero end

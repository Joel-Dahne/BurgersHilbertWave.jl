# This file contains code for evaluation of the approximate solution
# in different ways

function (u0::FractionalKdVAnsatz)(x, ::Ball)
    res = zero(u0.α)

    for j = 0:u0.N0
        s = 1 - u0.α + j * u0.p0
        res += u0.a[j] * clausencmzeta(x, s)
    end

    for n = 1:u0.N1
        res += u0.b[n] * (cos(n * x) - 1)
    end

    return res
end

function H(u0::FractionalKdVAnsatz, ::Ball)
    return x -> begin
        res = zero(u0.α)

        for j = 0:u0.N0
            s = 1 - 2u0.α + j * u0.p0
            res -= u0.a[j] * clausencmzeta(x, s)
        end

        for n = 1:u0.N1
            res -= u0.b[n] * n^u0.α * (cos(n * x) - 1)
        end

        return res
    end
end

"""
    D(u0::FractionalKdVAnsatz, xs::AbstractVector)
Returns a function such that D(u0, xs)(a, b) computes D(u0)(x) on the
points x ∈ xs with u0.a and u0.b set to the given values. Does this in
an efficient way by precomputing as much as possible.
"""
function D(u0::FractionalKdVAnsatz, xs::AbstractVector)
    u0_xs_a_precomputed = zeros(length(xs), u0.N0 + 1)
    u0_xs_b_precomputed = zeros(length(xs), u0.N1)
    Hu0_xs_a_precomputed = zeros(length(xs), u0.N0 + 1)
    Hu0_xs_b_precomputed = zeros(length(xs), u0.N1)

    for i in eachindex(xs)
        x = xs[i]
        for j = 0:u0.N0
            u0_xs_a_precomputed[i, j+1] = clausencmzeta(x, 1 - u0.α + j * u0.p0)
            Hu0_xs_a_precomputed[i, j+1] = -clausencmzeta(x, 1 - 2u0.α + j * u0.p0)
        end
        for n = 1:u0.N1
            u0_xs_b_precomputed[i, n] = cos(n * x) - 1
            Hu0_xs_b_precomputed[i, n] = -n^u0.α * (cos(n * x) - 1)
        end
    end

    return (a, b) -> begin
        return (
            (u0_xs_a_precomputed * a .+ u0_xs_b_precomputed * b) .^ 2 ./ 2 .+
            (Hu0_xs_a_precomputed * a .+ Hu0_xs_b_precomputed * b)
        )
    end
end

"""
    D(u0::FractionalKdVAnsatz, evaltype::Symbolic; M::Integer = 5)

Return a function such that `D(u0, evaltype, N)(a)` computes the
coefficients in the asymptotic expansion with indices `3` to `u0.N0 +
1` using the values from `a`.

This is used in [`_findas`](@ref) for numerically finding values for
`a`.
"""
function D(u0::FractionalKdVAnsatz{T}, ::Symbolic; M::Integer = 5) where {T}
    # Given a key get its exponent
    key_exponent = ((i, j, m),) -> -i * u0.α + j * u0.p0 + m

    # Precompute for u0
    u0_precomputed = OrderedDict{NTuple{3,Int},OrderedDict{Int,T}}()
    for j = 0:u0.N0
        s = u0.α - j * u0.p0
        u0_precomputed[(1, j, 0)] = OrderedDict(j => gamma(s) * cospi(s / 2))
    end

    for m = 1:M-1
        u0_precomputed[(0, 0, 2m)] = OrderedDict(
            j => (-1)^m * zeta(1 - u0.α + j * u0.p0 - 2m) / factorial(2m) for j = 0:u0.N0
        )
    end

    # Precompute H(u0)
    Hu0_precomputed = OrderedDict{NTuple{3,Int},OrderedDict{Int,T}}()

    for j = 0:u0.N0
        s = 2u0.α - j * u0.p0
        if s == -1
            # -gamma(s) * cospi(s / 2) has a removable singularity at
            # -s = -1, where it takes the value π / 2
            Hu0_precomputed[(2, j, 0)] = OrderedDict(j => π / 2)
        else
            Hu0_precomputed[(2, j, 0)] = OrderedDict(j => -gamma(s) * cospi(s / 2))
        end
    end

    for m = 1:M-1
        Hu0_precomputed[(0, 0, 2m)] = OrderedDict(
            j => -(-1)^m * zeta(1 - 2u0.α + j * u0.p0 - 2m) / factorial(2m) for j = 0:u0.N0
        )
    end

    # We want to find the value of the maximum exponent we return. To
    # do this we first compute ALL keys we will encounter and then
    # sort and find the value for the maximum one we return.
    all_keys = OrderedDict{NTuple{3,Int},T}()
    for key in keys(Hu0_precomputed)
        if !haskey(all_keys, key)
            all_keys[key] = key_exponent(key)
        end
    end
    for key1 in keys(u0_precomputed)
        for key2 in keys(u0_precomputed)
            key = key1 .+ key2
            if !haskey(all_keys, key)
                all_keys[key] = key_exponent(key)
            end
        end
    end
    # Find the keys we care about
    returned_keys = collect(sort(all_keys, by = Float64 ∘ key_exponent))[3:u0.N0+2]
    # Find the largest key/exponent we care about
    maxkey, maxexponent = returned_keys[end]

    # Check that M is large enough
    @assert maxexponent < 2M

    # Filter out any keys larger than the max exponent
    filter!(keyvalue -> !(key_exponent(keyvalue[1]) > maxexponent), u0_precomputed)
    filter!(keyvalue -> !(key_exponent(keyvalue[1]) > maxexponent), Hu0_precomputed)

    # Sort the dictionaries by exponent
    sort!(u0_precomputed, by = Float64 ∘ key_exponent)
    sort!(Hu0_precomputed, by = Float64 ∘ key_exponent)

    # Function to compute the dictionaries u0_res and H0_res
    sum_dict(precomputed, a, S) = begin
        res = empty(precomputed, S)
        @inbounds for (key, dict) in precomputed
            for (j, v) in dict
                res[key] = get(res, key, zero(S)) + v * a[j]
            end
        end
        return res
    end

    return a -> begin
        S = promote_type(T, eltype(a))

        u0_res = sum_dict(u0_precomputed, a, S)
        Hu0_res = sum_dict(Hu0_precomputed, a, S)

        # Compute u0^2/2
        res = empty(u0_res)
        u0_res_vector = collect(u0_res)
        @inbounds for (i, (key1, y1)) in enumerate(u0_res_vector)
            key = 2 .* key1
            key_exponent(key) > maxexponent && continue
            res[key] = get(res, key, zero(S)) + y1^2 / 2
            for j = i+1:length(u0_res_vector)
                (key2, y2) = u0_res_vector[j]
                key = key1 .+ key2
                key_exponent(key) > maxexponent && break
                res[key] = get(res, key, zero(S)) + y1 * y2
            end
        end

        # Compute u0^2/2 + H(u0)
        merge!(+, res, Hu0_res)

        return collect(values(sort(res, by = Float64 ∘ key_exponent)))[3:u0.N0+2]
    end
end

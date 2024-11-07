module InitialData

using ForwardDiff
using LambertW

function rstar(r::T, M::T) where {T}
    return r + 2M * log(abs(r / (2M) - one(T)))
end

function r(rstar::T, M::T) where {T}
    return 2(M + M * lambertw(exp(-one(T) + rstar / (2M))))
end

@inline function Gaussian1D(t::Real, r::Real, σ::Real, r0::Real)
    return exp(-(r - r0 + t)^2 / σ^2)
end

@inline function dtGaussian1D(t::Real, r::Real, σ::Real, r0::Real)
    return ForwardDiff.derivative(t1 -> Gaussian1D(t1, r, σ, r0), t)
end

@inline function drGaussian1D(t::Real, r::Real, σ::Real, r0::Real)
    return ForwardDiff.derivative(r1 -> Gaussian1D(t, r1, σ, r0), r)
end

#creating general functions
@inline function Lamb(l::T, r::T, M::T) where {T}
    return (l - one(T)) * (l - 2one(T)) - 6M / r
end
@inline function f(r::T, M::T) where {T}
    return one(T) - 2M / r
end
@inline function mul(l::T) where {T}
    return (l - one(T)) * (l - 2one(T))
end

#defining the function
@inline function H2(t::Real, r::Real, M::Real, P::Real, L::Real)
    return t * r * M
end
@inline function G(t::Real, r::Real, M::Real, P::Real, L::Real)
    return t * r * M
end
@inline function K(t::Real, r::Real, M::Real, P::Real, L::Real)
    return t * r * M
end
@inline function h1(t::Real, r::Real, M::Real, P::Real, L::Real)
    return t * r * M
end

#defining the relevant derivatives for the gauge invariant stuff
@inline function drG(t::Real, r::Real, M::Real, P::Real, L::Real)
    return ForwardDiff.derivative(r1 -> G(t, r1, M, P, L), r)
end

@inline function drdrG(t::Real, r::Real, M::Real, P::Real, L::Real)
    return ForwardDiff.derivative(r1 -> drGx(t, r1, M, P, L), r)
end

@inline function dtG(t::Real, r::Real, M::Real, P::Real, L::Real)
    return ForwardDiff.derivative(t1 -> G(t1, r, M, P, L), r)
end

@inline function drh1(t::Real, r::Real, M::Real, P::Real, L::Real)
    return ForwardDiff.derivative(r1 -> h1(t, r1, M, P, L), r)
end

#I need to write down the christoffel symbols to write down the covariant derivative
@inline function Christ(params, r::Real, M::Real)
    if params == ["r", "r", "r"]
        return -(M / r^2) / (1 - 2 * M / r)
    end
    if params == ["r", "t", "t"]
        return (M / r^2) * (1 - 2 * M / r)
    end
    if params == ["t", "r", "t"]
        return (M / r^2) / (1 - 2 * M / r)
    end
    if params == ["t", "t", "r"]
        return (M / r^2) / (1 - 2 * M / r)
    else
        return 0
    end
end

#here i define the gauge invariant quantities
@inline function khat(l::Real, t::Real, r::Real, M::Real, P::Real, L::Real)
    return K(t, r, M, P, L) + l * (l + 1) / 2 * G(t, r, M, P, L) -
           (2 / r) * h1(t, r, M, P, L) + r * +r * drG(t, r, M, P, L)
end

@inline function fhatrr(t::Real, r::Real, M::Real, P::Real, L::Real)
    return H2(t, r, M, P, L) - 2 * drh1(t, r, M, P, L) + 4 * drG(t, r, M, P, L) +
           r^2 *
           (drdrG(t, r, M, P, L) - Christ(["t", "r", "r"], r, M) * dtG(t, r, M, P, L) -
            Christ(["r", "r", "r"], r, M) * drG(t, r, M, P, L))
end

@inline function drKhat(l::Real, t::Real, r::Real, M::Real, P::Real, L::Real)
    return ForwardDiff.derivative(r1 -> khat(l, t, r1, M, P, L), r)
end

#variable for the extrensic curvature for l=2
@inline function KK(r::Real, M::Real, P::Real, L::Real)
    return P * L / r^3 - P * L^3 * 96 / (r^3 * 7 * (sqrt(r - 2 * M) + sqrt(r))^4)
end

@inline function KG(r::Real, M::Real, P::Real, L::Real)
    return P * L / r^3 - P * L^3 * 32 / (r^3 * 21 * (sqrt(r - 2 * M) + sqrt(r))^4)
end
@inline function Kh1(r::Real, M::Real, P::Real, L::Real)
    return -P * L^3 * 32 / (r^(3 / 2) * 7 * sqrt(r - 2 * M) * (sqrt(r - 2 * M) + sqrt(r))^4)
end
@inline function KH2(r::Real, M::Real, P::Real, L::Real)
    return -4 * P * L / r^3 - P * L^3 * 96 / (r^3 * 7 * (sqrt(r - 2 * M) + sqrt(r))^4)
end
#l , t, r, M, P, L

#writing down the Zerilli functions
@inline function InitialPsiValues(t::Real, r::Real, M::Real, P::Real, L::Real, l::Real)
    return 2 * r / (l * (l + 1)) * (khat(l, t, r, M, P, L) +
                                    2 / (Lamb(l, r, M)) * f(r, M)^2 * fhatrr(t, r, M, P, L) -
                                    r * f(r, M) * drKhat(l, t, r, M, P, L))
end

@inline function drInitialPsiValues(t::Real, r::Real, M::Real, P::Real, L::Real)
    return ForwardDiff.derivative(r1 -> InitialPsiValues(l, t, r1, M, P, L), r)
end

@inline function dtInitialPsiValues(t::Real, r::Real, M::Real, P::Real, L::Real)
    return t * r * M * P * L
end

# Potential for the function
# THIS IS EVALUATED at r and not rstar
@inline function Vpot(l::Real, r::Real, M::Real)
    return f(r, M) / (Lamb(l, r, M)^2) * (mul(l)^2 * ((mul(l) + 2) / r^2 + 6 * M / r^3) +
                                          36 * M^2 / r^4 * (mul(l) + 2 * M / r))
end

@inline function poschl_teller(l::Real, r::Real)
    # return -(l * (l + 1) / 2) * sech(r)^2
    return -sech(r)^2
end

end #end of module

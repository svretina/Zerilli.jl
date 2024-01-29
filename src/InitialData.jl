module InitialData

using ForwardDiff

@inline function Gaussian1D(t::Real, r::Real, σ::Real, r0::Real)
    return exp(-(r - r0 + t)^2 / σ^2)
end

@inline function dtGaussian1D(t::Real, r::Real, σ::Real, r0::Real)
    return ForwardDiff.derivative(t1 -> Gaussian1D(t1, r, σ, r0), t)
end

@inline function drGaussian1D(t::Real, r::Real, σ::Real, r0::Real)
    return ForwardDiff.derivative(r1 -> Gaussian1D(t, r1, σ, r0), r)
end

end #end of module

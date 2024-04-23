module Run

using ..InitialData
using ..Integrator
using ..ODE
using GridFunctions

export zerilli


function rstar(r, M)
    return r + 2M * log(abs(r - 2M))
end

# you need to provide the transformations
function r(rstar, M) end

function zerilli(ncells::Integer, tf::Real, cfl::Real;
    boundary_type::Symbol=:radiative, folder="", save_every=1)

    @assert boundary_type === :radiative || boundary_type === :reflective
    M = 1.0
    horizon = 2M
    grid = UniformGrid([horizon, 20M + horizon], ncells)

    N = ncells + 1

    nt = ceil(Int64, tf / (cfl * spacing(grid)))
    t = UniformGrid([0.0, tf], nt)

    # Initial Data
    σ = 2
    r0 = 20M / 2 + horizon
    Φ = InitialData.Gaussian1D.(0, coords(grid), σ, r0)
    Π = InitialData.dtGaussian1D.(0, coords(grid), σ, r0)
    Ψ = InitialData.drGaussian1D.(0, coords(grid), σ, r0)

    l = 0.0
    params = (h=spacing(grid), N=N, bc=boundary_type, t=t, ti=coords(t),
        dt=spacing(t), save_every=save_every, M=M, grid=grid,
        folder=folder, cfl=cfl, rcoord=coords(grid), l=l)
    statevector = hcat(Φ, Π, Ψ)
    Integrator.solve(ODE.rhs!, statevector, params)
    return nothing
end

end # end of module

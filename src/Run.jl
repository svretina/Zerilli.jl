module Run

using ..InitialData
using ..Integrator
using ..ODE
using GridFunctions

export zerilli

function zerilli(ncells::Integer, tf::Real, cfl::Real,
    boundary_type::Symbol=:radiative, folder="", overwrite=true, save_every=1)

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

    params = (h=spacing(grid), N=N, bc=boundary_type, ti=coords(t), dt=spacing(t))
    statevector = hcat(Φ, Π, Ψ)
    Integrator.solve(ODE.rhs!, statevector, params, t, grid, M, save_every, folder, overwrite)
    return nothing
end

end # end of module

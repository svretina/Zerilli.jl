module Run

using ..InitialData
using ..Integrator
using ..ODE
using GridFunctions

export zerilli

function zerilli(domain::Vector, ncells::Integer,
    tf::Real, cfl::Real, M::Real=1.0, l::Real=0;
    boundary_type::Symbol=:radiative,
    folder="", save_every=1)

    @assert boundary_type === :radiative || boundary_type === :reflective

    grid = UniformGrid(domain, ncells)

    N = ncells + 1
    nt = ceil(Int64, tf / (cfl * spacing(grid)))
    t = UniformGrid([0.0, tf], nt)

    # Initial Data
    σ = 6
    r0 = 20# 20M / 2 + horizon
    # Φ = InitialData.Gaussian1D.(0, coords(grid), σ, r0)
    # Π = InitialData.dtGaussian1D.(0, coords(grid), σ, r0)
    # Ψ = InitialData.drGaussian1D.(0, coords(grid), σ, r0)

    Φ = rand(N)
    Π = rand(N)
    Ψ = rand(N)

    params = (h=spacing(grid), N=N, bc=boundary_type, t=t, ti=coords(t),
        dt=spacing(t), save_every=save_every, M=M, grid=grid,
        folder=folder, cfl=cfl, rcoord=coords(grid), l=l)
    statevector = hcat(Φ, Π, Ψ)
    Integrator.solve(ODE.rhs!, statevector, params)
    return nothing
end

end # end of module

module Run

using ..InitialData
# using ..Integrator
using ..ODE
using GridFunctions
using OrdinaryDiffEqTsit5

export zerilli

function zerilli(; domain::Vector, ncells::Integer,
                 tf::Real, cfl::Real, M::Real=1.0, l::Real=0,
                 boundary_type::Symbol=:radiative)
    @assert boundary_type === :radiative || boundary_type === :reflective

    grid = UniformGrid(domain, ncells)

    N = ncells + 1
    nt = ceil(Int64, tf / (cfl * spacing(grid)))
    t = UniformGrid([0.0, tf], nt)
    tspan = (0.0, tf)
    dt = cfl * spacing(grid)

    # Initial Data

    σ = 6
    r0 = 0 # 20# 20M / 2 + horizon
    Φ = InitialData.Gaussian1D.(0.0, coords(grid), σ, r0)
    Π = InitialData.dtGaussian1D.(0.0, coords(grid), σ, r0)
    Ψ = InitialData.drGaussian1D.(0.0, coords(grid), σ, r0)

    # Φ = rand(N)
    # Π = rand(N)
    # Ψ = rand(N)

    params = (h=spacing(grid), N=N, bc=boundary_type, t=t, ti=coords(t),
              dt=spacing(t), M=M, grid=grid, L=domain[2],
              cfl=cfl, rcoord=coords(grid), l=l)
    statevector = hcat(Φ, Π, Ψ)

    ode = ODEProblem{true}(ODE.rhs!, statevector, tspan, params)

    sol = solve(ode,
                Tsit5();
                adaptive=false,
                dt=dt, calck=false)
    return (t, grid, sol, params)
end

end # end of module

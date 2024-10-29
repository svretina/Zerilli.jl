module ODE

using ..InitialData

@inline function rhs!(du::AbstractArray{<:Real}, U::AbstractArray{<:Real}, params, t::Real)
    h = params.h
    N = params.N
    l = params.l
    M = params.M
    rstar = params.rcoord
    boundary_type = params.bc
    @fastmath @inbounds begin
        h2 = 2h
        N1 = N - 1
        Φ = @view U[:, 1]
        Π = @view U[:, 2]
        Ψ = @view U[:, 3]

        dtΦ = @view du[:, 1]
        dtΠ = @view du[:, 2]
        dtΨ = @view du[:, 3]

        for i in 1:N
            #calculate RHS everywhere
            if i == 1
                drΨ = (Ψ[2] - Ψ[1]) / h
                drΠ = (Π[2] - Π[1]) / h
            elseif i == N
                drΨ = (Ψ[end] - Ψ[N1]) / h
                drΠ = (Π[end] - Π[N1]) / h
            else
                drΠ = (Π[i+1] - Π[i-1]) / h2
                drΨ = (Ψ[i+1] - Ψ[i-1]) / h2
            end
            dtΦ[i] = Π[i]
            # Apply boundary conditions
            # absorbing/radiative
            if i == 1
                if boundary_type === :radiative
                    ## Variable change to characteristic vars
                    a = drΨ + drΠ # speed = -1
                    b = drΨ - drΠ # speed = +1
                    ## Set incoming characteristic to 0
                    b = 0
                    ## Change back to primitive variables
                    drΨ = (a + b) / 2
                    drΠ = (a - b) / 2
                    ## ODE
                    dtΠ[i] = drΨ
                    dtΨ[i] = drΠ
                elseif boundary_type === :reflective
                    ## Set ingoing to outgoing characteristic
                    dtΠ[i] = 0 # dxΨx
                    dtΨ[i] = drΠ
                else
                    throw("boundary_type variable can only be :radiative or :reflective, you provided
                           $boundary_type")
                end
            elseif i == N
                if boundary_type === :radiative
                    ## Variable change to characteristic vars
                    a = drΨ + drΠ
                    b = drΨ - drΠ
                    ## Set incoming characteristic to 0
                    a = 0
                    ## Change back to primitive variables
                    drΨ = (a + b) / 2
                    drΠ = (a - b) / 2
                    # ODE
                    dtΠ[i] = drΨ
                    dtΨ[i] = drΠ
                elseif boundary_type === :reflective
                    dtΠ[i] = 0 # dxΨx
                    dtΨ[i] = drΠ
                else
                    throw("boundary_type variable can only be :radiative or :reflective, you provided
                           $boundary_type")
                end
            else
                # ri = InitialData.r(rstar[i], M) # transform rstar to r for potential
                # @assert ri > 2M
                # dtΠ[i] = drΨ + InitialData.Vpot(l, rstar[i], M) * Φ[i]
                dtΠ[i] = drΨ + InitialData.poschl_teller(l, rstar[i]) * Φ[i]
                dtΨ[i] = drΠ
            end
        end
    end
    return nothing
end

end

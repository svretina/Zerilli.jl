module Plot

using Plots
using GridFunctions
using ..InitialData

function plot_ti(u, ti)
    rs = GridFunctions.coords(u[2])

    Φ = u[3][ti][:, 1]
    Π = u[3][ti][:, 2]
    Ψ = u[3][ti][:, 3]

    p = plot(rs, Φ; label="Φ")
    plot!(rs, Π; label="Π")
    plot!(rs, Ψ; label="Ψ")
    return p
end

function trapz(f, h)
    return (h / 2) * (f[1] + f[end] + 2sum(f[2:(end - 1)]))
end

function field_energy(u, ti)
    h = u[4].h
    l = u[4].l
    L = u[4].L
    rs = GridFunctions.coords(u[2])
    Φ = u[3][ti][:, 1]
    Π = u[3][ti][:, 2]
    Ψ = u[3][ti][:, 3]
    # energy = @. 0.5 * (Π * Π + Ψ * Ψ + InitialData.poschl_teller((l), rs) * Φ)
    energy = @. 0.5 * (Π * Π + Ψ * Ψ - cos((π / 2L) * rs) * Φ)

    return trapz(energy, h)
end

function plot_energy(u; show=true)
    t = u[3].t
    n = length(t)
    e = Vector{Float64}(undef, n)

    for i in 1:n
        e[i] = field_energy(u, i)
    end

    show == true && (p = plot(t, e); display(p))
    return t, e
end

function plot_dtE(u)
    t, e = plot_energy(u; show=false)
    dt = t[2] - t[1]
    h = u[4].h
    l = u[4].l
    L = u[4].L
    n = length(t)
    r = GridFunctions.coords(u[2])
    dedt = similar(e)
    pt = similar(e)

    dedt[1] = (e[2] - e[1]) / dt
    dedt[end] = (e[end] - e[end - 1]) / dt
    for i in 1:n
        i > 1 && i < n && (dedt[i] = (e[i + 1] - e[i - 1]) / (2dt))
        Π = u[3][i][:, 2]
        Ψ = u[3][i][:, 3]
        # term = @. 0.5 * Π * InitialData.poschl_teller.(l, r)
        term = @. -0.5 * Π * cos((π / 2L) * r)
        pt[i] = trapz(term, h) + Π[end] * Ψ[end] - Π[1] * Ψ[1]
    end

    p = plot(t, dedt; label="numerical")
    plot!(p, t, pt; ls=:dash, label="semi-analytic")
    display(p)
    return nothing
end

end

module Plot

using Plots
using GridFunctions

function plot_ti(sim_folder, ti)
    rs = GridFunctions.coords(get_grid_from_yaml(sim_folder))

    Φ = h5read(dataset, "Phi")
    Π = h5read(dataset, "Pi")
    Ψ = h5read(dataset, "Psi")

    p = plot(rs, Φ)
    return p
end

function trapz(f, h)
    return (h / 2) * (f[1] + f[end] + 2sum(f[2:(end - 1)]))
end

function field_energy()
    energy = @. 0.5 * (ppi^2 + psi^2 + InitialData.poschl_teller(l, rs) * phi)
    return trapz(energy, h)
end

function plot_energy(sim_folder)
    sims = readdir(string(sim_folder, "/sims"); join=true)
    r = get_grid_from_yaml(sim_folder)
    h = spacing(r)
    n = length(sims)
    e = Vector{Float64}(undef, n)
    t = GridFunctions.coords(get_time_grid_from_yaml(sim_folder))
    for (i, dataset) in enumerate(sims)
        Π = h5read(dataset, "Pi")
        Ψ = h5read(dataset, "Psi")
        e[i] = field_energy(h, Π, Ψ)
    end
    p = plot(t, e)
    return p
end

end

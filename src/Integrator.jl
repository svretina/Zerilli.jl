module Integrator

using HDF5
using LoopVectorization

const proj_path = pkgdir(Integrator)
const output_path = string(proj_path, "/output")

function get_iter_str(i, N)
    i = string(i)
    while length(i) < length(string(N))
        i = string("0", i)
    end
    if length(i) <= 2
        return string("0", i)
    else
        return i
    end
end

function base_path(folder)
    return string(output_path, folder)
end

function sims_path(folder)
    return string(output_path, folder, "/sims/")
end

@inline function RK4(rhs!::F, dt::Real, reg1, reg2, reg3, reg4, params, t::Real,
    indices::CartesianIndices) where {F}
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + reg3[i] / 2
        reg2[i] = reg3[i]
    end
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + (reg3[i] - reg2[i]) / 2
    end
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] - reg3[i] / 2
        reg1[i] = reg1[i] + reg3[i]
        reg2[i] = reg2[i] / 6 - reg3[i]
    end
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] + reg3[i] + reg3[i]
        reg1[i] = reg1[i] + reg2[i] + reg3[i] / 6
    end
    return nothing
end

function solve(rhs, statevector, params, t, grid, M, save_every, folder, overwrite)
    tnpoints = t.ncells + 1
    istr = get_iter_str(0, tnpoints)
    base_dir = base_path(folder)
    sims_dir = sims_path(folder)

    base_str = string(sims_dir, "output_M=", M, "_L=", Int64(grid.domain[2]), "_nc=$(grid.ncells)_ti=")
    dataset = string(base_str, istr, ".h5")

    if !isdir(sims_dir)
        mkdir(base_dir)
        mkdir(sims_dir)
    end

    if overwrite == true && isfile(dataset)
        rm(dataset)
    end
    print("Writing initial data at: ", dataset)
    h5write(dataset, "Phi", statevector[:, 1])
    h5write(dataset, "Pi", statevector[:, 2])
    h5write(dataset, "Psi", statevector[:, 3])
    println("✅")
    println("=========Starting time integration=========")

    print("Allocating registers for time integrator...")
    reg2 = similar(statevector)
    reg3 = similar(statevector)
    reg4 = similar(statevector)
    println("✅")

    println("saving every=", save_every)
    nt = t.ncells + 1
    dt = params.dt
    indices = CartesianIndices(statevector)
    write_metadata(base_dir, params, grid, M, save_every, overwrite, dataset)

    for (i, ti) in enumerate(params.ti)
        if ti == 0.0
            continue
        end
        println("Iteration = ", i, "/", nt)
        @time RK4(rhs, dt, statevector, reg2, reg3, reg4, params, ti, indices)
        if i % save_every == 0
            istr = get_iter_str(i, tnpoints)
            dataset = string(base_str, istr, ".h5")
            @time h5open(dataset, "w") do file
                write(file, "Phi", statevector[:, 1])
                write(file, "Pi", statevector[:, 2])
                write(file, "Psi", statevector[:, 3])
            end
        end
    end
    return nothing
end

function write_metadata(md_path, params, grid, M, save_every, overwrite, dataset)
    metadata_filename = string(md_path, "/metadata.md")
    data = "## Grid
domain:       [$(grid.domain[1]), $(grid.domain[2])]
ncells:       $(grid.ncells)
spacing:      $(params.h)

## Time
CFL:          $(params.dt/params.h)
domain:       [0.0, $(params.ti[end])]
ncells:       $(length(params.ti)-1)
dt:           $(params.dt)

## Black Hole
mass:         $M

## Output
output path:  $dataset
save every:   $save_every
overwrite:    $overwrite
\n
"
    println(data)
    open(metadata_filename, "w") do file
        write(file, data)
    end
end

end #end of module

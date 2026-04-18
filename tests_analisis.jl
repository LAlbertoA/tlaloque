using VTKLegacy
using CairoMakie
#cd("tests/Sod")

function analyze(name)
    if name == "Sod"
        _sod_analisis()
    elseif name == "Sedov-Taylor"
    elseif name == "Radiative Sedov-Taylor"
    elseif name == "Gravitational Potential Accuracy Test (Multigrid)"
    elseif name == "Evrard's Collapse"
    elseif name == "Truelove Collapse Test"
    end
end

function _sod_analisis()
    files = readdir("./DATA")
    c = 0
    for i in files
        if i != "output.txt"
            c += 1
            m = LoadVTK("DATA/"*i);
            halfsize_index = floor.(Int,[m.nx, m.ny, m.nz]./2)

            f, ax, h = heatmap(m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], axis = (title = "Slice at $(halfsize_index[3]) in z",))
            cb = Colorbar(f[1,2], h)
            
            save("Sod"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]

            f, ax = lines(r, rho, axis = (title = "Diagonal profile",))

            save("Sod_profile"*string(c,base=10,pad=3)*".png", f)
        end
    end
end
#sod_analisis()
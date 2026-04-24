using VTKLegacy
using CairoMakie
#cd("tests/Sod")

function analyze(name)
    if name == "Sod"
        _sod_analisis()
    elseif name == "Sedov-Taylor"
        _sedov_analisis()
    elseif name == "Radiative Sedov-Taylor"
        _radiative_sedov_analisis()
    elseif name == "Gravitational Potential Accuracy Test (Multigrid)"
        _gravitational_potential_analisis()
    elseif name == "Point Potential Test"
        _point_potential_analisis()
    elseif name == "Evrard's Collapse"
        _evrard_collapse_analisis()
    elseif name == "Truelove Collapse Test"
        _truelove_collapse_analisis()
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

            f, ax, h = heatmap(m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], axis = (title = "Slice at $(halfsize_index[3]) in z",
                xlabel = L"x", ylabel = L"y"))
            cb = Colorbar(f[1,2], h, label = L"\rho")
            
            save("Sod"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]

            f, ax = lines(r, rho, axis = (title = "Diagonal profile", xlabel = L"r", ylabel = L"\rho"))

            save("Sod_profile"*string(c,base=10,pad=3)*".png", f)
        end
    end
end

function _sedov_analisis()
    files = readdir("./DATA")
    c = 0
    for i in files
        if i != "output.txt"
            c += 1
            m = LoadVTK("DATA/"*i);
            halfsize_index = floor.(Int,[m.nx, m.ny, m.nz]./2)

            f, ax, h = heatmap(m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], axis = (title = "Slice at $(halfsize_index[3]) in z",
                xlabel = L"x", ylabel = L"y"))
            cb = Colorbar(f[1,2], h, label = L"\rho")
            
            save("Sedov"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]

            f, ax = lines(r, rho, axis = (title = "Diagonal profile", xlabel = L"r", ylabel = L"\rho"))

            save("Sedov_profile"*string(c,base=10,pad=3)*".png", f)
        end
    end
end

function _radiative_sedov_analisis()
    files = readdir("./DATA")
    c = 0
    for i in files
        if i != "output.txt"
            c += 1
            m = LoadVTK("DATA/"*i);
            halfsize_index = floor.(Int,[m.nx, m.ny, m.nz]./2)

            f, ax, h = heatmap(m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], axis = (title = "Slice at $(halfsize_index[3]) in z",
                xlabel = L"x", ylabel = L"y"))
            cb = Colorbar(f[1,2], h, label = L"\rho")
            
            save("Radiative_Sedov"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]

            f, ax = lines(r, rho, axis = (title = "Diagonal profile", xlabel = L"r", ylabel = L"\rho"))

            save("Radiative_Sedov_profile"*string(c,base=10,pad=3)*".png", f)
        end
    end
end

function _evrard_collapse_analisis()
    files = readdir("./DATA")
    c = 0
    for i in files
        if i != "output.txt"
            c += 1
            m = LoadVTK("DATA/"*i);
            halfsize_index = floor.(Int,[m.nx, m.ny, m.nz]./2)

            f, ax, h = heatmap(m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], axis = (title = "Slice at $(halfsize_index[3]) in z",
                xlabel = L"x", ylabel = L"y"))
            cb = Colorbar(f[1,2], h, label = L"\rho")
            
            save("Evrard_Collapse"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]

            f, ax = lines(r, rho, axis = (title = "Diagonal profile", xlabel = L"r", ylabel = L"\rho"))

            save("Evrard_Collapse_profile"*string(c,base=10,pad=3)*".png", f)
        end
    end
end

function _truelove_collapse_analisis()
    files = readdir("./DATA")
    c = 0
    for i in files
        if i != "output.txt"
            c += 1
            m = LoadVTK("DATA/"*i);
            halfsize_index = floor.(Int,[m.nx, m.ny, m.nz]./2)

            f, ax, h = heatmap(m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], axis = (title = "Slice at $(halfsize_index[3]) in z",
                xlabel = L"x", ylabel = L"y"))
            cb = Colorbar(f[1,2], h, label = L"\rho")
            
            save("Truelove_Collapse"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]

            f, ax = lines(r, rho, axis = (title = "Diagonal profile", xlabel = L"r", ylabel = L"\rho"))

            save("Truelove_Collapse_profile"*string(c,base=10,pad=3)*".png", f)
        end
    end
end
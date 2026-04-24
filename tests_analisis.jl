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
            vel = magnitude(m.pointData[3:5,:,:,:])
            Uin = m.pointData[2,:,:,:]./(((5/3)-1).*m.pointData[1,:,:,:])

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Slices at $(halfsize_index[3]) in z", xlabel = L"x", ylabel = L"y")
            h1 = heatmap!(ax1, m.x, m.y, m.pointData[1,:,:,halfsize_index[3]])
            Colorbar(f[1,2], h1, label = L"\rho")

            ax2 = Axis(f[1,3], xlabel = L"x", ylabel = L"y")
            h2 = heatmap!(ax2, m.x, m.y, m.pointData[2,:,:,halfsize_index[3]])
            Colorbar(f[1,4], h2, label = L"P")

            ax3 = Axis(f[2,1], xlabel = L"x", ylabel = L"y")
            h3 = heatmap!(ax3, m.x, m.y, vel[:,:,halfsize_index[3]])
            Colorbar(f[2,2], h3, label = L"|v|")

            ax4 = Axis(f[2,3], xlabel = L"x", ylabel = L"y")
            h4 = heatmap!(ax4, m.x, m.y, Uin[:,:,halfsize_index[3]])
            Colorbar(f[2,4], h4, label = L"U")

            save("Sod"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]
            pres = [m.pointData[2,i,i,i] for i in 1:m.nx]
            rvel = [vel[i,i,i] for i in 1:m.nx]
            rUin = [Uin[i,i,i] for i in 1:m.nx]

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Diagonal profiles", xlabel = L"r", ylabel = L"\rho")
            lines!(ax1, r, rho)

            ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"P")
            lines!(ax2, r, pres)

            ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v|")
            lines!(ax3, r, rvel)

            ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"U")
            lines!(ax4, r, rUin)

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
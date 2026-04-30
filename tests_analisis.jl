using VTKLegacy
using CairoMakie

include("utils_analisis.jl")
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

            save("Sod_slice"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]
            pres = [m.pointData[2,i,i,i] for i in 1:m.nx]
            rvel = [vel[i,i,i] for i in 1:m.nx]
            rUin = [Uin[i,i,i] for i in 1:m.nx]

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Diagonal profiles", xlabel = L"r", ylabel = L"\rho", limits = (-0.25,0.25, nothing, nothing))
            lines!(ax1, r, rho, label = "Numerical")

            ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"P", limits = (-0.25,0.25, nothing, nothing))
            lines!(ax2, r, pres, label = "Numerical")

            ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v|", limits = (-0.25,0.25, nothing, nothing))
            lines!(ax3, r, rvel, label = "Numerical")

            ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"U", limits = (-0.25,0.25, nothing, nothing))
            lines!(ax4, r, rUin, label = "Numerical")

            save("Sod_profile"*string(c,base=10,pad=3)*".png", f)

            if c == length(files)-1
                xex = collect(0.00025:0.0005:0.5).-0.25
                SODXS = Array{Float64,2}(undef,5,1000)
                open("../../solutions/sod_e1rpex.out") do f
                    line = 1
                    while ! eof(f)
                        SODXS[:,line] = parse.(Float64,split(readline(f)))
                        line += 1
                    end
                end
                lines!(ax1, xex, SODXS[2,:], color=:red, label = "Exact")
                lines!(ax2, xex, SODXS[4,:], color=:red, label = "Exact")
                lines!(ax3, xex, SODXS[3,:], color=:red, label = "Exact")
                lines!(ax4, xex, SODXS[5,:], color=:red, label = "Exact")
                axislegend(ax1, position = :rt)

                save("Sod_profile"*string(c,base=10,pad=3)*".png", f)

                errd, maped,rmsped,L1d,L2d,Linfd = _all_error_metrics(xex, SODXS[2,:], r, rho)
                errp, mapep,rmspep,L1p,L2p,Linfp = _all_error_metrics(xex, SODXS[4,:], r, pres)
                errv, mapev,rmspev,L1v,L2v,Linfv = _all_error_metrics(xex, SODXS[3,:], r, rvel)
                errU, mapeU,rmspeU,L1U,L2U,LinfU = _all_error_metrics(xex, SODXS[5,:], r, rUin)

                idx = findfirst(x->x>=xex[1],r)
                fdx = findlast(x->x<=xex[end],r)

                f = Figure(size = (1200,1200), fontsize = 20)
                ax1 = Axis(f[1,1], title = "Diagonal profiles error", xlabel = L"r", ylabel = L"\rho_{exact} - \rho_{numerical}|/|\rho_{exact}|")
                lines!(ax1, r[idx:fdx], errd, color=:blue, label = "Density")

                ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"|P_{exact} - P_{numerical}|/|P_{exact}|")
                lines!(ax2, r[idx:fdx], errp, color=:blue, label = "Pressure")

                ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v_{exact} - v_{numerical}|/|v_{exact}|")
                lines!(ax3, r[idx:fdx], errv, color=:blue, label = "Velocity")

                ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"|U_{exact} - U_{numerical}|/|U_{exact}|")
                lines!(ax4, r[idx:fdx], errU, color=:blue, label = "U")

                save("Sod_profile_error"*string(c,base=10,pad=3)*".png", f)

                open("Sod_error_metrics.txt", "w") do f
                    println(f, "Density:")
                    println(f, "MAPE: ", maped)
                    println(f, "RMSPE: ", rmsped)
                    println(f, "L1: ", L1d)
                    println(f, "L2: ", L2d)
                    println(f, "Linf: ", Linfd)
                    println(f, "\nPressure:")
                    println(f, "MAPE: ", mapep)
                    println(f, "RMSPE: ", rmspep)
                    println(f, "L1: ", L1p)
                    println(f, "L2: ", L2p)
                    println(f, "Linf: ", Linfp)
                    println(f, "\nVelocity:")
                    println(f, "MAPE: ", mapev)
                    println(f, "RMSPE: ", rmspev)
                    println(f, "L1: ", L1v)
                    println(f, "L2: ", L2v)
                    println(f, "Linf: ", Linfv)
                    println(f, "\nU:")
                    println(f, "MAPE: ", mapeU)
                    println(f, "RMSPE: ", rmspeU)
                    println(f, "L1: ", L1U)
                    println(f, "L2: ", L2U)
                    println(f, "Linf: ", LinfU)
                end
            end
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
            
            save("Sedov_slice"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]
            pres = [m.pointData[2,i,i,i] for i in 1:m.nx]
            rvel = [vel[i,i,i] for i in 1:m.nx]
            rUin = [Uin[i,i,i] for i in 1:m.nx]

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Diagonal profiles", xlabel = L"r", ylabel = L"\rho", limits = (-5,5, nothing, nothing))
            lines!(ax1, r, rho, label = "Numerical")

            ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"P", limits = (-5,5, nothing, nothing))
            lines!(ax2, r, pres, label = "Numerical")

            ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v|", limits = (-5,5, nothing, nothing))
            lines!(ax3, r, rvel, label = "Numerical")

            ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"U", limits = (-5,5, nothing, nothing))
            lines!(ax4, r, rUin, label = "Numerical")

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
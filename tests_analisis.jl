using VTKLegacy
using CairoMakie

include("constants.jl")
include("utils_analisis.jl")

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
    elseif name == "Evrard Collapse"
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
            theta = collect(0.001:0.001:2*pi)
            vel = magnitude(m.pointData[3:5,:,:,:])
            Uin = m.pointData[2,:,:,:]./(((5/3)-1).*m.pointData[1,:,:,:])

            f = Figure(size = (1400,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Slices at $(halfsize_index[3]) in z", xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h1 = heatmap!(ax1, m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[1,2], h1, label = L"\rho")

            ax2 = Axis(f[1,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h2 = heatmap!(ax2, m.x, m.y, m.pointData[2,:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[1,4], h2, label = L"P")

            ax3 = Axis(f[2,1], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h3 = heatmap!(ax3, m.x, m.y, vel[:,:,halfsize_index[3]])
            Colorbar(f[2,2], h3, label = L"|v|")

            ax4 = Axis(f[2,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h4 = heatmap!(ax4, m.x, m.y, Uin[:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[2,4], h4, label = L"U")

            if c != 1
                time_sedov = (c-1)*0.1/(length(files)-2)
                radius_shock = _radius_sedov(time_sedov, 1e5, 1)
                xrshock, yrshock = radius_shock.*cos.(theta), radius_shock.*sin.(theta)
                lines!(ax1, xrshock, yrshock, color = :red, label = "Shock expected position")
                lines!(ax2, xrshock, yrshock, color = :red, label = "Shock expected position")
                lines!(ax3, xrshock, yrshock, color = :red, label = "Shock expected position")
                lines!(ax4, xrshock, yrshock, color = :red, label = "Shock expected position")
            end

            save("Sedov_slice"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]
            pres = [m.pointData[2,i,i,i] for i in 1:m.nx]
            rvel = [vel[i,i,i] for i in 1:m.nx]
            rUin = [Uin[i,i,i] for i in 1:m.nx]

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Diagonal profiles", xlabel = L"r", ylabel = L"\rho", yscale = log10)
            lines!(ax1, r, rho, label = "Numerical")

            ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"P", yscale = log10)
            lines!(ax2, r, pres, label = "Numerical")

            ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v|")
            lines!(ax3, r, rvel, label = "Numerical")

            ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"U", yscale = log10)
            lines!(ax4, r, rUin, label = "Numerical")

            if c != 1
                time_sedov = (c-1)*0.1/(length(files)-2)
                radius_shock = _radius_sedov(time_sedov, 1e5, 1)
                lines!(ax1, [radius_shock, radius_shock], [minimum(rho), maximum(rho)], color = :black, 
                    label = "Shock expected position")
                lines!(ax2, [radius_shock, radius_shock], [minimum(pres), maximum(pres)], color = :black, 
                    label = "Shock expected position")
                lines!(ax3, [radius_shock, radius_shock], [minimum(rvel), maximum(rvel)], color = :black, 
                    label = "Shock expected position")
                lines!(ax4, [radius_shock, radius_shock], [minimum(rUin), maximum(rUin)], color = :black, 
                    label = "Shock expected position")
                axislegend(ax1, position = :lt)
            end

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
            theta = collect(0.001:0.001:2*pi)
            vel = magnitude(m.pointData[3:5,:,:,:])
            Uin = m.pointData[2,:,:,:]./(((5/3)-1).*m.pointData[1,:,:,:])

            f = Figure(size = (1400,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Slices at $(halfsize_index[3]) in z", xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h1 = heatmap!(ax1, m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[1,2], h1, label = L"\rho")

            ax2 = Axis(f[1,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h2 = heatmap!(ax2, m.x, m.y, m.pointData[2,:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[1,4], h2, label = L"P")

            ax3 = Axis(f[2,1], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h3 = heatmap!(ax3, m.x, m.y, vel[:,:,halfsize_index[3]])
            Colorbar(f[2,2], h3, label = L"|v|")

            ax4 = Axis(f[2,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h4 = heatmap!(ax4, m.x, m.y, Uin[:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[2,4], h4, label = L"U")

            if c != 1
                time_sedov = (c-1)*50*YEAR*1e3/(length(files)-2)
                radius_shock = _radius_sedov(time_sedov, 0.5e51, 1*AMU)
                xrshock, yrshock = radius_shock.*cos.(theta), radius_shock.*sin.(theta)
                lines!(ax1, xrshock, yrshock, color = :red, label = "Shock expected position")
                lines!(ax2, xrshock, yrshock, color = :red, label = "Shock expected position")
                lines!(ax3, xrshock, yrshock, color = :red, label = "Shock expected position")
                lines!(ax4, xrshock, yrshock, color = :red, label = "Shock expected position")
            end
            
            save("Radiative_Sedov_slice"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]
            pres = [m.pointData[2,i,i,i] for i in 1:m.nx]
            rvel = [vel[i,i,i] for i in 1:m.nx]
            rUin = [Uin[i,i,i] for i in 1:m.nx]

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Diagonal profiles", xlabel = L"r", ylabel = L"\rho", limits = (-30*PC,30*PC, nothing, nothing), yscale = log10)
            lines!(ax1, r, rho, label = "Numerical")

            ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"P", limits = (-30*PC,30*PC, nothing, nothing), yscale = log10)
            lines!(ax2, r, pres, label = "Numerical")

            ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v|", limits = (-30*PC,30*PC, nothing, nothing))
            lines!(ax3, r, rvel, label = "Numerical")

            ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"U", limits = (-30*PC,30*PC, nothing, nothing), yscale = log10)
            lines!(ax4, r, rUin, label = "Numerical")

            if c != 1
                time_sedov = (c-1)*50*YEAR*1e3/(length(files)-2)
                radius_shock = _radius_sedov(time_sedov, 0.5e51, 1*AMU)
                lines!(ax1, [radius_shock,radius_shock], [minimum(rho), maximum(rho)], color = :black, 
                    label = "Shock expected position")
                lines!(ax2, [radius_shock,radius_shock], [minimum(pres), maximum(pres)], color = :black, 
                    label = "Shock expected position")
                lines!(ax3, [radius_shock,radius_shock], [minimum(rvel), maximum(rvel)], color = :black, 
                    label = "Shock expected position")
                lines!(ax4, [radius_shock,radius_shock], [minimum(rUin), maximum(rUin)], color = :black, 
                    label = "Shock expected position")
                axislegend(ax1, position = :lt)
            end

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
            vel = magnitude(m.pointData[4:6,:,:,:])

            f = Figure(size = (1400,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Slices at $(halfsize_index[3]) in z", xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h1 = heatmap!(ax1, m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[1,2], h1, label = L"\rho")

            ax2 = Axis(f[1,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h2 = heatmap!(ax2, m.x, m.y, m.pointData[2,:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[1,4], h2, label = L"P")

            ax3 = Axis(f[2,1], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h3 = heatmap!(ax3, m.x, m.y, vel[:,:,halfsize_index[3]])
            Colorbar(f[2,2], h3, label = L"|v|")

            ax4 = Axis(f[2,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h4 = heatmap!(ax4, m.x, m.y, m.pointData[3,:,:,halfsize_index[3]])
            Colorbar(f[2,4], h4, label = L"\Phi")
            
            save("Evrard_Collapse_slice"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]
            pres = [m.pointData[2,i,i,i] for i in 1:m.nx]
            rvel = [vel[i,i,i] for i in 1:m.nx]
            rphi = [m.pointData[3,i,i,i] for i in 1:m.nx]

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Diagonal profiles", xlabel = L"r", ylabel = L"\rho", limits = (-2,2, nothing, nothing), yscale = log10)
            lines!(ax1, r, rho, label = "Numerical")

            ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"P", limits = (-2,2, nothing, nothing), yscale = log10)
            lines!(ax2, r, pres, label = "Numerical")

            ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v|", limits = (-2,2, nothing, nothing))
            lines!(ax3, r, rvel, label = "Numerical")

            ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"\Phi", limits = (-2,2, nothing, nothing))
            lines!(ax4, r, rphi, label = "Numerical")

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
            vel = magnitude(m.pointData[4:6,:,:,:])

            if c == 1
                cscale = identity
            else
                cscale = log10
            end

            f = Figure(size = (1400,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Slices at $(halfsize_index[3]) in z", xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h1 = heatmap!(ax1, m.x, m.y, m.pointData[1,:,:,halfsize_index[3]], colorscale = log10)
            Colorbar(f[1,2], h1, label = L"\rho")

            ax2 = Axis(f[1,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h2 = heatmap!(ax2, m.x, m.y, m.pointData[2,:,:,halfsize_index[3]], colorscale = cscale)
            Colorbar(f[1,4], h2, label = L"P")

            ax3 = Axis(f[2,1], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h3 = heatmap!(ax3, m.x, m.y, vel[:,:,halfsize_index[3]])
            Colorbar(f[2,2], h3, label = L"|v|")

            ax4 = Axis(f[2,3], xlabel = L"x", ylabel = L"y", yticklabelrotation = pi/2)
            h4 = heatmap!(ax4, m.x, m.y, m.pointData[3,:,:,halfsize_index[3]])
            Colorbar(f[2,4], h4, label = L"\Phi")
            
            save("Truelove_Collapse_slice"*string(c,base=10,pad=3)*".png", f)

            r = [sign(i)*sqrt(i^2 + i^2 + i^2) for i in m.x]
            rho = [m.pointData[1,i,i,i] for i in 1:m.nx]
            pres = [m.pointData[2,i,i,i] for i in 1:m.nx]
            rvel = [vel[i,i,i] for i in 1:m.nx]
            rphi = [m.pointData[3,i,i,i] for i in 1:m.nx]

            f = Figure(size = (1200,1200), fontsize = 20)
            ax1 = Axis(f[1,1], title = "Diagonal profiles", xlabel = L"r", ylabel = L"\rho", limits = (-535.0*AU,535.0*AU, nothing, nothing), yscale = log10)
            lines!(ax1, r, rho, label = "Numerical")

            ax2 = Axis(f[1,2], xlabel = L"r", ylabel = L"P", limits = (-535.0*AU,535.0*AU, nothing, nothing), yscale = cscale)
            lines!(ax2, r, pres, label = "Numerical")

            ax3 = Axis(f[2,1], xlabel = L"r", ylabel = L"|v|", limits = (-535.0*AU,535.0*AU, nothing, nothing))
            lines!(ax3, r, rvel, label = "Numerical")

            ax4 = Axis(f[2,2], xlabel = L"r", ylabel = L"\Phi", limits = (-535.0*AU,535.0*AU, nothing, nothing))
            lines!(ax4, r, rphi, label = "Numerical")

            save("Truelove_Collapse_profile"*string(c,base=10,pad=3)*".png", f)
        end
    end
end
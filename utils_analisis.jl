function _radius_sedov(t, E0, rho0)
    return 1.15*(E0*t^2/rho0)^(1/5)
end
function _interpol(p1,p2)
    x1, y1 = p1
    x2, y2 = p2
    
    ls = [x2-x1,y2-y1]
    m = ls[2]/ls[1]
    b = (x2*y1-y2*x1)/(x2-x1)

    return m, b
end

function _neighbor_search(p, arr, start_idx)
    for i in start_idx:length(arr)
        if arr[i] > p
            return i-1,i
        end
    end
    return nothing, nothing
end
function _all_error_metrics(xex, yex, xnum, ynum)
    if xex[1] > xnum[1] 
        idx = findfirst(x->x>=xex[1],xnum)
    else
        idx = 1
    end

    if xex[end] < xnum[end]
        fdx = findlast(x->x<=xex[end],xnum)
    else
        fdx = length(xnum)
    end

    j = 1
    k = 2
    err = []
    mapes = 0
    rmspe = 0
    L1n = 0
    L1d = 0
    L2n = 0
    L2d = 0
    scale = maximum(abs.(yex))
    for i in idx:fdx
        j, k = _neighbor_search(xnum[i], xex, k)
        if j == nothing || k == nothing
            return nothing, nothing, nothing, nothing, nothing, nothing
        end
        m, b = _interpol((xex[j],yex[j]),(xex[k],yex[k]))
        yex_interp = m*xnum[i] + b
        y_err_abs = abs(yex_interp - ynum[i])
        if abs(yex_interp) > 1e-12
            y_err_rel = y_err_abs / abs(yex_interp)
        else
            y_err_rel = y_err_abs / scale
        end
        push!(err, y_err_rel)
        mapes += y_err_rel
        rmspe += y_err_rel^2
        L1n += y_err_abs
        L1d += abs(yex_interp)
        L2n += y_err_abs^2
        L2d += yex_interp^2
    end
    mapes = mapes/(length(err))
    rmspe = sqrt(rmspe/(length(err)))
    L1 = L1n/L1d
    L2 = sqrt(L2n/L2d)
    return err, mapes, rmspe, L1, L2, maximum(err)
end

function _percent_error(yex, ynum)
    return 100*abs.(yex .- ynum) ./ abs.(yex)
end

function _norm_L1(yex, ynum)
    numerator = sum(abs.(yex .- ynum))
    denominator = sum(abs.(yex))
    return numerator/denominator
end

function _norm_L2(yex, ynum)
    numerator = sum((yex .- ynum).^2)
    denominator = sum(yex.^2)
    return sqrt(numerator / denominator)
end

function _norm_Linf(yex, ynum)
    numerator = maximum(abs.(yex .- ynum))
    denominator = maximum(abs.(yex))
    return numerator/denominator
end
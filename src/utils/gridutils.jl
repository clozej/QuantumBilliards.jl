using LazyGrids
#used for plotting etc.
function lazy_grid(xlim, ylim, grd::Tuple)
    x = collect(LinRange(xlim...,grd[1]))
    y = collect(LinRange(ylim...,grd[2]))
    (xg, yg) = ndgrid(x, y)
    return x, y, xg, yg
end

#x, y, xg, yg = lazy_grid((0,1), (0,1), (100,100))

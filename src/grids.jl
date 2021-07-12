"""
    struct Grid1D

Returns a one-dimensional staggered `grid`.

$(TYPEDFIELDS)
"""
struct Grid1D
    "Number of points in x-direction"
    nx::Int
    "Grid spacing in x-direction"
    dx::Float64
    "Domain extent in x-direction"
    Lx::Float64
    "U-grid in x-direction"
    xu::Vector
    "T-grid in x-direction"
    xt::Vector
end

function Grid1D(nx, x_start, x_end)
    Lx = x_end - x_start
    dx = Lx/nx
    xu = range(x_start, stop = x_end-dx, length = nx)
    xt = xu .+ dx/2

    return Grid1D(nx, dx, Lx, xu, xt)
end

"""
    struct Grid2D

Returns a two-dimensional `grid`.

$(TYPEDFIELDS)
"""
struct Grid2D
    "Number of points in x-direction"
    nx::Int
    "Number of points in y-direction"
    ny::Int
    "Grid spacing in x-direction"
    dx::Float64
    "Grid spacing in y-direction"
    dy::Float64
    "Domain extent in x-direction"
    Lx::Float64
    "Domain extent in y-direction"
    Ly::Float64
    "U-grid in x-direction"
    xu::Vector
    "T-grid in x-direction"
    xt::Vector
    "U-grid in y-direction"
    yu::Vector
    "T-grid in y-direction"
    yt::Vector
end

function Grid2D(nx, ny, x_start, x_end, y_start, y_end)
    Lx = x_end - x_start
    Ly = y_end - y_start
    dx = Lx/nx
    dy = Ly/ny
    xu = range(x_start, stop = x_end - dx, length = nx)
    xt = xu .+ dx/2
    yu = range(y_start, stop = y_end - dy, length = ny)
    yt = yu .+ dy/2
    
    return Grid2D(nx, ny, dx, dy, Lx, Ly, xu, xt, yu, yt)
end

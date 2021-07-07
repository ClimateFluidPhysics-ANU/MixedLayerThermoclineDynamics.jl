"""
    struct Grid1D

A one-dimensional `grid`.

$(TYPEDFIELDS)
"""
struct Grid1D
    "number of points in x"
    nx::Int
    "grid spacing in x"
    dx::Float64
    "domain extent in x"
    Lx::Float64
    "to be filled..."
    xu::Vector
    "blah blah..."
    xt::Vector
end

function Grid1D(nx, x_beg, x_end)
    dx = (x_end - x_beg)/(nx + 0.5)
    xu = LinRange(x_beg, x_end - dx/2, nx)
    xt = LinRange(x_beg + dx/2, x_end, nx)

    return Grid1D(nx, dx, Lx, xu, xt)
end

struct Grid2D
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    Lx::Float64
    Ly::Float64
    xu::Vector
    xt::Vector
    yu::Vector
    yt::Vector
end

function Grid2D(nx, ny, x_beg, x_end, y_beg, y_end)
    dx = (x_end - x_beg)/(nx + 0.5)
    dy = (y_end - y_beg)/(ny + 0.5)
    xu = LinRange(x_beg, x_end - dx/2, nx)
    xt = LinRange(x_beg + dx/2, x_end, nx)
    yu = LinRange(y_beg, y_end - dy/2, ny)
    yt = LinRange(y_beg + dy/2, y_end, ny)
    
    return Grid2D(nx, ny, dx, dy, Lx, Ly, xu, xt, yu, yt)
end
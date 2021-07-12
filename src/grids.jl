"""
    struct Grid1D <: AbstractGrid

Returns a one-dimensional staggered `grid`.

$(TYPEDFIELDS)
"""
struct Grid1D <: AbstractGrid
    "Number of points in x-direction"
    nx::Int
    "Grid spacing in x-direction"
    dx::Float64
    "Domain extent in x-direction"
    Lx::Float64
    "Cell faces in x-direction"
    xF::Vector
    "Cell centres in x-direction"
    xC::Vector
end

function Grid1D(nx, x_start, x_end)
    Lx = x_end - x_start
    dx = Lx/nx
    xF = range(x_start, stop = x_end-dx, length = nx)
    xC = xF .+ dx/2

    return Grid1D(nx, dx, Lx, xF, xC)
end

"""
    struct Grid2D <: AbstractGrid

Returns a two-dimensional staggered `grid`.

$(TYPEDFIELDS)
"""
struct Grid2D <: AbstractGrid
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
    "Cell faces in x-direction"
    xF::Vector
    "Cell centres in x-direction"
    xC::Vector
    "Cell faces in y-direction"
    yF::Vector
    "Cell centres in y-direction"
    yC::Vector
end

function Grid2D(nx, ny, x_start, x_end, y_start, y_end)
    Lx = x_end - x_start
    Ly = y_end - y_start
    dx = Lx/nx
    dy = Ly/ny
    xF = range(x_start, stop = x_end - dx, length = nx)
    xC = xF .+ dx/2
    yF = range(y_start, stop = y_end - dy, length = ny)
    yC = yF .+ dy/2
    
    return Grid2D(nx, ny, dx, dy, Lx, Ly, xF, xC, yF, yC)
end

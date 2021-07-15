"""
    struct Grid1D{Tx<:AbstractTopology} <: AbstractGrid

Returns a one-dimensional staggered `grid` with topology `Tx`.

$(TYPEDFIELDS)
"""
struct Grid1D{Tx<:AbstractTopology} <: AbstractGrid
    "Number of points in x-direction"
    nx::Int
    "Number of halo points in x-direction"
    hx::Int
    "Grid spacing in x-direction"
    dx::Float64
    "Domain extent in x-direction"
    Lx::Float64
    "Cell faces in x-direction"
    xF::Vector
    "Cell centres in x-direction"
    xC::Vector
end

function Grid1D(Tx, nx, x_start, x_end; hx=0)
    Lx = x_end - x_start
    dx = Lx/nx
    xF = construct_faces(Tx, nx, hx, dx, Lx, x_start)
    xC = xF .+ dx/2

    return Grid1D{typeof(Tx)}(nx, hx, dx, Lx, xF, xC)
end

"""
    Grid2D{Tx<:AbstractTopology, Ty<:AbstractTopology} <: AbstractGrid

Returns a two-dimensional staggered `grid` with topologies `{Tx, Ty}`.

$(TYPEDFIELDS)
"""
struct Grid2D{Tx<:AbstractTopology, Ty<:AbstractTopology} <: AbstractGrid
    "Number of points in x-direction"
    nx::Int
    "Number of points in y-direction"
    ny::Int
    "Number of halo points in x-direction"
    hx::Int
    "Number of halo points in y-direction"
    hy::Int
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

function Grid2D(Tx, Ty, nx, ny, x_start, x_end, y_start, y_end; hx=0, hy=0)
    Lx = x_end - x_start
    Ly = y_end - y_start
    dx = Lx/nx
    dy = Ly/ny
    xF = construct_faces(Tx, nx, hx, dx, Lx, x_start)
    xC = xF .+ dx/2
    yF = construct_faces(Ty, ny, hy, dy, Ly, y_start)
    yC = yF .+ dy/2
    
    return Grid2D{typeof(Tx), typeof(Ty)}(nx, ny, hx, hy, dx, dy, Lx, Ly, xF, xC, yF, yC)
end

##
# helper functions

"""
    construct_faces(T::AbstractTopology, n, h, d, L)

Returns the face locations for a dimension with topology `T`, number of grid points `n` and `h`
halo points, grid spacing `d`, and extent `L`.
"""
function construct_faces(T::AbstractTopology, n, h, d, L, startface)
    
    end_face = isa(T, Periodic) ? L - d : L
    total_faces = isa(T, Periodic) ? n : n + 1
    
    return range(startface - h*d, stop = end_face + h*d, length = total_faces)
end

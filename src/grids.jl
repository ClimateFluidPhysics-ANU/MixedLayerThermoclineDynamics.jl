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
    xF::AbstractArray
    "Cell centres in x-direction"
    xC::AbstractArray
end

"""
    Grid1D(Tx, nx, x_start, x_end; hx=0)
Construct a one-dimensional staggered `grid` on domain `x ∈ [x_start, x_end]` with topology
`Tx`, with `nx` interior grid points, and `hx` halo points.
"""
function Grid1D(Tx, nx, x_start, x_end; hx=1)

    if hx == 0; throw(error("Number of halo points cannot be zero")); end

    Lx = x_end - x_start
    
    dx = Lx/nx
    
    xF = construct_faces(Tx, nx, hx, dx, Lx, x_start)
    xC = construct_centres(Tx, nx, hx, dx, Lx, x_start)

    return Grid1D{typeof(Tx)}(nx, hx, dx, Lx, xF, xC)
end

"""
    struct Grid2D{Tx<:AbstractTopology, Ty<:AbstractTopology} <: AbstractGrid

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
    xF::AbstractArray
    "Cell centres in x-direction"
    xC::AbstractArray
    "Cell faces in y-direction"
    yF::AbstractArray
    "Cell centres in y-direction"
    yC::AbstractArray
end

"""
    Grid2D(Tx, Ty, nx, ny, x_start, x_end, y_start, y_end; hx=0, hy=0)

Construct a two-dimensional staggered `grid` on domain `(x, y) ∈ [x_start, x_end] x [y_start, y_end]`
with topologies `{Tx, Ty}`, with `{nx, ny}` interior grid points, and `{hx, hy}` halo points.
"""
function Grid2D(Tx, Ty, nx, ny, x_start, x_end, y_start, y_end; hx=1, hy=1)

    if hx == 0; throw(error("Number of halo points in x cannot be zero")); end
    if hy == 0; throw(error("Number of halo points in y cannot be zero")); end

    Lx = x_end - x_start
    Ly = y_end - y_start
    
    dx, dy = Lx/nx, Ly/ny
    
    xF = construct_faces(Tx, nx, hx, dx, Lx, x_start)
    xC = construct_centres(Tx, nx, hx, dx, Lx, x_start)
    yF = construct_faces(Ty, ny, hy, dy, Ly, y_start)
    yC = construct_centres(Ty, ny, hy, dy, Ly, y_start)
    
    return Grid2D{typeof(Tx), typeof(Ty)}(nx, ny, hx, hy, dx, dy, Lx, Ly, xF, xC, yF, yC)
end

##
# helper functions

"""
    construct_faces(T::AbstractTopology, n, h, d, L)

Returns the face locations for a dimension with topology `T`, number of grid points `n` and `h`
halo points, grid spacing `d`, and extent `L`.
"""
function construct_faces(T::AbstractTopology, n, h, d, L, start)
    start_face  = start
      end_face  = isa(T, Periodic) ? start + L - d : start + L
    total_faces = isa(T, Periodic) ? n             : n + 1
    
    F = range(start_face - h*d, stop = end_face + h*d, length = total_faces + 2h)
    
    return OffsetArray(F, -h)
end

"""
    construct_centres(T::AbstractTopology, n, h, d, L, start_center)

Returns the cell centers locations for a dimension with topology `T`, number of grid points
`n` and `h` halo points, grid spacing `d`, and extent `L`.
"""
function construct_centres(T::AbstractTopology, n, h, d, L, start)
    start_center = start + d/2
      end_center = start + L - d/2
    
    C = range(start_center - h*d, stop = end_center + h*d, length = n + 2h)
    
    return OffsetArray(C, -h)
end

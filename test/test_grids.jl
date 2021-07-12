function test_dx(grid, dx)
    return grid.dx ≈ dx
end

function test_dy(grid::Grid2D, dy)
    return grid.dy ≈ dy
end

function test_xF(grid, xF)
    return grid.xF ≈ xF
end

function test_xC(grid, xC)
    return grid.xC ≈ xC
end

function test_yF(grid::Grid2D, yF)
    return grid.yF ≈ yF
end

function test_yC(grid::Grid2D, yC)
    return grid.yC ≈ yC
end

function xdomain_length(grid, Lx)
    return grid.Lx ≈ Lx
end

function ydomain_length(grid::Grid2D, Ly)
    return grid.Ly ≈ Ly
end
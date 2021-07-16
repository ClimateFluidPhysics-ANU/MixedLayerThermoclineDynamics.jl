function test_interpolation(func, grid::Grid1D, rtol_interpolate)
  
    input_field = Field(Centre, func.(grid.xC), grid) 
    output_field = Field(Face, 0*grid.xF, grid)
    
    interpolate!(output_field, input_field)
    
    return isapprox(output_field.data, func.(grid.xF), rtol=rtol_interpolate)
end
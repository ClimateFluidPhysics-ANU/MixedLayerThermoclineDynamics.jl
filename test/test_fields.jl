function test_halos(input_field, array)
    return isapprox(input_field.data, array, rtol = 1e-1)
end
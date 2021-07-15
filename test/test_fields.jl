function location_centre(h1D)
    return typeof(h1D) == Field1D{Centre}
end

function location_face(u1D)
    return typeof(u1D) == Field1D{Face}
end

function test_grid(field_value, grid)
    centre_bounds = false
    face_bounds = false
    length_eval =false

    centre_bounds = (field_value.grid.xC ≈ grid.xC)
    face_bounds = (field_value.grid.xF ≈ grid.xF)

    length_eval = (field_value.grid.Lx ≈ grid.Lx)

    if(centre_bounds & face_bounds & length_eval)
        return true
    else
        return false
    end
end

function test_data(field_value, data)
    return field_value.data ≈ data
end
using MixedLayerThermoclineDynamics

function test_𝐼x(actual::Field1D, test::Field1D, input::Field1D)
    MixedLayerThermoclineDynamics.𝐼x!(test, input)
    return isapprox(test.data, actual.data; rtol = 1e-2)
end

function test_𝐼x(actual::Field2D, test::Field2D, input::Field2D)
    MixedLayerThermoclineDynamics.𝐼x!(test, input)
    return isapprox(test.data, actual.data; rtol = 1e-2)
    #return test.data ≈ actual.data
end

function test_𝐼y(actual::Field2D, test::Field2D, input::Field2D)
    MixedLayerThermoclineDynamics.𝐼y!(test, input)
    return isapprox(test.data, actual.data; rtol = 1e-2)
    #return test.data ≈ actual.data
end
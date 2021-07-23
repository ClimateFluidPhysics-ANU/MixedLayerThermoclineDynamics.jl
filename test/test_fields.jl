using MixedLayerThermoclineDynamics

function test_𝐼x(actual::Field1D, test::Field1D, input::Field1D)
    MixedLayerThermoclineDynamics.𝐼x!(test, input)
    return isapprox(test.data, actual.data; rtol = 1e-2)
end

function test_𝐼x(actual::Field2D, test::Field2D, input::Field2D)
    MixedLayerThermoclineDynamics.𝐼x!(test, input)
    return isapprox(test.data, actual.data; rtol = 1e-2)
end

function test_𝐼y(actual::Field2D, test::Field2D, input::Field2D)
    MixedLayerThermoclineDynamics.𝐼y!(test, input)
    return isapprox(test.data, actual.data; rtol = 1e-2)
end

function test_∂x(actual::Field1D, test::Field1D, input::Field1D)
    MixedLayerThermoclineDynamics.∂x!(test, input)
    return isapprox(test.data, actual.data; rtol = 1e-2)
end

function test_∂x(actual::Field2D, test::Field2D, input::Field2D)
    MixedLayerThermoclineDynamics.∂x!(test, input)
    return isapprox(test.data, actual.data, rtol = 1e-2)
end

function test_∂y(actual::Field2D, test::Field2D, input::Field2D)
    MixedLayerThermoclineDynamics.∂y!(test, input)
    return isapprox(test.data, actual.data, rtol = 1e-2)
end

module MixedLayerThermoclineDynamics

export hello, domath

"""
    hello(who::String)

Return "Hello, `who`".
"""
hello(who::String) = "Hello, $who"

"""
    domath(x::Number)

Return `x + 5`.
"""
domath(x::Number) = x + 5

"""
    compute_square(x)

Compute the square of something.
"""
compute_square(x) = x^2

end # module

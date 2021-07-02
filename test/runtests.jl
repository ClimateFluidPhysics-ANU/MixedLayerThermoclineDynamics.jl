using Test, MixedLayerThermoclineDynamics

@test hello("Julia") == "Hello, Julia"
@test domath(2.0) ≈ 7.0
@test MixedLayerThermoclineDynamics.compute_square(2.0) ≈ 4.0

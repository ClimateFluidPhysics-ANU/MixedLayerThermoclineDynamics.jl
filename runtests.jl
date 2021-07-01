using Test, MixedLayerThermoclineDynamics

@test hello("Julia") == "Hello, Julia"
@test domath(2.0) ≈ 7.0
@test ToyModel.compute_square(3) == 9

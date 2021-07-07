using Documenter, MixedLayerThermoclineDynamics

makedocs(modules = [MixedLayerThermoclineDynamics], sitename = "MixedLayerThermoclineDynamics.jl")

deploydocs(        repo = "github.com/ClimateFluidPhysics-ANU/MixedLayerThermoclineDynamics.jl.git",
               versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
           push_preview = true
)

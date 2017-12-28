# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Documenter, FEMBasis

cp("../README.md", "src/index.md"; remove_destination=true)

makedocs(modules=[FEMBasis],
         format = :html,
         checkdocs = :all,
         sitename = "FEMBasis.jl",
         pages = [
                  "Introduction" => "index.md",
                  "API documentation" => "api.md",
                 ])

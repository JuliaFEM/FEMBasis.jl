# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBasis.jl/blob/master/LICENSE

using Documenter, FEMBasis

makedocs(modules=[FEMBasis],
         format = :html,
         sitename = "FEMBasis",
         pages = [
                  "FEMBasis" => "index.md",
                 ])

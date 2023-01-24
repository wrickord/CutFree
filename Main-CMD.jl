import Pkg

Pkg.activate("cutfree-venv")
Pkg.instantiate()

include("CutFree.jl")
include("CutFreeRL.jl")

"""
main()

    Description:
    Run simulations using your command line.
"""
function CutFreeCMD(starting_oligo::String, restriction_sites::Vector{String}, min_blocks::Int64, increase_diversity::Bool)
    CutFree_out = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    print("Time: ", CutFree_out.time, " seconds")
end

starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
restriction_sites = ["GGTCTC", "GGCCGG"]
min_blocks = 1
increase_diversity = true

CutFreeCMD(starting_oligo, restriction_sites, min_blocks, increase_diversity)
import Pkg

Pkg.activate("cutfree-venv")
Pkg.instantiate()

include("CutFree.jl")
include("CutFreeRL.jl")

"""
main()

    Description:
    Run simulations.
"""
function main()
    starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
    restriction_sites = ["GGTCTC", "GGCCGG"]
    min_blocks = 1
    increase_diversity = true

    CutFree_out = @timed CutFree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    println("Time: ", CutFree_out.time, " seconds")
    
    # CutFreeRL_out = @timed CutFreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    # println("Time: ", CutFreeRL_out.time, " seconds")
end

main()
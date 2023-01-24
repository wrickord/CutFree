import Pkg

Pkg.activate("cutfree-venv")
Pkg.instantiate()

include("CutFree.jl")
include("CutFreeRL.jl")

"""
main()

    Description:
    Run CutFree simulations.
"""
function main()
    starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
    restriction_sites = ["GGTCTC", "GGCCGG", "GATGTC", "ATCGAGCCGG", "TCAGTCTC", "GGCG"]
    min_blocks = 1
    increase_diversity = true

    cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    println("Time: ", cutfree_output.time, " seconds")
    println("Randomer: ", cutfree_output.value)
    
    # cutfreeRL_output = @timed CutFreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    # println("Time: ", cutfreeRL_output.time, " seconds")
    # println("Randomer: ", cutfreeRL_output.value)
end

main()
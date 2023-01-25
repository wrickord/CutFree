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
    restriction_sites = ["GGTCTC", "GGCCGG"]
    min_blocks = 1
    increase_diversity = true

    cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    cutfreeRL_output = @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)

    println("\nCutFree Time: ", cutfree_output.time, " seconds")
    println("CutFree Randomer: ", cutfree_output.value)
    println("CutFree Degeneracy: ", degeneracy(cutfree_output.value))

    println("\nCutFreeRL Time: ", cutfreeRL_output.time, " seconds")
    println("CutFreeRL Randomer: ", cutfreeRL_output.value)
    println("CutFreeRL Degeneracy: ", degeneracy(cutfreeRL_output.value))

end

main()
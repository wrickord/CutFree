import Pkg

Pkg.activate("cutfree-venv")
Pkg.instantiate()

using Suppressor

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

    # compile functions while supressing outputs
    @suppress_out begin
        @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
        @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    end

    cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    cutfreeRL_output = @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)

    println("\nCutFree Randomer: ")
    print_oligo_block(cutfree_output.value)
    println("CutFree Degeneracy: ", degeneracy(cutfree_output.value))
    println("CutFree Time: ", cutfree_output.time, " seconds")

    println("\nCutFreeRL Randomer: ")
    print_oligo_block(String(cutfreeRL_output.value))
    println("CutFreeRL Degeneracy: ", degeneracy(String(cutfreeRL_output.value)))
    println("CutFreeRL Time: ", cutfreeRL_output.time, " seconds")
end

main()
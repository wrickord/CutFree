import Pkg

Pkg.activate("cutfree-venv")
Pkg.instantiate()

using Distributed, Suppressor

@everywhere include("cutfree-algorithms/CutFree.jl")
@everywhere include("cutfree-algorithms/CutFreeRL.jl")

"""
parse_commandline()

    Description:
    Retrieve information from command line.
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--starting_oligo"
            help = "starting options for oligo"
            arg_type = String
            default = "NNNNNNNNNNNNNNNNNNNN"
        "--restriction_sites"
            help = "restriction enzyme sequences to be restricted"
            arg_type = String
            default = "GGTCTC,GGCCGG"
        "--min_blocks"
            help = "minimum number of desired blocks"
            arg_type = Int64
            default = 1
        "--increase_diversity"
            help = "increase the diversity of the randomer while maintaining its degeneracy"
            arg_type = Bool
            default = true
    end
    return parse_args(s)
end

"""
main()

    Description:
    Run simulations using your command line.
"""
function main()
    @show parsed_args = parse_commandline()

    starting_oligo = parsed_args["starting_oligo"]
    restriction_sites = split(parsed_args["restriction_sites"], ",")
    min_blocks = parsed_args["min_blocks"]
    increase_diversity = parsed_args["increase_diversity"]

    # for manual entry of inputs
    #=
    starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
    restriction_sites = ["GGTCTC", "GGCCGG"]
    min_blocks = 1
    increase_diversity = true
    =#

    # compile functions while supressing outputs
    @suppress_out begin
        @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
        @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    end

    cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    cutfreeRL_output = @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)

    println("\nCutFree Randomer: ")
    print_oligo_block(String(cutfree_output.value))
    println("CutFree Degeneracy: ", get_degeneracy(String(cutfree_output.value)))
    println("CutFree Time: ", cutfree_output.time, " seconds")

    println("\nCutFreeRL Randomer: ")
    print_oligo_block(String(cutfreeRL_output.value))
    println("CutFreeRL Degeneracy: ", get_degeneracy(String(cutfreeRL_output.value)))
    println("CutFreeRL Time: ", cutfreeRL_output.time, " seconds")
end

main()
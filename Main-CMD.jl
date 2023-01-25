import Pkg

Pkg.activate("cutfree-venv")
Pkg.instantiate()

include("CutFree.jl")
include("CutFreeRL.jl")

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

    cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    cutfreeRL_output = @timed CutFreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)

    println("\nCutFree Time: ", cutfree_output.time, " seconds")
    println("CutFree Randomer: ", cutfree_output.value)
    println("CutFree Degeneracy: ", degeneracy(cutfree_output.value))

    println("\nCutFreeRL Time: ", cutfreeRL_output.time, " seconds")
    println("CutFreeRL Randomer: ", cutfreeRL_output.value)
    println("CutFreeRL Degeneracy: ", degeneracy(cutfreeRL_output.value))
end

main()
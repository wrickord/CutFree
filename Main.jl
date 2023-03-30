import Pkg

setup_environment = true

if setup_environment
    Pkg.activate("cutfree-venv")
    Pkg.instantiate()
end

using ArgParse, PyCall, Conda

include("./cutfree-algorithms/CutFree.jl")
include("./cutfree-algorithms/CutFreeRL.jl")

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

    py"""
    import pickle 

    def get_algorithm(starting_oligo, restriction_sites):
        with open("./cutfree-models/cutfree_model.pkl", "rb") as f:
            model = pickle.load(f)
        return model.predict([[len(starting_oligo), len(restriction_sites), len(restriction_sites[0])]])
    """

    algorithm_choice = py"get_algorithm"(starting_oligo, restriction_sites)

    if algorithm_choice[1] == 0 || min_blocks > 1
        println("Optimizing...")
        algorithm_name = "CutFree"
        cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    elseif algorithm_choice[1] == 1
        println("Optimizing...")
        algorithm_name = "CutFreeRL"
        cutfree_output = cutfree_output = @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    end

    # compile functions while supressing outputs
    # @suppress_out begin
    #     @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    #     @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    # end

    # cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    # cutfreeRL_output = @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)

    println("\nAlgorithm: ", algorithm_name)
    println("Randomer: ", String(cutfree_output.value))
    println("Degeneracy: ", get_degeneracy(String(cutfree_output.value)))
    println("Time: ", cutfree_output.time, " seconds")
    
    return [algorithm_name, String(cutfree_output.value), get_degeneracy(String(cutfree_output.value)), cutfree_output.time]
end

main()
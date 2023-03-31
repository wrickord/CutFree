import Pkg

setup_environment = false

if setup_environment
    Pkg.activate("cutfree-venv")
    Pkg.instantiate()
end

using ArgParse, PyCall, Conda, Suppressor

if setup_environment
    Conda.add("scikit-learn")
end

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
            default = false
        "--algorithm"
            help = "force solve with a specific algorithm (CutFree or CutFreeRL)"
            arg_type = String
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
    forced_algorithm = parsed_args["algorithm"]

    if forced_algorithm == "CutFree"
        println("Optimizing...")
        algorithm_name = "CutFree"
        cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    elseif forced_algorithm == "CutFreeRL"
        println("Optimizing...")
        algorithm_name = "CutFreeRL"
        cutfree_output = cutfree_output = @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    else
        py"""
        import pickle 

        def get_algorithm_classification(starting_oligo, restriction_sites):
            with open("./cutfree-models/cutfree_algorithm_classification_model.pkl", "rb") as f:
                model = pickle.load(f)
            return model.predict([[len(starting_oligo), len(restriction_sites), len(restriction_sites[0])]])
        """

        choose_RL = 0
        algorithm_choice = py"get_algorithm_classification"(starting_oligo, restriction_sites)

        if algorithm_choice[1] == 0 || min_blocks > 1 || increase_diversity == true
            py"""
            import pickle 

            def accept_classification(starting_oligo, restriction_sites):
                with open("./cutfree-models/cutfree_runtime_classification_model.pkl", "rb") as f:
                    model = pickle.load(f)
                return model.predict([[len(starting_oligo), len(restriction_sites), len(restriction_sites[0])]])
            """
            println("Checking if CutFree is a good choice...")
            accept_choice = py"accept_classification"(starting_oligo, restriction_sites)

            if accept_choice[1] == 1
                println("Accepted!")
                println("Optimizing...")
                algorithm_name = "CutFree"
                cutfree_output = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
            else
                println("Rejected! Using CutFreeRL instead...")
                choose_RL = 1
            end
        end

        if algorithm_choice[1] == 1 || choose_RL == 1
            println("Optimizing...")
            algorithm_name = "CutFreeRL"
            cutfree_output = cutfree_output = @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
        end
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
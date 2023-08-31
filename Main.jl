# import packages
import Pkg

setup_environment = true
if setup_environment
    Pkg.activate("venv")
    Pkg.instantiate()

    # run(`conda create -n cutfree python=3.9 conda tensorflow=2.12 transformers=4.18`)
    
    USERNAME = "wrick"
    ENV["CONDA_JL_HOME"] = "C:/Users/$USERNAME/anaconda3/envs/cutfree"
    Pkg.build("Conda")

    # ENV["PYTHON"] = ""
    # Pkg.build("PyCall")
end

using ArgParse, PyCall, Conda
Conda.add(["transformers", "tensorflow", "numpy"])

# import julia files
include("./cutfree-algorithms/CutFree.jl")
include("./cutfree-algorithms/CutFreeRL.jl")


# helper functions
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
            help = "increase randomer diversity while maintaining degeneracy"
            arg_type = Bool
            default = false
        "--algorithm"
            help = "force solve with specific algorithm (CutFree or CutFreeRL)"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end


"""
main()

    Description:
    Run CutFree (or CutFreeRL) algorithm on command line arguments.
"""
function main(;
    starting_oligo::String="NNNNNNNNNNNNNNNNNNNN",
    restriction_sites::Any=["GGTCTC"],
    min_blocks::Int64=1,
    increase_diversity::Bool=false,
    forced_algorithm::String="auto",
    command_line::Bool=true
)

    if command_line
        @show parsed_args = parse_commandline()

        starting_oligo = parsed_args["starting_oligo"]
        restriction_sites = split(parsed_args["restriction_sites"], ",")
        min_blocks = parsed_args["min_blocks"]
        increase_diversity = parsed_args["increase_diversity"]
        forced_algorithm = lowercase(parsed_args["algorithm"])
    end

    if forced_algorithm == "cutfree"
        println("Optimizing...")
        algorithm_name = "CutFree"
        cutfree_output = @timed cutfree(
            starting_oligo, 
            restriction_sites, 
            min_blocks, 
            increase_diversity
        )
    elseif forced_algorithm == "cutfreerl"
        println("Optimizing...")
        algorithm_name = "CutFreeRL"
        cutfree_output = @timed cutfreeRL(
            starting_oligo, 
            restriction_sites, 
            simulate=simulate_random,
            nsims=1000
        )
    else
        oligo_length = length(starting_oligo)
        restriction_sites = "GGTCTC"

        py"""
        import tensorflow as tf
        import numpy as np
        from transformers import DistilBertTokenizer, TFDistilBertModel, \
            AdamWeightDecay
        
        def choose_algorithm(oligo_length, restriction_sites):
            model = tf.keras.models.load_model(
                "algorithm-classifier/models/AlgorithmClassifier-V2/" 
                + "checkpoints/best_model.h5",
                custom_objects={
                    "TFDistilBertModel": TFDistilBertModel, 
                    "AdamWeightDecay": AdamWeightDecay
                }
            )

            model_name = "distilbert-base-uncased"
            tokenizer = DistilBertTokenizer.from_pretrained(model_name)

            text_attributes = tf.convert_to_tensor(
                ([oligo_length])
            )

            tokens = tokenizer.batch_encode_plus(
                [restriction_sites], 
                padding="max_length",
                truncation=True,
                max_length=100,
                return_tensors="tf"
            )

            input_ids = tokens["input_ids"]
            attention_mask = tokens["attention_mask"]
            texts = (input_ids, attention_mask, text_attributes)

            return np.argmax(model.predict(texts), axis=1)
        """

        algorithm_choice = py"choose_algorithm"(
            oligo_length, 
            restriction_sites
        )

        algorithm_choice = [1]

        if (
            algorithm_choice[1] == 0 
            || min_blocks > 1 
            || increase_diversity == true
        )
            println("Optimizing...")
            algorithm_name = "CutFree"
            cutfree_output = @timed cutfree(
                starting_oligo, 
                restriction_sites, 
                min_blocks, 
                increase_diversity
            )
        else
            println("Optimizing...")
            algorithm_name = "CutFreeRL"
            cutfree_output = @timed cutfreeRL(
                starting_oligo, 
                restriction_sites, 
                simulate=simulate_random, 
                nsims=1000
            )
        end
    end

    println("\nAlgorithm: ", algorithm_name)
    println("Randomer: ", String(cutfree_output.value))
    println("Degeneracy: ", get_degeneracy(String(cutfree_output.value)))
    println("Time: ", cutfree_output.time, " seconds")
    
    return [
        algorithm_name, 
        String(cutfree_output.value), 
        get_degeneracy(String(cutfree_output.value)), 
        cutfree_output.time
    ]
end


# test case
starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
restriction_sites = ["GGTCTC"] # BsaI

main(
    starting_oligo=starting_oligo, 
    restriction_sites=restriction_sites,
    command_line=false
)
import Pkg

include("CutFree.jl")
include("CutFreeRL.jl")

"""
setup_packages(setup)

    Description:
    Setup packages required for CutFree and CutFreeRL.
"""
function setup_packages(setup=false)
    packages = ["Clp", "Gurobi", "GLPK", "JuMP", "DataStructures", 
        "NamedArrays", "DataFrames", "CSV", "BioSequences", "XLSX", "ArgParse"]
    if setup
        for package in packages
            Pkg.add(package)
            Pkg.build(package)
            Pkg.update(package)
        end
    end
end

"""
main()

    Description:
    Run simulations.
"""
function main()
    setup_packages(false) # change to true for your first run

    starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
    restriction_sites = ["GGTCTC", "GGCCGG"]
    min_blocks = 1
    increase_diversity = true

    cutfree_out = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    print("Time: ", cutfree_out.time, " seconds")
    
    # cutfreeRL_out = @timed CutFreeRL(starting_oligo, restriction_sites)
    # print("Time: ", cutfreeRL_out.time, " seconds")
end

main()
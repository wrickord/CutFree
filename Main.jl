import Pkg

"""
main()

    Description:
    Run simulations.
"""
function main()
    Pkg.activate("cutfree-venv")
    Pkg.instantiate()
    
    include("CutFree.jl")
    include("CutFreeRL.jl")

    starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
    restriction_sites = ["GGTCTC", "GGCCGG"]
    min_blocks = 1
    increase_diversity = true

    CutFree_out = @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
    print("Time: ", CutFree_out.time, " seconds")
    
    CutFreeRL_out = @timed CutFreeRL(starting_oligo, restriction_sites)
    print("Time: ", CutFreeRL_out.time, " seconds")
end

main()






# old function, no longer needed
"""
setup_packages(setup)

    Description:
    Setup packages required for CutFree and CutFreeRL.
    Not needed if you activate the virtual enviornment.
"""
function setup_packages(setup=true)
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
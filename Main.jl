setup = false # change to true for your first run

if setup
    import Pkg
    
    packages = ["Clp", "Gurobi", "GLPK", "JuMP", "DataStructures", 
    "NamedArrays", "DataFrames", "CSV", "BioSequences", "XLSX", "ArgParse"]
    for package in packages
        Pkg.add(package)
        Pkg.build(package)
        Pkg.update(package)
    end
end

include("CutFree.jl")
# include("CutFreeRL.jl")

function main()
    starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
    restriction_sites = ["GGTCTC", "GGCCGG"]

    cutfree(starting_oligo, restriction_sites)
end

main()
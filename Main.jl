if update == true
    import Pkg
    
    Pkg.add("Clp")
    Pkg.add("Gurobi")
    Pkg.add("GLPK")
    Pkg.add("JuMP")
    Pkg.add("DataStructures")
    Pkg.add("NamedArrays")

    Pkg.build("Clp")
    Pkg.build("Gurobi")
    Pkg.build("GLPK")
    Pkg.build("JuMP")
    Pkg.build("DataStructures")
    Pkg.build("NamedArrays")

    Pkg.update("Clp")
    Pkg.update("Gurobi")
    Pkg.update("GLPK")
    Pkg.update("JuMP")
    Pkg.update("DataStructures")
    Pkg.update("NamedArrays")

starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
restriction_sites = ["GGTCTC", "GGCCGG"]

CutFree.CutFree(starting_oligo, restriction_sites)
# CutFree

CutFree is a tool created to design pools of DNA barcodes compatible with restriction enzymes.

Feel free to check out CutFree through the following link: https://jensenlab.shinyapps.io/cutfree/

Make sure to add and build the correct packages before running the program:

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

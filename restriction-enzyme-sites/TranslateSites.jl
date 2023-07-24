import Pkg

Pkg.activate("cutfree-venv")
Pkg.instantiate()

using CSV, DataFrames, StatsBase, Suppressor

"""
translate_blocking_sites(file)

    Description:
    Creates dictionary to check if input is a valid enzyme, and if so, returns the sequence of the blocking site.
"""
function translate_blocking_sites(enzyme_input)
    df = CSV.read("cutfree-blockingsites\\rebase_data.csv", DataFrame)

    enzymes = []
    for enz in df[:,"prototype"]
        push!(enzymes, String(enz))
    end

    blocking_sites = []
    for seq in df[:,"base_sequence"]
        push!(blocking_sites, String(seq))
    end

    translation = Dict{String, String}()
    for i in eachindex(enzymes)
        enzyme = lowercase(enzymes[i])
        blocking_site = blocking_sites[i]
        translation[enzyme] = blocking_site
    end

    try
        return translation[lowercase(enzyme_input)]
    catch
        return "INVALID"
    end
end

"""
main()

    Description:
    Translate input to valid sequence.
"""
function main()
    print(translate_blocking_sites("hahaha"))
end

main()
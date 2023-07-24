import Pkg

Pkg.activate("venv")
Pkg.instantiate()

using CSV, DataFrames, StatsBase, Suppressor

include("../cutfree-algorithms/CutFree.jl")
include("../cutfree-algorithms/CutFreeRL.jl")

"""
read_restriction_sites(file)

    Description:
    Filter restriction enzymes data for the commericially available restriction enzymes.
    Adjust results based on length and palindromic arguments.
"""
function read_restriction_sites(file="../restriction-enzyme-sites/rebase_data.csv", enzyme_length=6, is_palindromic=false)
    df = CSV.read(file, DataFrame)
    commercial_df = df[df[:,"commercial_source"].!="NA",:] # clean to only include commercially available restriction enzymes
    re_df = commercial_df[commercial_df[:,"length"].==enzyme_length,:] # clean based on length specification
    # re_df = temp_df[temp_df[:,"is_palindrome"].==is_palindromic,:] # clean based on palindromic specification

    blocking_sites = []
    for seq in re_df[:,"base_sequence"]
        push!(blocking_sites, String(seq))
    end

    return blocking_sites
end

"""
get_oligos(oligo_min_length, oligo_max_length)

    Description:
    Get all possible oligo lengths for simulation.
    All oligos will be repeats of "N".
"""
function get_oligos(oligo_min_length=6, oligo_max_length=20)
    oligos = []
    for i in oligo_min_length:oligo_max_length
        seq = repeat("N", i)
        push!(oligos, seq)
    end

    return oligos
end

"""
write_output(oligo, sites)

    Description:
    Record all simulation results and output to CSV file.
"""
function write_output(oligo, sites, clear_csv)
    cutfree_output = @timed cutfree(oligo, sites, 1, false)
    cutfreeRL_output = @timed cutfreeRL(oligo, sites, simulate=simulate_random, nsims=1000)

    file = "runtime_data.csv"
    output = DataFrame(
        Oligo = oligo,
        OligoLength = length(oligo),
        Sites = [sites],
        TotalSites = length(sites),
        SiteLength = length(sites[1]),
        CutFree_Randomer = String(cutfree_output.value),
        CutFree_Degeneracy = get_degeneracy(String(cutfree_output.value)),
        CutFree_Time = cutfree_output.time,
        CutFreeRL_Randomer = String(cutfreeRL_output.value),
        CutFreeRL_Degeneracy = get_degeneracy(String(cutfreeRL_output.value)),
        CutFreeRL_Time = cutfreeRL_output.time,
    )

    if clear_csv
        CSV.write(file, output, append=false)
    else
        CSV.write(file, output, append=true)
    end
end

"""
read_input(oligo, sites)

    Description:
    Read all simulation data from CSV.
"""
function read_input(file="runtime_data.csv")
    df = CSV.read(file, DataFrame)
    return df
end

"""
main()

    Description:
    Run simulations.
"""
function main()
    # compile functions while supressing outputs
    starting_oligo = "NNNNNNNNNNNNNNNNNNNN"
    restriction_sites = ["GGTCTC", "GGCCGG"]
    min_blocks = 1
    increase_diversity = false
    @suppress_out begin
        @timed cutfree(starting_oligo, restriction_sites, min_blocks, increase_diversity)
        @timed cutfreeRL(starting_oligo, restriction_sites, simulate=simulate_random, nsims=1000)
    end

    # run simulations
    clear_csv = false
    num_samples = 1
    oligos = get_oligos(6, 40)
    for _ in 1:num_samples
        for oligo in oligos
            for site_length in 4:8
                for site_num in 1:10
                    blocking_sites = unique(read_restriction_sites("../restirction-enzyme-sites/rebase_data.csv", site_length, true))
                    write_output(oligo, sample(blocking_sites, site_num, replace=false), clear_csv)
                    clear_csv = false
                end
            end
        end
    end
end

main()
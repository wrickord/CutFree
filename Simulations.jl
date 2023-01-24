"""
read_restriction_sites(file)

    Description:
    Filter restriction enzymes data for the commericially available restriction enzymes.
    Adjust results based on length and palindromic arguments.
"""
function read_restriction_sites(file="rebase_data.csv", enzyme_length=6, is_palindromic=true)
    df = DataFrame(CSV.File(file))
    commercial_df = df[df[:,"commercial_source"].!="NA",:] # clean to only include commercially available restriction enzymes
    temp_df = commercial_df[commercial_df[:,"length"].==enzyme_length,:] # clean based on length specification
    re_df = temp_df[temp_df[:,"is_palindrome"].==is_palindromic,:] # clean based on palendromic specification

    blocking_sites = []
    for seq in re_df[:,"base_sequence"]
        push!(blocking_sites, seq)
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
main()

    Description:
    Run simulations.
"""
function main()
    blocking_sites = read_restriction_sites("rebase_data.csv", 6, true)
    oligos = get_oligos(6, 20)

    print(oligos)
end

main()
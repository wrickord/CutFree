# import packages
import Base: *, sort
using Random, Statistics, DataFrames, CSV 
using BioSequences, XLSX, ArgParse


# helper functions
"""
*(seq, base)

    Description:
    Overload the concatenation operator to combine a sequence and a single 
    base, i.e. dna"AGCGTGC" * DNA_T.
"""
function *(seq::LongSequence{DNAAlphabet{4}}, base::DNA)
    seq * LongDNA{4}([base])
end


"""
sort_dna(seq; rev)

    Description:
    Create sort function to order DNA codes by their degeneracy 
    (i.e., A,G,C,T,M,R,W,S,Y,K,V,H,D,B,N).
"""
function sort_dna(seq::LongSequence{DNAAlphabet{4}}; rev::Bool=false)
    sort(collect(seq), by=get_degeneracy_rl, rev=rev)
end


"""
isvalid(seq, sites)

    Description:
    Check if any of the sites appear in the sequence.
"""
function isvalid(seq, sites)
    for site in sites
        if occursin(ExactSearchQuery(site,iscompatible), seq)
            return false
        end
    end

    return true
end


"""
find_valid_randomer(prefix, bases, sites)

    Description:
    Randomly select bases to add the prefix, ensuring no site appears in 
    the randomer.
"""
function find_valid_randomer(prefix, bases, sites)
    n = length(bases)
    pre_n = length(prefix)
    randomer = prefix * (dna"-" ^ n)
    for i = 1:n
        candidates = shuffle(bases[i])
        for cand in candidates
            randomer[pre_n+i] = cand
            if isvalid(randomer[1:pre_n+i], sites)
                break
            else
                randomer[pre_n+i] = DNA_Gap
            end
        end
        if randomer[pre_n+i] == DNA_Gap
            break
        end
    end
    return randomer
end


"""
get_degeneracy_rl(base; uselog)

    Description:
    Calculate the degeneracy of a single base. The degeneracy is the  
    number of bases in the code, so:
        degeneracy(DNA_N) == 4
    and
        degeneracy(DNA_A) == 1

    We score using the log of the degeneracy. Because degeneracy(DNA_Gap) = 0, 
    we approximate the log transformation as 2^-100. 
"""
function get_degeneracy_rl(base::DNA; uselog=true)
    deg = 0.0
    if base == DNA_N
        deg = 4.0
    elseif base == DNA_A || base == DNA_C || base == DNA_G || base == DNA_T
        deg = 1.0
    elseif base == DNA_B || base == DNA_D || base == DNA_H || base == DNA_V
        deg = 3.0
    elseif (base == DNA_K || base == DNA_Y || base == DNA_S || base == DNA_W 
            || base == DNA_R || base == DNA_M)
        deg = 2.0
    else # DNA_Gap
        deg = 0.0
    end
    
    if uselog
        if base == DNA_Gap
            deg = 2^-100
        end

        return log(deg)
    else
        return deg
    end
end


"""
get_degeneracy_rl(seq; uselog)

    Description:
    Calculate the degeneracy of a DNA sequence.
"""
function get_degeneracy_rl(seq::LongSequence{DNAAlphabet{4}}; uselog=true)
    if uselog
        return sum(get_degeneracy_rl.(seq, uselog=true))
    else
        return prod(get_degeneracy_rl.(seq; uselog=false))
    end
end


"""
simulate_random(prefix, bases, sites; nsims, horizon)

    Description:
    Run Monte Carlo simulations using a random base policy over a fixed horizon.
"""
function simulate_random(
    prefix, 
    bases, 
    sites; 
    nsims=1000, 
    horizon=length(bases)
)
    degs = zeros(nsims)
    true_horizon = min(horizon, length(bases)) # not beyond randomer length
    for i = 1:nsims
        degs[i] = get_degeneracy_rl(
            find_valid_randomer(prefix, bases[1:true_horizon], sites)
        )
    end

    return mean(degs)
end


"""
simulate_greedy(prefix, bases, sites; horizon)

    Description:
    Run a 1-step greedy lookahead policy.
"""
function simulate_greedy(prefix, bases, sites; horizon=length(bases))
    true_horizon = min(horizon, length(bases)) # not beyond randomer length
    for h = 1:true_horizon
        # iterate through the bases by decending degeneracy, so the greedy 
        # approach is to stop as soon as we've found one.
        found = false
        for base in sort_dna(bases[h], rev=true)
            found = false
            if isvalid(prefix * base, sites)
                prefix *= base
                found = true
                break
            end
        end

        if ~found
            # no valid bases, so stop and leave the DNA_Gap
            break
        end
    end

    return get_degeneracy_rl(prefix)
end


"""
cutfreeRL(sequence, sites; simulate, kwargs)

    Description:
    Run the rollout algorithm to solve the CutFree MDP.

    Arguments:
    - `sequence`: Starting DNA sequence that should be blocked from containing 
    restriction sites. To generate a set of barcodes with the highest diversity, 
    start with a string of N's the length of your oligo. Option is given because 
    some companies restrict which degernate bases are allowed.

    - `sites`: An array of restriciton enzyme recognition sequences to be 
    blocked in the random barcode. 

    - `simulate`: The choice of policy for rollout simulations. 
        simulate_random will use a random rollout policy. 
        simulate_greedy will use a greedy 1 step lookahead policy. 

    - `nsims`: (Optional) The number of rollout simulations per action.
"""
function cutfreeRL(
    sequence::String, 
    sites::Any;
    simulate=simulate_random, 
    kwargs...
)
    # convert strings to biosequence objects 
    sequence = LongSequence{DNAAlphabet{4}}(sequence)
    sites = LongSequence{DNAAlphabet{4}}.(sites)

    allowedbasedict = Dict(
        DNA_A => dna"A", 
        DNA_T => dna"T", 
        DNA_C => dna"C", 
        DNA_G => dna"G",
        DNA_R => dna"AG",
        DNA_Y => dna"CT",
        DNA_S => dna"GC",
        DNA_W => dna"AT",
        DNA_K => dna"GT",
        DNA_M => dna"AC",
        DNA_B => dna"CGTYSKB",
        DNA_D => dna"AGTRWKD",
        DNA_H => dna"ACTYWMH",
        DNA_V => dna"ACGRSMV",
        DNA_N => dna"ACGTRYSWKMBDHVN",
        DNA_Gap => dna"-"
    )

    bases = [dna"-" for _=1:length(sequence)]
    for i in eachindex(sequence)
        bases[i] = allowedbasedict[sequence[i]]
    end 

    n = length(bases)
    randomer = dna"-" ^ n

    # if the restriction site is non-palindromic, we also need to block
    # the reverse complement.
    for i = 1:length(sites)
        if ~ispalindromic(sites[i])
            push!(sites, reverse_complement(sites[i]))
        end
    end

    # find the longest restriction site and use this as the horizon.
    max_len = map(length, sites) |> maximum
    horizon = max_len

    for i = 1:n
        best_base = DNA_Gap
        best_deg = -1.0
        start = max(i - max_len, 1) # only look back max_len to check oligos
        for base in bases[i]
            if ~isvalid(randomer[start:i-1] * base, sites)
                # adding this bases creates a restriction site; skip it
                continue
            end

            deg = simulate(
                randomer[start:i-1] * base, 
                bases[i+1:end], 
                sites; 
                horizon=horizon, 
                kwargs...
            )

            # update use this base if:
            #   1. the mean randomer degeneracy is higher than the 
            #      previous best, OR
            #   2. the randomer degeneracy is tied with the previous best, 
            #      but the new base has higher individual degeneracy
            if deg > best_deg || (
                deg == best_deg 
                && get_degeneracy_rl(base) > get_degeneracy_rl(best_base)
            )
                best_deg = deg
                best_base = base
            end
        end

        if best_base == DNA_Gap
            # sequence has terminated; there are no bases that can be added 
            # without making a restriction site.
            break
        else
            randomer[i] = best_base
        end
    end

    return randomer
end

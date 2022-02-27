IUB_CODES = Dict(
    "A" => ["A"], 
    "C" => ["C"], 
    "G" => ["G"], 
    "T" => ["T"], 
    "R" => ["A", "G"], 
    "Y" => ["C", "T"], 
    "M" => ["A", "C"], 
    "K" => ["G", "T"], 
    "S" => ["C", "G"],
    "W" => ["A", "T"], 
    "H" => ["A", "C", "T"], 
    "B" => ["C", "G", "T"], 
    "V" => ["A", "C", "G"], 
    "D" => ["A", "G", "T"], 
    "N" => ["A", "C", "G", "T"]
)

#= IUB_CODES Description
    Dictionary consisting of keys and values (IUB codes)

    Keys:
        Input to call a value of the Dictionary
    Values:
        Output of the IUB codes used by CutFree
=#

# IUB_CODES Testing
    # println(IUB_CODES["R"][2])
# 

function get_blocking_codes(code, names)
    blocking_codes = Vector()

    # Constant Time: O(1)
    for (key, value) in IUB_CODES
        if filter(x -> x != code, value) == IUB_CODES[key] && names
            push!(blocking_codes, value)
        end

        if filter(x -> x != code, value) == IUB_CODES[key] && !names
            push!(blocking_codes, key)
        end
    end

    return blocking_codes
end

#= get_blocking_codes Description
    Find codes that do not contain a code.

    Example: 
    If code is A, 
        Then blocking codes are C,G,T,Y,K,S,B (all codes without A)
    If names==TRUE, 
        return the codes (character vector of bases);
    otherwise return the code names.
=#

# get_blocking_codes Testing
    # println(get_blocking_codes("A", false))
# 

function str_to_vector(str)
    # Linear Time: O(n)?
    return collect(str)
end

#= str_to_vector Description
    Turn a string into a vector of single characters.

    Example: 
    If input str is ACGT, 
    Then returned vector is ['A', 'C', 'G', 'T']
=#

# str_to_vector Testing
    # println(str_to_vector("CGTACCCGG"))
# 

function make_oligo_block(oligo)
    output = ""
    oligo_bases = str_to_vector(oligo)
    bases = Dict("A" => false, "C" => false, "G" => false, "T" => false)

    # Linear Time: O(n)
    for code in oligo_bases
        # Constant Time: O(1)
        for base in IUB_CODES[string(code)]
            if (!bases[base])
                bases[base] = true
            end
        end
    end

    # Constant Time: O(1)
    for (key, value) in bases
        if value
            output *= key * ": present \n"
        else
            output *= key * ": not present \n"
        end
    end

    return output
end

#= make_oligo_block Description
    creates a text block showing which bases are allowed
    by the oligo. It returns four strings denoting the presence
    of each base.
=#

# make_oligo_block Testing
    # make_oligo_block("ACY")
# 

function print_oligo_block(oligo)
    println(oligo * "\n" * make_oligo_block(oligo))
end

#= print_oligo_block Description
    prints the oligo and the output of make_oligo_block
=#

# print_oligo_block Testing
    # print_oligo_block("AY")
# 




# Overall time test
    #=
    @time begin
        println(IUB_CODES["R"][2])
        println(get_blocking_codes("A", true))
        println(str_to_vector("CGTACCCGG"))
        print_oligo_block("AY")
    end
    #=#
#
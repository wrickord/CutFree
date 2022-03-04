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

function get_blocking_codes(code, names)
    blocking_codes = Vector()

    for (key, value) in IUB_CODES
        if filter(x -> x != code, value) == IUB_CODES[key] && !names
            push!(blocking_codes, value)
        end

        if filter(x -> x != code, value) == IUB_CODES[key] && names
            push!(blocking_codes, key)
        end
    end

    return blocking_codes
end

function str_to_vector(str)
    return string.(collect(str))
end

function vector_to_str(vctr)
    return join(vctr)
end

function make_oligo_block(oligo)
    oligo_bases = str_to_vector(oligo)
    strings = Dict("A" => "", "C" => "", "G" => "", "T" => "")

    i = 0
    for code in oligo_bases
        for base in IUB_CODES[code]
            for (key, value) in strings
                if (key == base)
                    strings[key] *= key
                end
            end
        end

        i += 1
        for (key, value) in strings
            if (sizeof(value) != i)
                strings[key] *= "-"
            end
        end
    end

    return collect(values(strings))
end

function print_oligo_block(oligo)
    return print(oligo * "\n" * make_oligo_block(oligo)[1] * "\n" * make_oligo_block(oligo)[2]
    * "\n" * make_oligo_block(oligo)[3] * "\n" * make_oligo_block(oligo)[4])
end

function xiny(x, y)
    return (x ∈ y)
end

function subcodes(code)
    codes = ""

    for (key, value) in IUB_CODES
        if (issubset(value, IUB_CODES[code]))
            codes *= key
        end
    end
    
    return str_to_vector(codes)
end

function complement(oligo)
    codes = ""
    oligo_bases = str_to_vector(oligo)
    base_comps = Dict(
        "A" => "T", 
        "C" => "G", 
        "G" => "C", 
        "T" => "A",
    )

    for code in oligo_bases
        temp = ""
        for base in IUB_CODES[code]
            temp *= base_comps[base]
        end

        for (key, value) in IUB_CODES
            if (issetequal(str_to_vector(temp), value))
                codes *= key
            end
        end
    end

    return str_to_vector(codes)
end

function reverse_complement(oligo)
    return reverse(complement(oligo))
end

function degeneracy(oligo)
    value = 1
    oligo_bases = str_to_vector(oligo)

    for code in oligo_bases
        value *= length(IUB_CODES[code])
    end

    return value
end





# Overall time test
    # @time begin
    #     println(IUB_CODES["R"][2])
    #     println(get_blocking_codes("A", true))
    #     println(str_to_vector("CGTACCCGG"))
    #     print_oligo_block("AY")
    # end
#
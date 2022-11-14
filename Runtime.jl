include("./CutFreeRL.jl")
using StatsBase

rebase_data=CSV.read("./rebase_data.csv",DataFrame)

rebase_data=filter(x -> x[:length]>=4 && x[:length]<=6,rebase_data)


function runtime_dataset(rebase_data)
    randomer_length_values=collect(6:20) #make a randomer between length 6 and 20
    restriction_site_length_values=[4,5,6]
    restriction_site_values=collect(1:12)
    nsims_values=[100,500,1000]
    algorithms=[simulate_greedy,simulate_random]
    algorithm_dict=Dict(simulate_greedy => "greedy",simulate_random => "RL")
    
    Random.seed!(1111)

    n_experiments=2000
    data=DataFrame(randomer_length=[],site_length=[],n_sites=[],sites=[],nsims=[],algorithm=[],sequence=[],time=[],degeneracy=[])
    for i=1:n_experiments
        randomer_length=sample(randomer_length_values,1)[1]
        site_length=sample(restriction_site_length_values,1)[1]
        n_sites=sample(restriction_site_values,1)[1]
        nsims_value=sample(nsims_values,1)[1]
        algorithm=sample(algorithms,1)[1]
        all_ns=dna"N"^randomer_length
        sites=LongSequence{DNAAlphabet{4}}.(sample(filter(x -> x[:length]==site_length,rebase_data)[:base_sequence],n_sites))
        string_sites=join(String.(sites),",")
        stats=@timed CutFreeRL(all_ns,sites; simulate=algorithm,nsims=nsims_value)
        string_value=String(stats.value)
        deg=degeneracy(stats.value;uselog2=false)
        algorithm_type=algorithm_dict[algorithm]
        exp_data=DataFrame(randomer_length=randomer_length,site_length=site_length,n_sites=n_sites,sites=string_sites,nsims=nsims_value,algorithm=algorithm_type,sequence=string_value,time=stats.time,degeneracy=deg)
        append!(data,exp_data)
        CSV.write("runtime_dataset.csv",data)
    end 
end 


runtime_dataset(rebase_data)





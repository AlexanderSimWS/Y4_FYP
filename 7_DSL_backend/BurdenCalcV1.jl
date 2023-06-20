# -------------------- Burden Calculator -------------------- %
# 1. Implementation of ATP (phosphate bond) cost of biosynthesis of genes/proteins based on gene sequence (from Lynch, M. & Marinov, G.K. (2015) The bioenergetic costs of a gene. Proceedings of the National Academy of Sciences of the United States of America. 112 (51), 15690–15695. doi:https://doi.org/10.1073/pnas.1514974112.)
# 2. Implementation of "nodes" (mRNA to Protein) with specified reaction rates and initial conditions (used Julia macro for easy definition)
# 3. Implementation of function to create ODE system (using Catalyst to create reaction system) by calculating resource demand coefficients (j_term) and k_eff
# 4. "Total cost" of biosynthesis determined by final concentrations * ATP cost per gene
# ----------------------------------------------------------- %

## Import packages
# using Catalyst, Latexify, DifferentialEquations, Plots
using BioSequences

## Calculate transcription cost of 1 transcript using average NT biosynthesis cost
function txcost(x::BioSequences.NucSeq)
	num_nucs = count(iscertain, x)
	tx_cost = num_nucs * (46+2) # 46 P for synthesis, 2 P for polymerisation
	return tx_cost
end

## Calculate translation cost of 1 protein using individual AA biosynthesis costs
### Count AAs
function countaas(x::BioSequences.AASeq)
	num_A = 0
	num_C = 0
	num_D = 0
	num_E = 0
	num_F = 0
	num_G = 0
	num_H = 0
	num_I = 0
	num_K = 0
	num_L = 0
	num_M = 0
	num_N = 0
	num_P = 0
	num_Q = 0
	num_R = 0
	num_S = 0
	num_T = 0
	num_V = 0
	num_W = 0
	num_Y = 0
	
	for aa in x
		if aa == AA_A
			num_A +=1
		elseif aa == AA_C
			num_C +=1
		elseif aa == AA_D
			num_D +=1
		elseif aa == AA_E
			num_E +=1
		elseif aa == AA_F
			num_F +=1
		elseif aa == AA_G
			num_G +=1
		elseif aa == AA_H
			num_H +=1
		elseif aa == AA_I
			num_I +=1
		elseif aa == AA_K
			num_K +=1
		elseif aa == AA_L
			num_L +=1
		elseif aa == AA_M
			num_M +=1
		elseif aa == AA_N
			num_N +=1
		elseif aa == AA_P
			num_P +=1
		elseif aa == AA_Q
			num_Q +=1
		elseif aa == AA_R
			num_R +=1
		elseif aa == AA_S
			num_S +=1
		elseif aa == AA_T
			num_T +=1
		elseif aa == AA_V
			num_V +=1
		elseif aa == AA_W
			num_W +=1
		elseif aa == AA_Y
			num_Y +=1
		end
	end
	return num_A, num_C, num_D, num_E, num_F, num_G, num_H, num_I, num_K, num_L, num_M, num_N, num_P, num_Q, num_R, num_S, num_T, num_V, num_W, num_Y
end
### Map AA counts to costs
function tlsyncost(x::NTuple{20,Int64})
    AA_cost_dict = Dict("A"=>11.7, "C"=>24.7, "D"=>12.7, "E"=>15.3, 
					"F"=>52.0, "G"=>11.7, "H"=>38.3, "I"=>32.3, 
					"K"=>30.3, "L"=>27.3, "M"=>34.3, "N"=>14.7, 
					"P"=>20.3, "Q"=>16.7, "R"=>27.3, "S"=>11.7,
					"T"=>18.7, "V"=>23.3, "W"=>74.3, "Y"=>50.0)
	cost_A = x[1] * AA_cost_dict["A"]
	cost_C = x[2] * AA_cost_dict["C"]
	cost_D = x[3] * AA_cost_dict["D"]
	cost_E = x[4] * AA_cost_dict["E"]
	cost_F = x[5] * AA_cost_dict["F"]
	cost_G = x[6] * AA_cost_dict["G"]
	cost_H = x[7] * AA_cost_dict["H"]
	cost_I = x[8] * AA_cost_dict["I"]
	cost_K = x[9] * AA_cost_dict["K"]
	cost_L = x[10] * AA_cost_dict["L"]
	cost_M = x[11] * AA_cost_dict["M"]
	cost_N = x[12] * AA_cost_dict["N"]
	cost_P = x[13] * AA_cost_dict["P"]
	cost_Q = x[14] * AA_cost_dict["Q"]
	cost_R = x[15] * AA_cost_dict["R"]
	cost_S = x[16] * AA_cost_dict["S"]
	cost_T = x[17] * AA_cost_dict["T"]
	cost_V = x[18] * AA_cost_dict["V"]
	cost_W = x[19] * AA_cost_dict["W"]
	cost_Y = x[20] * AA_cost_dict["Y"]

	tl_syn_cost = cost_A + cost_C + cost_D + cost_E + cost_F + cost_G + cost_H + cost_I + cost_K + cost_L + cost_M + cost_N + cost_P + cost_Q + cost_R + cost_S + cost_T + cost_V + cost_W + cost_Y

	return tl_syn_cost
end
### Cost of translation
function tlcost(x::BioSequences.NucSeq)
	aa_seq = BioSequences.translate(x, code=BioSequences.standard_genetic_code) # translate codons based on e. coli
	num_aas = count(iscertain, aa_seq)
	tl_syn_cost = tlsyncost(countaas(aa_seq))
	tl_cost = (tl_syn_cost + 4*num_aas)
	return tl_cost
end

# ## Circuit Simulation
# ### Node Definition 
# struct NodeStruct
#     name::Symbol
#     k_plus::Float64
#     k_minus::Float64
#     k_cat::Float64
#     init_m::Float64
#     init_x::Float64
# end

# ### Compiled Nodes Definition
# struct CompiledNodes
#     idx::Int8
#     names::Vector{Symbol}
#     k_cat_arr::Vector{Float64}
#     k_m_arr::Vector{Float64}
#     init_m_arr::Vector{Float64}
#     init_x_arr::Vector{Float64}
# end

# ### Node Constructor Macro 
# macro node(sym, k_plus, k_minus, k_cat, init_m, init_x)
#     NodeStruct(sym, k_plus, k_minus, k_cat, init_m, init_x)
# end

# ### Function to Compile Nodes
# function compilenodes(nodes...)
#     idx = 0
#     var_names = []
#     k_cat_arr = []
#     k_m_arr = []
#     init_m_arr = []
#     init_x_arr = []

#     for node in nodes
#         idx +=1 # Node counter
#         k_m = node.k_plus/(node.k_minus + node.k_cat) # Calculate k_m for each node
#         push!(var_names, node.name)
#         push!(k_cat_arr, node.k_cat)
#         push!(k_m_arr, k_m)
#         push!(init_m_arr, node.init_m)
#         push!(init_x_arr, node.init_x)
#     end

#     compiled_nodes = CompiledNodes(idx,var_names,k_cat_arr,k_m_arr,init_m_arr,init_x_arr)
#     return compiled_nodes
# end


# ### Function to Create Reaction System
# function runsystem(compiled_nodes::CompiledNodes, r_tot::Real, t_span::Tuple)
# 	# Use of Catalyst package
#     @parameters t ρ_tot
#     @species m(t)[1:compiled_nodes.idx] x(t)[1:compiled_nodes.idx]
#     rxns = []
#     u0::Vector{Pair{Num, Int64}} = []
#     p = [ρ_tot=>r_tot]

# 	# Maps species names specified in node creation to names in the model
#     legend = Dict()

#     # Calculate resource demand coefficient
#     j_term = 0
#     for i = 1:compiled_nodes.idx
#         j_term += compiled_nodes.k_m_arr[i]*x[i]
#     end

#     # k_eff calculation function
#     k_effs(k_m, k_cat) = k_cat*r_tot*k_m/(1+j_term)

#     # Add reactions with k_eff into reaction vector, append initial conidtions to u0 vector
#     for i = 1:compiled_nodes.idx
#         push!(rxns, Reaction(k_effs(compiled_nodes.k_m_arr[i], compiled_nodes.k_cat_arr[i]), [m[i]], [m[i], x[i]]))
#         append!(u0, [m[i]=>compiled_nodes.init_m_arr[i], x[i]=>compiled_nodes.init_x_arr[i]])
#         get!(legend, "mRNA/Protein "*string(compiled_nodes.names[i]), (m[i], x[i]))
#     end

#     # Create ODE system
#     @named rxn_sys = ReactionSystem(rxns)
#     ode_sys = convert(ODESystem, rxn_sys)
# 	render(latexify(rxn_sys))

#     # Simulation
#     prob = ODEProblem(ode_sys, u0, t_span, p)
#     sol = solve(prob)

#     # Print Legend to REPL
#     println(legend)
#     return sol
# end

# ## Overall Burden (NOT COMPLETE, only takes in 2 dna sequences)
# function calcburden(ode_sol::ODESolution, dna1::BioSequences.NucSeq, dna2::BioSequences.NucSeq)
#     gene_1_cost = txcost(dna1) * ode_sol.u[length(ode_sol)][1]
#     protein_1_cost = tlcost(dna1) * ode_sol.u[length(ode_sol)][2]
#     gene_2_cost = txcost(dna2) * ode_sol.u[length(ode_sol)][3]
#     protein_2_cost = tlcost(dna2) * ode_sol.u[length(ode_sol)][4]

# 	x = ["mRNA 1", "Protein 1", "mRNA 2", "Protein 2"]
# 	y = [gene_1_cost, protein_1_cost, gene_2_cost, protein_2_cost]
# 	bar(x,y, label="ATP Cost")
# end

# # -------------------- Testing Code -------------------- #
# ## Creation of nodes
# node_a = @node(A, 5, 0.2, 10, 1000, 0) # Arguments: name, k+, k-, kcat, initial mRNA conc, initial protein conc.
# node_b = @node(B, 6, 0.3, 15, 100, 0)

# ## Compiles nodes
# compiled_nodes = compilenodes(node_a,node_b)

# ## Creation of ODE system
# system_sol = runsystem(compiled_nodes, 1000, (0,1000)) # Arguments: compiled nodes, total ribosomes, time span

# ## Plot
# plot(system_sol)

## Testing 2
dna_1 = LongDNA{2}("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG") # eGFP sequence
dna_2 = LongDNA{2}("gccctggacaccaactattgcttcagctccacggagaagaactgctgcgtgcgccagctgtacattgacttccgcaaggacctcggctggaagtggatccacgagccgaagggctaccatgccaacttctgcctcggcccgtgcccgtacatttggagcctggacacgcagtacagcaaggtcctggccctgtacaaccagcataacccgggcgcctcggcggcgccgagctgcgtgccgcaggcgctggagccgctgccgatcgtgtactacgtgggccgcaagccgaaggtggagcagctgtccaacatgatcgtgcgctcctgcaagtgcagc") # TGF-B1 sequence
# calcburden(system_sol, dna_1, dna_2)

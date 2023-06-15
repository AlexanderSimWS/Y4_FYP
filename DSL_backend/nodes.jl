# -------------------- Control Architecture Simulator -------------------- %
# Full implementation of control architecture simulation
# Node-based system where activation and repression reactions are stored
# in a modular way
# ------------------------------------------------------------------------ %

# Node
mutable struct Node
    ### Name
    name::Symbol # Alphabet name of node
    ### Transcription
    k_plus_1::Float64 # Binding rate of DNA polymerase to promoter
    k_min_1::Float64 # Dissociation rate of DNA polymerase from promoter
    k_cat_1::Float64 # Elongation constant of mRNA
    ### Translation
    k_plus_2::Float64 # Binding rate of ribosome to RBS
    k_min_2::Float64 # Dissociation rate of ribosome from RBS
    k_cat_2::Float64 # Elongation constant of protein
    ### Degradation
    d_1::Float64 # Degradation rate of mRNA
    ### Initial Conditions
    init_g::Float64 # Initial DNA concentration
end

# Links
abstract type Link end
abstract type TFLink <: Link end
abstract type InducerLink <: Link end
abstract type AnnhilationLink <: Link end

mutable struct Activation <: TFLink
    input_node::Node
    output_node::Node
    k_hat_plus::Float64
    k_hat_min::Float64
    cooperativity::Real
end

mutable struct Repression <: TFLink # TF link
    input_node::Node
    output_node::Node
    k_hat_plus::Float64
    k_hat_min::Float64
    cooperativity::Real
end

mutable struct TFActivation <: InducerLink # Inducer link
    tf_species::Symbol
    node::Node
    k_plus::Float64
    k_min::Float64
    cooperativity::Float64
    init_I::Float64
end

mutable struct TFInactivation <: InducerLink
    tf_species::Symbol
    node::Node
    k_plus::Float64
    k_min::Float64
    cooperativity::Float64
    init_I::Float64
end

mutable struct Annhilation <: AnnhilationLink
    node_1::Node
    node_2::Node
    type::String # Only "mRNA" or "protein"
    k_plus::Float64
    k_min::Float64
end

# Constructors
function createnode(name::Symbol, tx_rates::Vector, tl_rates::Vector, dg_rate::Real, init_g::Real)
    # Convert rate vectors to Float64 type
    float_tx_rates = float(tx_rates)
    float_tl_rates = float(tl_rates)
    float_dg_rate = float(dg_rate)
    float_init_g = float(init_g)

    # Instantiate instance of node 
    basic_node = Node(name, float_tx_rates[1], float_tx_rates[2], float_tx_rates[3],
    float_tl_rates[1], float_tl_rates[2], float_tl_rates[3],
    float_dg_rate,
    float_init_g)
    
    return basic_node
end

function activate(input_node::Node, output_node::Node, rates::Vector, cooperativity::Real)
    # Convert rate vector to Float64 type
    float_rates = float(rates)

    # Instantiate activation link
    link = Activation(input_node, output_node, float_rates[1], float_rates[2], float(cooperativity))
    
    return link
end

function repress(input_node::Node, output_node::Node, rates::Vector, cooperativity::Real)
    # Convert rate vector to Float64 type
    float_rates = float(rates)

    # Instantiate activation link
    link = Repression(input_node, output_node, float_rates[1], float_rates[2], cooperativity)
    
    return link
end

function tfactivation(tf_species::Symbol, node::Node, rates::Vector, cooperativity::Real, init_I::Real)
    # Convert rate vector to Float64 type
    float_rates = float(rates)
    float_init = float(init_I)

    # Instantiate activation link
    link = TFActivation(tf_species, node, float_rates[1], float_rates[2], float(cooperativity), float_init)
    
    return link
end

function tfinactivation(tf_species::Symbol, node::Node, rates::Vector, cooperativity::Real, init_I::Real)
    # Convert rate vector to Float64 type
    float_rates = float(rates)
    float_init = float(init_I)

    # Instantiate activation link
    link = TFInactivation(tf_species, node, float_rates[1], float_rates[2], float(cooperativity), float_init)
    
    return link
end

function annhilation(node_1::Node, node_2::Node, type::String, rates::Vector)
    k_plus = rates[1]
    k_min = rates[2]
    if (type in ["mRNA", "protein"])
        link = Annhilation(node_1, node_2, type, float(k_plus), float(k_min))
    else
        error("Error in annhilation initiator, 'type' field only excepts \"mRNA\" or \"protein\"")
    end
end

# # Macros
# function getreactions(lines::Expr)
#     @assert lines.head == :block
#     name::Symbol = :Z
#     tx_rates::Vector{Float64} = []
#     tl_rates::Vector{Float64} = []
#     d::Real = 0
#     init_g::Real = 0
#     for line in lines.args
#         if (typeof(line) == Expr)
#             if line.args[1] == :~
#                 init_g = line.args[3].args[2]
#                 name = line.args[2]
#             else
#                 (rates, rxn) = line.args
#                 @assert line.head == :tuple
#                 @assert rates.head == :vect
#                 @assert rxn.head in (:-->, :call)

#                 if rxn.head == :call
#                     @assert rates.args[1].head == Symbol("=")
#                     @assert rxn.args[1] == Symbol("<-->")
#                     # fwd_sym = rates.args[1].args[1]
#                     k_plus = rates.args[1].args[2]
#                     # bwd_sym = rates.args[2].args[1]
#                     k_min = rates.args[2].args[2]
#                     # return [k_plus, k_min]
#                     if (:p in rxn.args[2].args)
#                         push!(tx_rates, float(k_plus))
#                         push!(tx_rates, float(k_min))
#                     elseif (:r in rxn.args[2].args)
#                         push!(tl_rates, float(k_plus))
#                         push!(tl_rates, float(k_min))
#                     end
#                 elseif (:∅ == rxn.args[2]) && (:m == rxn.args[1])
#                     d = rates.args[1].args[2]
#                 elseif (rxn.head == :-->) && (length(rxn.args[2].args) == 4)
#                     @assert rxn.args[2].args[1] == :+
#                     @assert rates.args[1].head == Symbol("=")
#                     kcat = rates.args[1].args[2]
#                     if (:C == rxn.args[1])
#                         push!(tx_rates, float(kcat))
#                     elseif (:X == rxn.args[1])
#                         push!(tl_rates, float(kcat))
#                     end
#                 end
#             end
#         end
#     end
#     return createnode(name, tx_rates, tl_rates, d, init_g)
# end

# macro createnode(rxns)
#     getreactions(rxns)
# end

# node_1 = @createnode begin
#     A ~ g = 0.5
#     [k1⁺₁ = 4, k1⁻₁ = 1], g₁ + p <--> C₁
#     [k1θ₁ = 5], C --> g₁ + p + m₁
#     [k2⁺₁ = 4, k2⁻₁ = 1], m₁ + r <--> X₁
#     [k2θ₁ = 5], X₁ --> m₁ + r + x₁
#     [d₁ = 0.04], m₁ --> ∅
# end

# node_2 = @createnode begin
#     B ~ g = 0.5
#     [k1⁺₂ = 4, k1⁻₂ = 1], g₂ + p <--> C₂
#     [k1θ₂ = 5], C₂ --> g₂ + p + m₂
#     [k2⁺₂ = 4, k2⁻₂ = 1], m₂ + r <--> X₂
#     [k2θ₂ = 5], X₂ --> m₂ + r + x₂
#     [d₂ = 0.04], m₂ --> ∅
# end

# node_3 = @createnode begin
#     C ~ g = 0.5
#     [k1⁺₃ = 4, k1⁻₃ = 1], g₃ + p <--> C₃
#     [k1θ₃ = 5], C₃ --> g₃ + p + m₃
#     [k2⁺₃ = 4, k2⁻₃ = 1], m₃ + r <--> X₃
#     [k2θ₃ = 5], X₃ --> m₃ + r + x₃
#     [d₃ = 0.04], m₃ --> ∅
# end

# edge_1 = @activate begin
#     A -> B
#     [k'⁺₂ = 100, k'⁻₂ = 1], g₂ + 3x₁ <--> c₂
# end

# edge_2 = @repress begin
#     B -| A
#     [k'⁺₂ = 100, k'⁻₂ = 1], g₁ + 3x₂ <--> c
# end

# edge_3 = @combi begin
#     A -> C, B -| C
#     [k⁺ₐ₃ = 100, k⁻ₐ₃ = 1], g₃ + x₁ <--> c₃
#     [k⁺ᵣ₃ = 100, k⁻ᵣ₃ = 1], g₃ + 3x₂ <--> c₃
# end

# edge_4 = @annihilate begin
#     A -X- B
#     [k⁺₁₂ = 1000, k⁻₁₂ = 0.1], m₁ + m₂ --> m∅
# end

# edge_5 = @TFactivate begin
#     I₁ -> A
#     [κ⁺₁ = 1000, κ⁺₂ = 1], x⁰₁ + 2 I₁ <--> x⁺₁
# end
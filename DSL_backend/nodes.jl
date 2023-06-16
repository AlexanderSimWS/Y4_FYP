# ------------------------- Node Definitions ------------------------- #
# Definition of Nodes and Edges (links)
# Constructor functions to handle user inputs
# -------------------------------------------------------------------- #
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
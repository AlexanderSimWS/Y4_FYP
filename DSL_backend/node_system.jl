include("nodes.jl")



# Node System
mutable struct NodeSystem
    num_nodes::Int # Total number of nodes
    num_tf_links::NTuple{2, Int} # Number of activation/repression links
    num_inducer_links::NTuple{2, Int} # Number of TF activation/inactivation links (direct activation, direct repression, tf activation, tf inactivation)
    num_ann_links::NTuple{2, Int} # Number of annihilation links (mrna, protein)

    nodes::Vector{Node}
    tf_links::Vector{TFLink}
    inducer_links::Vector{InducerLink}
    ann_links::Vector{AnnhilationLink}

    km_arr::Vector{NTuple{2, Float64}}
    km_hat_arr::Vector{Float64}

    protein_species::Vector{Symbol}
    adj_matrix::Matrix{Int} # matrix for activation/repression of nodes by other nodes/TFs, 1 = activation, 2 = repression
    inducer_matrix::Matrix{Int} # matrix for activation/inactivation of proteins produced by nodes, 1 = activation, 2 = inactivation
    ann_matrix::Matrix{Int} # 1 = annhilation between nodes
end

# Compile Nodes
function compilenodes(nodes...)
    compiled_nodes = Vector{Node}(collect(nodes))
    return compiled_nodes
end

# Compile Links
function compilelinks(links...)
    compiled_links = Vector{Link}(collect(links))
    return compiled_links
end

# Node system constructor
function createnodesystem(nodes::Vector{Node}, links::Vector{Link})
    # Tally Number of Nodes and Links
    ### Number of nodes
    num_nodes = length(nodes)
    ### Number of activated/repressed links
    num_acti = 0
    num_repr = 0
    ### Total number of InducerLinks
    total_inducer_links = 0
    ### Individual number of InducerLinks
    num_tf_acti = 0
    num_tf_inacti = 0
    # Number of annihilation links
    num_mrna_annihilation = 0
    num_protein_annihilation = 0

    # Initialise vectors for TFLinks, InducerLinks, AnnhilationLinks
    tf_links::Vector{TFLink} = []
    inducer_links::Vector{InducerLink} = []
    ann_links::Vector{AnnhilationLink} = []
    
    # Initialise protein species and TF species arrays
    species_arr::Vector{Symbol} = []
    tf_arr::Vector{Symbol} = []

    # Count InducerLinks, append names to tf_arr
    for link in links
        if typeof(link) <: InducerLink
            total_inducer_links +=1
            push!(tf_arr, link.tf_species)
        end
    end

    # Initialise protein species array, adjacency matrix, inducer matrix, annhilation matrix
    adj_matrix::Matrix{Int} = zeros(Int8, num_nodes, num_nodes)
    inducer_matrix::Matrix{Int} = zeros(Int8, total_inducer_links, num_nodes)
    ann_matrix::Matrix{Int} = zeros(Int8, num_nodes, num_nodes)

    # Initialise km array (km_1, km_2)
    km_arr::Vector{NTuple{2, Float64}} = []
    km_hat_arr::Vector{Float64} = []

    # Loop through nodes
    for node in nodes
        # Calculation of km_1 and km_2
        km_1 = node.k_plus_1/(node.k_min_1 + node.k_cat_1)
        km_2 = node.k_plus_2/(node.k_min_2 + node.k_cat_2)

        # Push km values to km_arr
        rates = (km_1, km_2)
        push!(km_arr, rates)

        # Push node names to protein species array
        push!(species_arr, node.name)
    end

    # Loop through links
    for link in links
        if typeof(link) <: TFLink
            km_hat = link.k_hat_plus/link.k_hat_min
            if typeof(link) == Activation
                # Get node indices
                input_idx = indexin([link.input_node], nodes)[1]
                output_idx = indexin([link.output_node], nodes)[1]
                # Update adjacency matrix
                adj_matrix[input_idx, output_idx] = 1
                # Push to tf_links array
                push!(tf_links, link)
                push!(km_hat_arr, km_hat)
                # Increment counter
                num_acti +=1
            elseif typeof(link) == Repression
                # Get node indices
                input_idx = indexin([link.input_node], nodes)[1]
                output_idx = indexin([link.output_node], nodes)[1]
                # Update adjacency matrix
                adj_matrix[input_idx, output_idx] = 1
                adj_matrix[input_idx, output_idx] = 2
                # Push to tf_links array
                push!(tf_links, link)
                push!(km_hat_arr, km_hat)
                # Increment counter
                num_repr +=1
            end
        elseif typeof(link) <: InducerLink
            if typeof(link) == TFActivation
                tf_idx = indexin([link.tf_species], tf_arr)[1]
                node_idx = indexin([link.node], nodes)[1]
                inducer_matrix[tf_idx, node_idx] = 1
                num_tf_acti +=1
                push!(inducer_links, link)
            elseif typeof(link) == TFInactivation
                tf_idx = indexin([link.tf_species], tf_arr)[1]
                node_idx = indexin([link.node], nodes)[1]
                inducer_matrix[tf_idx, node_idx] = 2
                num_tf_inacti +=1
                push!(inducer_links, link)
            end
        elseif typeof(link) == Annhilation
            idx1 = indexin([link.node_1], nodes)[1]
            idx2 = indexin([link.node_2], nodes)[1]
            ann_matrix[idx1, idx2] = 1
            push!(ann_links, link)
            if link.type == "mRNA"
                num_mrna_annihilation +=1
            elseif link.type == "Protein"
                num_protein_annihilation +=1
            end
        end
    end

    # Tally links
    num_ann_links::NTuple{2, Int} = (num_mrna_annihilation, num_protein_annihilation)
    num_tf_links::NTuple{2, Int} = (num_acti, num_repr)
    num_inducer_links::NTuple{2, Int} = (num_tf_acti, num_tf_inacti)
    
    # Create node system
    node_system = NodeSystem(num_nodes, num_tf_links, num_inducer_links, num_ann_links, nodes, tf_links, inducer_links, ann_links, km_arr, km_hat_arr, species_arr, adj_matrix, inducer_matrix, ann_matrix)

    return node_system
end

## Testing
# name = :A
# tx = [1,1,1]
# tl = [1,1,1]
# dg = 0.5 
# init = [1,1,1]

# node_1 = createnode(:A, tx, tl, dg, init)
# link_1 = directactivate(:I1, node_1, [1,1], 1)
# node_2 = createnode(:B, tx, tl, dg, init)
# link_2 = activate(node_1, node_2, [1,1], 2)
# link_3 = tfactivation(:I2, node_2, [1,1], 2, 0.3)

# nodes = compilenodes(node_1, node_2)
# links = compilelinks(link_1, link_2, link_3)

# node_system = createnodesystem(nodes, links)

# println(node_system.num_inducer_links)
# println(sum(node_system.num_inducer_links))

## Rendering
# using Latexify
# render(latexify(vcat(permutedims(node_system.protein_species), node_system.inducer_matrix)))
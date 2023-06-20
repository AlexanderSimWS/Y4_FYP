include("node_system.jl")

# Functions to retrieve reaction data from DSL 
function getnodereactions(lines::Expr)
    @assert lines.head == :block
    name::Symbol = :Z
    tx_rates::Vector{Float64} = []
    tl_rates::Vector{Float64} = []
    d::Real = 0
    init_g::Real = 0
    counter = 1
    for line in lines.args
        if (typeof(line) == Expr)
            if line.args[1] == :~
                init_g = line.args[3].args[2]
                name = line.args[2]
            else
                (rates, rxn) = line.args
                @assert line.head == :tuple
                @assert rates.head == :vect
                @assert rxn.head in (:-->, :call)

                if rxn.head == :call
                    @assert rates.args[1].head == Symbol("=")
                    @assert rxn.args[1] == Symbol("<-->")
                    # fwd_sym = rates.args[1].args[1]
                    k_plus = rates.args[1].args[2]
                    # bwd_sym = rates.args[2].args[1]
                    k_min = rates.args[2].args[2]
                    # return [k_plus, k_min]
                    if (:p in rxn.args[2].args)
                        push!(tx_rates, float(k_plus))
                        push!(tx_rates, float(k_min))
                    elseif (:r in rxn.args[2].args)
                        push!(tl_rates, float(k_plus))
                        push!(tl_rates, float(k_min))
                    end
                elseif (:∅ in rxn.args) 
                    d = rates.args[1].args[2]
                elseif (rxn.head == :-->) && (length(rxn.args[2].args) == 4)
                    @assert rxn.args[2].args[1] == :+
                    @assert rates.args[1].head == Symbol("=")
                    kcat = rates.args[1].args[2]
                    if (counter == 1)
                        push!(tx_rates, float(kcat))
                        counter += 1
                    elseif (counter == 2)
                        push!(tl_rates, float(kcat))
                    end
                end
            end
        end
    end
    return createnode(name, tx_rates, tl_rates, d, init_g)
end

function getnode(sym::Symbol, nodes::Vector{Node})
    for Node in nodes
        if sym == Node.name
            return Node
        end
    end
end

function getedgereaction(lines::Expr)
    @assert lines.head == :block
    type = 0
    input_sym::Symbol = :Z
    output_sym::Symbol = :Y
    rxn_rates::Vector{Float64} = []
    cooperativity::Real = 0
    for line in lines.args
        if (typeof(line) == Expr)
            if (line.head == :->)
                type = 1
                input_sym = line.args[1]
                output_sym = line.args[2].args[2]
            elseif (line.head == :call) && (line.args[1] == :-)
                type = 2
                input_sym = line.args[2]
                output_sym = line.args[3].args[2]
            else
                (rates, rxn) = line.args
                @assert line.head == :tuple
                @assert rates.head == :vect
                push!(rxn_rates, float(rates.args[1].args[2].args[2]))
                push!(rxn_rates, float(rates.args[2].args[2].args[2]))

                @assert rxn.head == :call
                for arg in rxn.args[2].args[3].args
                    if typeof(arg) != Symbol
                        cooperativity = arg
                    end
                end
            end
        end
    end
    return ProtoEdge(type, input_sym, output_sym, rxn_rates[1], rxn_rates[2], cooperativity)
end

function getedges(lines::Expr)
    @assert lines.head == :block
    input1_sym::Symbol = :Z
    input2_sym::Symbol = :K
    output_sym::Symbol = :Y
    rates1::Vector{Float64} = []
    rates2::Vector{Float64} = []
    cooperativity1::Real = 1
    cooperativity2::Real = 1
    counter = 1

    for line in lines.args
        if (typeof(line) == Expr)
            @assert line.head == :tuple
            if line.args[1].head == :->
                input1_sym = line.args[1].args[1]
                output_sym = line.args[1].args[2].args[2]
                input2_sym = lines.args[2].args[2].args[2]
            elseif line.args[1].head == :vect
                (rate1, rate2) = line.args[1].args 
                if counter == 1
                    push!(rates1, rate1.args[2])
                    push!(rates1, rate2.args[2])
                    if typeof(line.args[2].args[2].args[3]) == Expr
                        cooperativity1 = (line.args[2].args[2].args[3].args[2])
                    end
                    counter +=1
                elseif counter == 2
                    push!(rates2, rate1.args[2])
                    push!(rates2, rate2.args[2])
                    if typeof(line.args[2].args[2].args[3]) == Expr
                        cooperativity2 = (line.args[2].args[2].args[3].args[2])
                    end
                    counter +=1
                end
            end
        end
    end
    return ProtoEdge(1, input1_sym, output_sym, rates1[1], rates1[2], cooperativity1), ProtoEdge(2, input2_sym, output_sym, rates2[1], rates2[2], cooperativity2)
end

function createedge(nodes::Vector{Node}, edge::Prototype)
    node_1 = 0
    node_2 = 0
    type = edge.type
    rates::Vector{Float64} = [edge.k_plus, edge.k_min]
    if typeof(edge) == ProtoEdge
        node_1 = getnode(edge.input_node, nodes)
        node_2 = getnode(edge.output_node, nodes)
    elseif typeof(edge) == ProtoAnnEdge
        node_1 = getnode(edge.node_1, nodes)
        node_2 = getnode(edge.node_2, nodes)
    end
    if type == 1
        return activate(node_1, node_2, rates, edge.cooperativity)
    elseif type == 2
        return repress(node_1, node_2, rates, edge.cooperativity)
    elseif typeof(type) == String
        return annhilation(node_1, node_2, type, rates)
    end
end

function createindedge(nodes::Vector{Node}, edge::ProtoIndEdge)
    tf_species = edge.tf_species
    node = getnode(edge.node, nodes)
    k_plus = edge.k_plus
    k_min = edge.k_min
    cooperativity = edge.cooperativiy
    init_I = edge.init_I
    if edge.type == 1
        return tfactivation(tf_species, node, [k_plus, k_min], cooperativity, init_I)
    elseif edge.type == 2
        return tfinactivation(tf_species, node, [k_plus, k_min], cooperativity, init_I)
    end
end

function edge(nodes::Vector{Node}, proto::Prototype)
    if typeof(proto) == ProtoEdge || typeof(proto) == ProtoAnnEdge
        createedge(nodes, proto)
    elseif typeof(proto) == ProtoIndEdge
        createindedge(nodes, proto)
    end
end

function edges(nodes::Vector{Node}, inputedges::Tuple{ProtoEdge, ProtoEdge})
    edge1 = 0
    edge2 = 0
    for i = 1:2
        if inputedges[i].type == 1
            edge1 = createedge(nodes, inputedges[i])
        elseif inputedges[i].type == 2
            edge2 = createedge(nodes, inputedges[i])
        end
    end
    return edge1, edge2
end

function getannedge(lines::Expr)
    @assert lines.head == :block
    node_1 = :Z
    node_2 = :Y
    rates::Vector{Float64} = []
    type::String = "Default"
    for line in lines.args
        if (typeof(line) == Expr)
            if (line.args[1] == :-)
                @assert line.args[2].args[3] == :X
                node_1 = line.args[2].args[2]
                node_2 = line.args[3]
            elseif (line.args[1].head == :vect)
                push!(rates,line.args[1].args[1].args[2])
                push!(rates,line.args[1].args[2].args[2])
                if ((occursin("m", string(line.args[2].args[1].args[2]))))
                    type = "mRNA"
                else
                    type = "Protein"
                end
            end
        end
    end
    return ProtoAnnEdge(node_1, node_2, type, rates[1], rates[2])
end

function getindedge(lines::Expr)
    @assert lines.head == :block
    tf_species = :Z
    node = :Y
    k_plus = 0
    k_min = 0
    cooperativity = 1
    init_I = 0
    type = 0
    for line in lines.args
        if typeof(line) == Expr
            if line.head == :->
                tf_species = line.args[1]
                node = line.args[2].args[2].args[2]
                init_I = line.args[2].args[2].args[3].args[2]
                type = 1
            elseif line.head == :call
                tf_species = line.args[2].args[2]
                node = line.args[2].args[3].args[2]
                init_I = line.args[3].args[2]
                type = 2
            elseif line.args[1].head == :vect
                k_plus = line.args[1].args[1].args[2]
                k_min = line.args[1].args[2].args[2]
                cooperativity = line.args[2].args[2].args[3].args[2]
            end
        end
    end
    return ProtoIndEdge(type, tf_species, node, float(k_plus), float(k_min), cooperativity, init_I)
end

# Macros to wrap functions
macro createnode(rxns)
    getnodereactions(rxns)
end

macro activate(lines::Expr) 
    getedgereaction(lines)
end

macro repress(lines::Expr) 
    getedgereaction(lines)
end

macro combi(lines::Expr)
    getedges(lines)
end

macro annihilate(lines::Expr)
    getannedge(lines)
end

macro TFactivate(lines::Expr)
    getindedge(lines)
end

macro TFinactivate(lines::Expr)
    getindedge(lines)
end

# DSL Showcase (Actual Use)
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

# nodes = compilenodes(node_1, node_2, node_3)

# edge1 = edge(nodes, @repress begin
#     A -! B
#     [k'⁺₂ = 100, k'⁻₂ = 1], g₂ + 3x₁ <--> c₂
# end)

# edge2_3 =  edges(nodes, @combi begin
#     A -> C, B -! C
#     [k⁺ₐ₃ = 100, k⁻ₐ₃ = 1], g₃ + x₁ <--> c₃
#     [k⁺ᵣ₃ = 99, k⁻ᵣ₃ = 2], g₃ + 3x₂ <--> c⁰₃
# end)

# edge4 = edge(nodes, @annhilate begin
#     A -X- B
#     [k⁺₁₂ = 1000, k⁻₁₂ = 0.1], m₁ + m₂ --> m∅
# end)

# edge5 = edge(nodes, @TFinactivate begin
#     I₁ -! A ~ I₁ = 30
#     [κ⁺₁ = 1000, κ⁺₂ = 1], x⁰₁ + 2I₁ <--> x⁺₁
# end)

## DSL Showcase (Extract data)
# @activate nodes (begin
#     A -> B
#     [k'⁺₂ = 100, k'⁻₂ = 1], g₂ + 3x₁ <--> c₂
# end)

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

# edge_5 = @TFinactivate begin
#     I₁ -! A
#     [κ⁺₁ = 1000, κ⁺₂ = 1], x⁰₁ + 2 I₁ <--> x⁺₁
# end
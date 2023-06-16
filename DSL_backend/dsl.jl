include("nodes.jl")

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
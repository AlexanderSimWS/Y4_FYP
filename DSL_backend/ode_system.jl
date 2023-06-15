include("node_system.jl")
using Catalyst, DifferentialEquations, Latexify, Plots, GlobalSensitivity, Statistics

# Result struct
mutable struct SimResult
    sol::ODESolution
    sensitivity::GlobalSensitivity.MorrisResult
end


# Run system
function runsystem(node_system::NodeSystem, p_tot::Real, r_tot::Real, λ::Real, t_span::NTuple{2, Real})
    # Load summarised information from node system
    num_nodes = node_system.num_nodes
    sum_inducer_links = sum(node_system.num_inducer_links)
    num_mrna_annihilation = node_system.num_ann_links[1]
    num_protein_annihilation = node_system.num_ann_links[2]
    nodes = node_system.nodes
    tf_links = node_system.tf_links
    inducer_links = node_system.inducer_links
    ann_links = node_system.ann_links
    km_arr = node_system.km_arr
    km_hat_arr = node_system.km_hat_arr
    adj_matrix = node_system.adj_matrix
    inducer_matrix = node_system.inducer_matrix
    ann_matrix = node_system.ann_matrix

    # Initialise paramters and species
    @parameters t ℼ_tot ρ_tot
    @species m(t)[1:num_nodes] x(t)[1:num_nodes] I(t)[1:sum_inducer_links] x_0(t)[1:sum_inducer_links] m0(t) x0(t)
    # Empty arrary of reactions
    rxns = []
    # Initial conditions
    u0::Vector{Pair{Num, Num}} = []
    # Parameter mapping
    p = [ℼ_tot=>p_tot, ρ_tot=>r_tot]
    t = collect(range(0, stop=5, length=200))

    # Vector to contain species names in order of creation
    ordered_species = []
    # Function to create new symbol that relates final protein species to associated DNA(g), mRNA(m), protein(x)
    createsym(sym::Symbol, substrate::String) = Symbol(string(substrate, "_", sym))

    # Activation function (alpha)
    α(km, x, n) = (km * x^n)/(1 + (km * x^n))

    # Repression function (beta)
    β(km, x, n) = 1/(1 + (km * x^n))

    # Combinatorial function (gamma)
    γ(km_acti, x_plus, n_plus, km_repr, x_min, n_min) = (km_acti * x_plus^n_plus)/(1 + (km_acti * x_plus^n_plus) + (km_repr * x_min^n_min))

    # Search for link idx
    function findlinkidx(node::Node)
        idx_arr = []
        for link in tf_links
            if link.output_node == node
                idx = indexin([link], tf_links)[1]
                push!(idx_arr, idx)
            end
        end
        return idx_arr
    end

    # Find input node idx
    function findinputidx(i::Int)
        input_idx_arr = []
        for j = 1:num_nodes
            if adj_matrix[j,i] != 0
                push!(input_idx_arr, j)
            end
        end
        return input_idx_arr
    end

    # Find inducer link idx
    function findinduceridx(i::Int)
        for k = 1:sum_inducer_links
            if inducer_matrix[k,i] != 0
                return k
            end
        end
    end

    # Find annhilation link idx
    function findannidx(j::Int, i::Int)
        for link in ann_links
            if (link.node_1 == nodes[j]) || (link.node_1 == nodes[i]) && (link.node_2 == nodes[j]) || (link.node_2 == nodes[i])
                return indexin([link], ann_links)[1]
            end
        end
    end

    # Calculate resource demand coefficient
    j1_term = 0
    j2_term = 0

    for i = 1:num_nodes # Scan through column

        km_1 = km_arr[i][1]
        km_2 = km_arr[i][2]
        g_T = nodes[i].init_g
        
        # Categorise nodes
        # 1 Basic node
        if (sum(adj_matrix[:,i]) == 0)
            j1_term += km_1 * g_T
        else
            idx = findlinkidx(nodes[i])[1]
            input_idx = findinputidx(i)[1]
            cooperativity = tf_links[idx].cooperativity
            # 2 Activated Node
            if (sum(adj_matrix[:,i]) == 1)
                j1_term += km_1 * g_T * α(km_hat_arr[idx], x[input_idx], cooperativity)
            # 3 Repressed Node
            elseif (sum(adj_matrix[:,i]) == 2)
                j1_term += km_1 * g_T * β(km_hat_arr[idx], x[input_idx], cooperativity)
            # 4 Combinatorial Node
            elseif (sum(adj_matrix[:,i]) == 3)
                idx2 = findlinkidx(nodes[i])[2]
                input_idx2 = findinputidx(i)[2]
                cooperativity2 = tf_links[idx2].cooperativity
                if typeof(tf_links[idx]) == Activation
                    j1_term += km_1 * g_T * γ(km_hat_arr[idx], x[input_idx], cooperativity, km_hat_arr[idx2], x[input_idx2], cooperativity2)
                elseif typeof(tf_links[idx]) == Repression
                    j1_term += km_1 * g_T * γ(km_hat_arr[idx2], x[input_idx2], cooperativity2, km_hat_arr[idx], x[input_idx], cooperativity)
                end
            end
        end
        j2_term += km_2 * m[i]
    end

    # Calculation of k_eff
    k_eff(km, kcat, j_term, R) = kcat*R*km/(1+j_term)

    # Append reactions
    # 1.1 TF produced is activated
    for i = 1:num_nodes
        kcat_1 = nodes[i].k_cat_1
        kcat_2 = nodes[i].k_cat_2
        k_eff1 = k_eff(km_arr[i][1], kcat_1, j1_term, ℼ_tot)
        k_eff2 = k_eff(km_arr[i][2], kcat_2, j2_term, ρ_tot)
        deg_rate = nodes[i].d_1
        g_T = nodes[i].init_g
        name = nodes[i].name

        push!(rxns, Reaction(deg_rate, [ m[i] ], nothing))
        push!(rxns, Reaction(λ, [ x[i] ], nothing))
        # Basic Node
        if (sum(adj_matrix[:,i]) == 0)
            # Transcription reaction
            push!(rxns, Reaction(k_eff1 * g_T, nothing, [ m[i] ]))
            append!(u0, [m[i]=>0, x[i]=>0])
            push!(ordered_species, createsym(name, "m"))
            # Product not activated/inactivated
            if (sum(inducer_matrix[:,i]) == 0)
                push!(rxns, Reaction(k_eff2, [ m[i] ], [ m[i], x[i] ]))
                push!(ordered_species, createsym(name, "x"))
            else
                inducer_idx = findinduceridx(i)
                k_plus = inducer_links[inducer_idx].k_plus
                k_min = inducer_links[inducer_idx].k_min
                n = inducer_links[inducer_idx].cooperativity
                init_I = inducer_links[inducer_idx].init_I
                append!(u0, [I[inducer_idx]=>init_I, x_0[inducer_idx]=>0])
                push!(rxns, Reaction(λ, [ x_0[inducer_idx] ], nothing))
                # Activated by inducer
                if (sum(inducer_matrix[:,i]) == 1)
                    push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x_0[inducer_idx]]))
                    push!(rxns, Reaction(k_plus, [ x_0[inducer_idx], I[inducer_idx] ], [ x[i] ], [1, n], [1]))
                    push!(rxns, Reaction(k_min, [ x[i] ], [ x_0[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                    push!(ordered_species, createsym(name, "x_0"))
                    push!(ordered_species, createsym(name, "I"))
                    push!(ordered_species, createsym(name, "x"))
                # Inactivated by inducer
                elseif (sum(inducer_matrix[:,i]) == 2)
                    push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x[i]]))
                    push!(rxns, Reaction(k_plus, [ x[inducer_idx], I[inducer_idx] ], [ x_0[i] ], [1, n], [1]))
                    push!(rxns, Reaction(k_min, [ x_0[i] ], [ x[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                    push!(ordered_species, createsym(name, "x"))
                    push!(ordered_species, createsym(name, "I"))
                    push!(ordered_species, createsym(name, "x_0"))
                end
            end
        else
            idx = findlinkidx(nodes[i])[1]
            input_idx = findinputidx(i)[1]
            cooperativity = tf_links[idx].cooperativity
            alpha = α(km_hat_arr[idx], x[input_idx], cooperativity)
            beta = β(km_hat_arr[idx], x[input_idx], cooperativity)
        
            # Activated Node
            if (sum(adj_matrix[:,i]) == 1)
                # Product not activated/inactivated
                push!(rxns, Reaction(k_eff1 * g_T * alpha, nothing, [m[i]]))
                push!(ordered_species, createsym(name, "m"))
                append!(u0, [m[i]=>0, x[i]=>0])
                if (sum(inducer_matrix[:,i]) == 0)
                    push!(rxns, Reaction(k_eff2, [ m[i] ], [ m[i], x[i] ]))
                    push!(ordered_species, createsym(name, "x"))
                else
                    inducer_idx = findinduceridx(i)
                    k_plus = inducer_links[inducer_idx].k_plus
                    k_min = inducer_links[inducer_idx].k_min
                    n = inducer_links[inducer_idx].cooperativity
                    init_I = inducer_links[inducer_idx].init_I
                    append!(u0, [I[inducer_idx]=>init_I, x_0[inducer_idx]=>0])
                    push!(rxns, Reaction(λ, [ x_0[inducer_idx] ], nothing))
                    # Activated by inducer
                    if (sum(inducer_matrix[:,i]) == 1)
                        push!(rxns, Reaction(k_eff2, [m[i]], [x_0[inducer_idx], m[i]]))
                        push!(rxns, Reaction(k_plus, [ x_0[inducer_idx], I[inducer_idx] ], [ x[i] ], [1, n], [1]))
                        push!(rxns, Reaction(k_min, [ x[i] ], [ x_0[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                        push!(ordered_species, createsym(name, "x_0"))
                        push!(ordered_species, createsym(name, "I"))
                        push!(ordered_species, createsym(name, "x"))
                    # Inactivated by inducer
                    elseif (sum(inducer_matrix[:,i]) == 2)
                        push!(rxns, Reaction(k_eff2 * g_T * alpha, [m[i]], [m[i], x[i]]))
                        push!(rxns, Reaction(k_plus, [ x[inducer_idx], I[inducer_idx] ], [ x_0[i] ], [1, n], [1]))
                        push!(rxns, Reaction(k_min, [ x_0[i] ], [ x[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                        push!(ordered_species, createsym(name, "x"))
                        push!(ordered_species, createsym(name, "I"))
                        push!(ordered_species, createsym(name, "x_0"))
                    end
                end        
            # Repressed Node
            elseif (sum(adj_matrix[:,i]) == 2)
                append!(u0, [m[i]=>0, x[i]=>0])
                push!(rxns, Reaction(k_eff1 * g_T * beta, nothing, [m[i]]))
                push!(ordered_species, createsym(name, "m"))
                # Product not activated/inactivated
                if (sum(inducer_matrix[:,i]) == 0)
                    push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x[i]]))
                    push!(ordered_species, createsym(name, "x"))
                else
                    inducer_idx = findinduceridx(i)
                    k_plus = inducer_links[inducer_idx].k_plus
                    k_min = inducer_links[inducer_idx].k_min
                    n = inducer_links[inducer_idx].cooperativity
                    init_I = inducer_links[inducer_idx].init_I
                    append!(u0, [I[inducer_idx]=>init_I, x_0[inducer_idx]=>0])
                    push!(rxns, Reaction(λ, [ x_0[inducer_idx] ], nothing))
                    # Activated by inducer
                    if (sum(inducer_matrix[:,i]) == 1)
                        push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x_0[inducer_idx]]))
                        push!(rxns, Reaction(k_plus, [ x_0[inducer_idx], I[inducer_idx] ], [ x[i] ], [1, n], [1]))
                        push!(rxns, Reaction(k_min, [ x[i] ], [ x_0[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                        push!(ordered_species, createsym(name, "x_0"))
                        push!(ordered_species, createsym(name, "I"))
                        push!(ordered_species, createsym(name, "x"))
                    # Inactivated by inducer
                    elseif (sum(inducer_matrix[:,i]) == 2)
                        push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x[i]]))
                        push!(rxns, Reaction(k_plus, [ x[inducer_idx], I[inducer_idx] ], [ x_0[i] ], [1, n], [1]))
                        push!(rxns, Reaction(k_min, [ x_0[i] ], [ x[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                        push!(ordered_species, createsym(name, "x"))
                        push!(ordered_species, createsym(name, "I"))
                        push!(ordered_species, createsym(name, "x_0"))
                    end
                end
            # Combinatorial node
            elseif (sum(adj_matrix[:,i]) == 3)
                idx2 = findlinkidx(nodes[i])[2]
                input_idx2 = findinputidx(i)[2]
                cooperativity2 = tf_links[idx2].cooperativity
                gamma = 0

                append!(u0, [m[i]=>0, x[i]=>0])
                if typeof(tf_links[idx]) == Activation
                    gamma = γ(km_hat_arr[idx], x[input_idx], cooperativity, km_hat_arr[idx2], x[input_idx2], cooperativity2)
                elseif typeof(tf_links[idx]) == Repression
                    gamma = γ(km_hat_arr[idx2], x[input_idx2], cooperativity2, km_hat_arr[idx], x[input_idx], cooperativity)
                end

                push!(rxns, Reaction(k_eff1 * g_T * gamma, nothing, [m[i]]))
                push!(ordered_species, createsym(name, "m"))
                if (sum(inducer_matrix[:,i]) == 0)
                    push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x[i]]))
                    push!(ordered_species, createsym(name, "x"))
                else
                    inducer_idx = findinduceridx(i)
                    k_plus = inducer_links[inducer_idx].k_plus
                    k_min = inducer_links[inducer_idx].k_min
                    n = inducer_links[inducer_idx].cooperativity
                    init_I = inducer_links[inducer_idx].init_I
                    append!(u0, [I[inducer_idx]=>init_I, x_0[inducer_idx]=>0])
                    push!(rxns, Reaction(λ, [ x_0[inducer_idx] ], nothing))
                    # Activated by inducer
                    if (sum(inducer_matrix[:,i]) == 1)
                        push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x_0[inducer_idx]]))
                        push!(rxns, Reaction(k_plus, [ x_0[inducer_idx], I[inducer_idx] ], [ x[i] ], [1, n], [1]))
                        push!(rxns, Reaction(k_min, [ x[i] ], [ x_0[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                        push!(ordered_species, createsym(name, "x_0"))
                        push!(ordered_species, createsym(name, "I"))
                        push!(ordered_species, createsym(name, "x"))
                    # Inactivated by inducer
                    elseif (sum(inducer_matrix[:,i]) == 2)
                        push!(rxns, Reaction(k_eff2, [m[i]], [m[i], x[i]]))
                        push!(rxns, Reaction(k_plus, [ x[inducer_idx], I[inducer_idx] ], [ x_0[i] ], [1, n], [1]))
                        push!(rxns, Reaction(k_min, [ x_0[i] ], [ x[inducer_idx], I[inducer_idx] ], [1], [1, n]))
                        push!(ordered_species, createsym(name, "x"))
                        push!(ordered_species, createsym(name, "I"))
                        push!(ordered_species, createsym(name, "x_0"))
                    end
                end 
            end        
        end
    end

    for i = 1:num_nodes
        for j = 1:num_nodes
            if ann_matrix[j,i] == 1
                ann_idx = findannidx(j,i)
                k_plus = ann_links[ann_idx].k_plus
                k_min = ann_links[ann_idx].k_min
                type = ann_links[ann_idx].type
                if type == "mRNA"
                    append!(u0, [m0=>0])
                    push!(rxns, Reaction(k_plus, [ m[i], m[j] ], [ m0 ]))
                    push!(rxns, Reaction(k_min, [ m0 ], [ m[i], m[j] ]))
                    push!(ordered_species, :m∅)
                elseif type == "protein"
                    append!(u0, [x0=>0])
                    push!(rxns, Reaction(k_plus, [ x[i], x[j] ],  [x0] ))
                    push!(rxns, Reaction(k_min, [ x0 ], [ x[i], x[j] ]))
                    push!(ordered_species, :x∅)
                end
            end
        end
    end


    @named rn = ReactionSystem(rxns)
    # ode_sys = convert(ODESystem, rn)
    ode_sys = convert(ODESystem, rn, combinatoric_ratelaws=false)
    prob = ODEProblem(ode_sys, u0, t_span, p)
    sol = solve(prob,Rosenbrock23())

    # f1 = function(p)
    #     prob1 = remake(prob;p=p)
    #     sol1 = solve(prob1, Rosenbrock23();saveat=t)
    #     return [maximum(sol1[10,:])]
    # end

    # m = gsa(f1, Morris(total_num_trajectory=1000, num_trajectory=150), [[0.1, 1], [0.1, 1]])
    # println(m.means)
    # println(m.variances)
    # labels = permutedims(string.(ordered_species))
    # plot(sol)
    # render(latexify(ode_sys; mathjax = true))
    # println(rxns)
    oldstd = stdout
    redirect_stdout(devnull)
    # res = SimResult(sol, m);
    redirect_stdout(oldstd)
    return sol
end

## Testing
# name = :A
# tx = [1,1,1]
# tl = [1,1,1]
# dg = 0.5 
# init = 1
# ann = [5,1]

# node_1 = createnode(:A, tx, tl, dg, init)
# link_1 = tfactivation(:I1, node_1, [1,1], 1, 2)

# node_2 = createnode(:B, tx, tl, dg, init)
# link_2 = activate(node_1, node_2, [1,1], 2)
# link_3 = tfactivation(:I2, node_2, [1,1], 2, 0.3)

# node_3 = createnode(:C, tx, tl, dg, init)
# link_4 = repress(node_2, node_3, [0.5,0.5], 4)

# node_4 = createnode(:D, tx, tl, dg, init)
# link_5 = activate(node_1, node_4, [2,2], 1)

# node_5 = createnode(:E, tx, tl, dg, init)
# link_6 = repress(node_5, node_4, [3,2], 1)

# link_7 = annhilation(node_4, node_3, "protein", ann)

# nodes = compilenodes(node_1, node_2, node_3, node_4, node_5)
# links = compilelinks(link_1, link_2, link_3, link_4, link_5, link_6, link_7)

# node_system = createnodesystem(nodes, links)

# # # render(latexify(vcat(permutedims(node_system.protein_species), node_system.adj_matrix)))
# runsystem(node_system, 100, 200, 0.2, (1,1000))

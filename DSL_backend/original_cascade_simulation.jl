include("ode_system.jl")

# Units = uM
tx1 = [100,5,95] # km = 1
tx2 = [17,2,100] # km = 1/6
tx3 = [100,5,95] # km = 1
tl1 = [7,1,104] # km = 1/15
tl2 = [7,1,104] # km = 1/15
tl3 = [7,1,104] # km = 1/15

act1 = [10,1] # k =  10 uM
act2 = [100,1] # k = 100 uM
act3 = [10,1] # k = 10 uM

ind1 = [1,252] # k_hat = 1/252
ind2 = [1,1000] # k_hat = 1/1000

dg = 0.0833 # degradation rate δ

vect = []
conc = [0.01, 0.02, 0.05, 0.08, 0.1, 0.5,1, 5, 10] # AHL concentrations

for i=1:length(conc)
node_1 = createnode(:A, tx1, tl1, dg, 0.06) # Gene concentration = 0.06uM 
link_1 = tfactivation(:I1, node_1, ind1, 1.4, conc[i])

node_2 = createnode(:B, tx2, tl2, dg, 0.06)
link_2 = activate(node_1, node_2, act1, 1) # cooperativity = 1
link_3 = tfactivation(:I2, node_2, ind2, 4, 100) # SAL conc. = 100uM

node_3 = createnode(:C, tx3, tl3, dg, 0.06)
link_4 = activate(node_2, node_3, act2, 1)

compiled_links = compilelinks(link_1, link_2, link_3, link_4)
compiled_nodes = compilenodes(node_1, node_2, node_3)

node_sys = createnodesystem(compiled_nodes, compiled_links)

res = runsystem(node_sys, 1, 1, 0.0067, (0,150)) # ribosome = rnap concentration = 1, λ = 0.0067

push!(vect, maximum(res[10,:]))
end
plot(conc, vect)
## Testing
using DataFrames
using JuMP
using LinearAlgebra
using Clp


function calculate_ptdf(nodes, lines)
    N = collect(1:length(nodes[:, :N]))
    SLACK = findfirst(n -> n.slack, eachrow(nodes))
    incidence = zeros(Int, (length(lines[:, :L]), length(N)))
    for l in 1:length(lines[:, :L])
        incidence[l, [findfirst(n -> n == lines[l, :from], nodes[:, :N]), 
                      findfirst(n -> n == lines[l, :to], nodes[:, :N])]] = [1,-1]
    end
    B = Diagonal(lines[:, :susceptance])
    Bl = B*incidence # Line suceptance Matrix
    Bn = (incidence'*B)*incidence # Nodes suceptance Matrix

    B_inv = zeros(length(N),length(N))
    B_inv[setdiff(N, SLACK), setdiff(N,SLACK)] = inv(Bn[setdiff(N,SLACK), setdiff(N,SLACK)])
    PTDF = Bl*B_inv
    return PTDF
end

plants = DataFrame(
    P = ["pv", "gas", "wind", "coal"],
    mc = [0, 50, 0, 30],
    gmax = [80, 150, 120, 300],
    node = ["N1", "N1", "N2", "N3"]
);

nodes = DataFrame(
    N = ["N1", "N2", "N3"],
    zone = ["Z1", "Z2", "Z2"],
    demand = [150, 30, 200],
    slack = [true, false, false]
);

lines = DataFrame(
    L = ["L1", "L2", "L3"],
    from = ["N2", "N1", "N2"],
    to = ["N1", "N3", "N3"],
    fmax = [40, 100, 100],
    susceptance = [1, 1, 2]
);
zones = DataFrame(
    Z = ["Z1", "Z2"],
);

P = plants[:, :P]
N = nodes[:, :N]
L = lines[:, :L]
Z = zones[:, :Z]

mc = Dict(k => v for (k,v) in eachrow(plants[:, [:P, :mc]]))
gmax = Dict(k => v for (k,v) in  eachrow(plants[:, [:P, :gmax]]))
map_pn = Dict(n => filter(row -> row.node == n, plants)[:, :P] for n in N)
map_zn = Dict(z => filter(row -> row.zone == z, nodes)[:, :N] for z in Z)
demand =  Dict(k => v for (k,v) in eachrow(nodes[:, [:N, :demand]]))

ptdf = calculate_ptdf(nodes, lines)

dispatch = Model(Clp.Optimizer)

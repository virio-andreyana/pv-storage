using Pkg
Pkg.add("Query")

using DataFrames, Query
using JuMP
using LinearAlgebra
using Clp
using CSV # readin of CSV files

df_prod = CSV.read("monthly_output_Wh.csv",DataFrame)
day_month = DataFrame(
    month = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],
    day = [31,28,31,30,31,30,31,31,30,31,30,31]
)




df
for a in df_prod[!,"Apr"]
    println(sum(a for a in df_prod[!,"Apr"]))
end

println(sum(a for a in df_prod[!,"Apr"]))

specs = DataFrame(
    store = ["li-on","pb-acid","redox"], # storage name
    charge = [0.95,0.86,0.83], # charge efficiency
    discharge = [0.95,0.86,0.83], # discharge efficiency
    losses = [0.999,0.96,0.95], # storage loss per hour
    cost = [310,200,500], # cost per kwh
    dod = [100,0.75,100] # depth of discharge
);

S = specs[:, :store]
S = ["li-on"]
H = 1:24

for h in H
    println(h)
end

sum(df_prod[!,"Apr"])/(24*31)


Charge = Dict(k => v for (k,v) in eachrow(specs[:, [:store, :charge]]))
Discharge = Dict(k => v for (k,v) in eachrow(specs[:, [:store, :discharge]]))
Losses = Dict(k => v for (k,v) in eachrow(specs[:, [:store, :losses]]))
Cost = Dict(k => v for (k,v) in eachrow(specs[:, [:store, :cost]]))
Dod = Dict(k => v for (k,v) in eachrow(specs[:, [:store, :dod]]))

sum(charge[s] for s in S)

Consumption = 24

EST = Model(Clp.Optimizer)

@variable(EST,StorageEnergyCapacity[s=S; Discharge[s]>0]>=0)
@variable(EST,StorageLevel[s=S, H; Discharge[s]>0]>=0)
@variable(EST,StorageCharge[s=S, H; Discharge[s]>0]>=0)
@variable(EST,StorageDischarge[s=S, H; Discharge[s]>0]>=0)
@variable(EST,StorageCost[s=S] >= 0)

# Generation must meet demand
@constraint(EST,DemandAdequacy[h in H],
    df_prod[h,"Aug"]/31 == Consumption + sum(StorageCharge[s,h] - StorageDischarge[s,h] for s in S)
)

# storage level depends on previous period's storage level and current period charge/discharge
@constraint(EST,DefineStorageLevel[s in S,h in H ; h>1 && Discharge[s]>0],
    StorageLevel[s,h] == StorageLevel[s,h-1]*Losses[s] + StorageCharge[s,h]*Charge[s] - StorageDischarge[s,h]/Discharge[s]
)

# storage level for first period does not depend on previous level but we set it to 50% energy capacity
@constraint(EST,DefineStorageInitialLevel[s in S,h in H ; h==1 && Discharge[s]>0],
    StorageLevel[s,h] == 0.5*StorageEnergyCapacity[s] + StorageCharge[s,h]*Charge[s] - StorageDischarge[s,h]/Discharge[s]
)

# storage level is limited by storage capacity
@constraint(EST, MaxStorageLevelFunction[s in S, h in H; Discharge[s]>0],
    StorageLevel[s,h] <= StorageEnergyCapacity[s])

@constraint(EST, StorageDailyBalanceFunction[s in S; Discharge[s]>0],
    sum(StorageCharge[s,h] for h in H)*Charge[s] - sum(StorageDischarge[s,h] for h in H) / Discharge[s] == 0)

@constraint(EST, CostFunction[s in S],
    StorageCost[s] == StorageEnergyCapacity[s]*Cost[s]
)

@objective(EST, Min, sum(StorageCost[s] for s in S))
#@objective(EST, Max, Consumption)


optimize!(EST)
# reading our objective value
objective_value(EST)

value.(StorageEnergyCapacity)
value.(StorageLevel)
value.(StorageCharge)
value.(StorageDischarge)
#q1 = @from i in specs begin
#    @where i.store == "li-on"
#    @select i.charge
#    @collect
#end
sum(value.(Curtailment))
sum(df_prod[10,"Aug"])/(31)



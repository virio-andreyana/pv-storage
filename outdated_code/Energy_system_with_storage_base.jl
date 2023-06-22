using JuMP # building models
using DataStructures # using dictionaries with a default value
using Clp # solver for the JuMP model
using CSV # readin of CSV files
using DataFrames # data tables
using Statistics # mean function
using Plots  # generate graphs
using StatsPlots # additional features for plots
using Plots.Measures
include(joinpath(@__DIR__, "colors.jl")) 


# read the csv files
df_demand_timeseries = CSV.read("demand_timeseries.csv", DataFrame,stringtype=String)
df_capacity_factor = CSV.read("capacity_factors_2018.csv", DataFrame,stringtype=String)

### read the new csv files
df_storages = CSV.read("storages.csv", DataFrame,stringtype=String)
df_investmentcoststorage = CSV.read("investmentcoststorage.csv", DataFrame,stringtype=String)
df_e2pratio = CSV.read("e2pratio.csv", DataFrame,stringtype=String)
df_storagechargeefficiency = CSV.read("storagechargeefficiency.csv", DataFrame,stringtype=String)
df_storagedischargeefficiency = CSV.read("storagedischargeefficiency.csv", DataFrame,stringtype=String)
df_storagelosses = CSV.read("storagelosses.csv", DataFrame,stringtype=String)


# readin function for parameters; this makes handling easier
readin(x; default=0,dims=1) = DefaultDict(default,Dict((dims > 1 ? Tuple(row[y] for y in 1:dims) : row[1]) => row[dims+1] for row in eachrow(x)))

# We define our sets from the csv files
hour = 1:120

### add storage here
storages = ["Battery","H2Tank"]#df_storages.storage

# Also, we read our input parameters via csv files
DemandProfile = readin(df_demand_timeseries,default=1/120,dims=1)
CapacityFactor = readin(df_capacity_factor,default=1,dims=1)

### add additional storage parameters
InvestmentCostStorage = readin(df_investmentcoststorage,dims=1)
E2PRatio = readin(df_e2pratio,dims=1)
StorageChargeEfficiency = readin(df_storagechargeefficiency,dims=1)
StorageDisChargeEfficiency = readin(df_storagedischargeefficiency,dims=1)
StorageLosses = readin(df_storagelosses, dims=1)

# our demand
Demand = 30

# define the dictionary for max capacities with specific default value
#MaxCapacity = readin(df_maxcapacity,default=999,dims=1)
#MaxCapacity["SolarPV"]=1000
#MaxCapacity["WindOnshore"]=1000
# instantiate a model with an optimizer

ESM = Model(Clp.Optimizer)

# this creates our variables
@variable(ESM,Production[hour] >= 0)
@variable(ESM,Capacity >=0)

### add variables
@variable(ESM,StorageEnergyCapacity[s=storages; StorageDisChargeEfficiency[s]>0]>=0)
@variable(ESM,StorageCharge[s=storages, hour; StorageDisChargeEfficiency[s]>0]>=0)
@variable(ESM,StorageDischarge[s=storages, hour; StorageDisChargeEfficiency[s]>0]>=0)
@variable(ESM,StorageLevel[s=storages, hour; StorageDisChargeEfficiency[s]>0]>=0)
@variable(ESM,TotalStorageCost[storages] >= 0)

## constraints
# Generation must meet demand
@constraint(ESM, DemandAdequacy[h in hour],
    Production[h] + sum(StorageDischarge[s,h] for s in storages if StorageDisChargeEfficiency[s]>0) == Demand*DemandProfile[h])

# for variable renewables, the production needs to be always at maximum
@constraint(ESM, ProductionFunction_res[h in hour],
    Capacity*CapacityFactor[h]*(1/120) == Production[h])

### Add storage constraints
# storage charge is limited by storage energy capacity and E2PRatio
@constraint(ESM, StorageChargeFunction[s in storages, h in hour; StorageDisChargeEfficiency[s]>0],
    StorageCharge[s,h] <= StorageEnergyCapacity[s]/E2PRatio[s])

# storage discharge is limited by storage energy capacity and E2PRatio
@constraint(ESM, StorageDischargeFunction[s in storages, h in hour; StorageDisChargeEfficiency[s]>0],
    StorageDischarge[s,h] <= StorageEnergyCapacity[s]/E2PRatio[s])

# storage level depends on previous period's storage level and current period charge/discharge
@constraint(ESM, StorageLevelFunction[s in storages, h in hour; h>1 && StorageDisChargeEfficiency[s]>0],
    StorageLevel[s,h] == StorageLevel[s,h-1]*StorageLosses[s] + StorageCharge[s,h]*StorageChargeEfficiency[s] - StorageDischarge[s,h]/StorageDisChargeEfficiency[s])

# storage level for first period does not depend on previous level but we set it to 50% energy capacity
@constraint(ESM, StorageLevelStartFunction[s in storages, h in hour; h==1 && StorageDisChargeEfficiency[s]>0],
    StorageLevel[s,h] == 0.5*StorageEnergyCapacity[s]+ StorageCharge[s,h]*StorageChargeEfficiency[s] - StorageDischarge[s,h]/StorageDisChargeEfficiency[s])

# storage level is limited by storage capacity
@constraint(ESM, MaxStorageLevelFunction[s in storages, h in hour; StorageDisChargeEfficiency[s]>0],
    StorageLevel[s,h] <= StorageEnergyCapacity[s])

# storage cost are the sum of all storage technology costs
@constraint(ESM, StorageCostFunction[s in storages],
    TotalStorageCost[s] == sum(StorageEnergyCapacity[s]*InvestmentCostStorage[s]))

# storage level at the end of a year has to equal storage level at the beginning of year
@constraint(ESM, StorageAnnualBalanceFunction[s in storages; StorageDisChargeEfficiency[s]>0],
    sum(StorageCharge[s,h] for h in hour)*StorageChargeEfficiency[s] - sum(StorageDischarge[s,h] for h in hour) / StorageDisChargeEfficiency[s] == 0)


# the objective function
@objective(ESM, Min, sum(TotalStorageCost[s] for s in storages))

# this starts the optimization
# the assigned solver (here Clp) will takes care of the solution algorithm
optimize!(ESM)
# reading our objective value
objective_value(ESM)

# some result analysis
value.(Production)
value.(Capacity)
value.(StorageEnergyCapacity)
value.(StorageDischarge)
value.(StorageLevel)
value.(StorageCharge)

value.(sum(StorageCharge[s,h] for h in hour)*StorageChargeEfficiency[s] for s in storages)
value.(sum(StorageDischarge[s,h] for h in hour) / StorageDisChargeEfficiency[s] for s in storages)

df_res_production = DataFrame(Containers.rowtable(value,Production; header = [:Hour, :value]))
insertcols!(df_res_production, :Technology => "SolarPV")
#df_res_capacity = DataFrame(Containers.rowtable(value,Capacity; header = [:value]))

df_storage_production = DataFrame(Containers.rowtable(value,StorageDischarge; header = [:Technology, :Hour, :value]))
df_storage_charge = DataFrame(Containers.rowtable(value,StorageCharge; header = [:Technology, :Hour, :value]))
df_storage_level = DataFrame(Containers.rowtable(value,StorageLevel; header = [:Technology, :Hour, :value]))

append!(df_res_production,df_storage_production)

transform!(df_res_production, "Technology" => ByRow(x-> colors[x]) => "Color")
#transform!(df_res_capacity, "Technology" => ByRow(x-> colors[x]) => "Color")

# and some plots
insertcols!(df_res_production, :Fuel => "Power")

groupedbar(
    df_res_production.Fuel,
    df_res_production.value,
    group=df_res_production.Technology,
    bar_position=:stack,
    title="Production by Technology",
    linewidth=0,
    color=df_res_production.Color,
    legend=false
)

if length(storages) > 1
    bar(
        df_res_production.Technology,
        df_res_production.value,
        title="Installed Capacity by Technology",
        color=df_res_production.Color,
        linewidth=0,
        rotation=90
    )
end

gdf_production_by_fuel = groupby(df_res_production, :Fuel)

sto_charge = Dict((row.Technology, row.Hour) => row.value for row in eachrow(df_storage_charge))

n_fuels = length(gdf_production_by_fuel)
plts = map(enumerate(pairs(gdf_production_by_fuel))) do (i,(k,v))
    p = groupedbar(
        v.Hour,
        v.value,
        group=v.Technology,
        bar_position=:stack,
        title="$(k[1])",
        linewidth=0,
        color=v.Color,
        legend=i == n_fuels ? (0.15,-0.5) : false,
        bottom_margin=i == n_fuels ? 20mm : 2mm,
        legend_column=5
    )

    d = [Demand*DemandProfile[h] for h in hour]
    #u = sum(value.(Use)[:,t, k[1]] for t in technologies)
    #c = [sum(get(sto_charge, (s, k[1], h), 0) for s in storages) for h in hour]
    #du = d .+ u.data
    #dus = d .+ c
    plot!(p, hour, d, color=:black, linewidth=2, label="Demand")
    #plot!(p, hour, du, color=:black, linestyle=:dash, linewidth=2, label="Demand + Use")
    #plot!(p, hour, d, color=:Blue, linestyle=:dot, linewidth=2, label="Demand + Storage")

    return p
end

plot(plts..., layout=(n_fuels,1), size=(1200,1200))


sto_prod = Dict((row.Technology, row.Hour) => row.value for row in eachrow(df_storage_production))
sto_lvl = Dict((row.Technology, row.Hour) => row.value for row in eachrow(df_storage_level))

plt_storage_lvl = map(storages) do s
    
   val = [get(sto_prod, (s, h), 0) for h in hour ]

   color = colors[s]
    p = bar(
        hour,
        val,
        title="$(s)",
        label="Production",
        color=:green,
        ylabel="GW",
        linewidth=0
    )

    val_charg = [get(sto_charge, (s, h), 0) for h in hour ]

    bar!(p,
        hour,
        -val_charg,
        title="$(s)",
        label="Charge",
        color=:red,
        linewidth=0,
        legend=:bottomleft
    )


    val_lvl = [get(sto_lvl, (s, h), 0) for h in hour ]


    p2 = twinx(p)

    plot!(
        p2,
        hour,
        val_lvl,
        label = "Level",
        color=:black,
        linewidth=3,
        ylabel="GWh",
        legend=:bottomright
    )

    

    return p
end

plot(plt_storage_lvl...)
### load modules ###
using Distributed, DifferentialEquations
using Statistics: mean
using Dates, TimeZones, Logging, DelimitedFiles

### prepare for file generation ###
# get current time to generate a file name
const starttime_str = Dates.format(now(tz"Asia/Tokyo"), "yymmdd-HHMMSS")
println("start time: ", starttime_str)

# change directory if needed
const wd = "src"
if splitdir(pwd())[end] != wd cd(wd) end
# -> now, we are inside "src" directory!

# make sure that output directory exists
const outdir = "output"
mkpath(outdir)

### logging ###
const io = open(joinpath(outdir, starttime_str * "_log.txt"), "a")
const logger = ConsoleLogger(io)
# set the global logger to logger
global_logger(logger)

#####################################
### model & simulation parameters ###
#####################################
const r = 0.05
const K_vec = [0, 0.1, 0.2]  #[collect(0:0.05:0.25); 10 .^(collect(-0.5:0.5:2))]
const N_vec = [50, 100]  #[50, 100, 200]
const D = 0.005
const ξ = 0.5  # threshold value

const first_seed = 1001
const num_sample = 15

const num_workers = 9
#####################################

# log parameter values
@info "Parameters set" r D ξ num_sample
@info "K_vec:" K_vec
@info "N_vec:" N_vec
flush(io)

# parameter constructor
function init_params(r, K, N::Int, σ, ξ)
    return tuple(r, K, N, σ, ξ, zeros(N))
end

### add workers ###
addprocs(num_workers)
@info "nprocs: $(nprocs())"
flush(io)

### preparation in all worker processes ###
@everywhere begin 
    # load modules in workers
    using DifferentialEquations
    using Statistics: mean

    # index of vector to record escapes
    const id_rec = 6

    ### define model equations ###
    function f_global!(du, u, p, t)
        r, K, N = p  # expand parameters
        du .= -u .* (u .- r) .* (u .- 1) .+ K * (sum(u) / N .- u)
    end

    function g!(du, u, p, t)
        du .= p[4]
    end

    ### prepare options for problems and solvers ###
    # give a particular seed for each trajectory
    const seeds = collect($first_seed:$first_seed + $num_sample)
    prob_func = (prob, i, repeat) -> remake(prob; seed=seeds[i])

    # output reduction: record the mean escape time & retcode
    output_func = (sol, i) -> ((sol.retcode, mean(sol.t[2:end-1])), false)

    ### define methods for escape time measurement ###
    # condition for escape time measurement
    function escape_cond(out, u, t, int)
        # up-crossing corresponds to an escape
        out .= 1.0  # initialise the outcome
        if_escaped = int.p[id_rec] .>= 1
        # define quantity that is positive after an escape
        out[.!if_escaped] .= u[.!if_escaped] .- int.p[5]
    end

    # function to record escape orders in p
    function escape_affect!(int, event_index)
        # record order of escapes in int.p
        int.p[id_rec][event_index] = 1 + maximum(int.p[id_rec])
    end

    ### define methods for termination ###
    # function of terminate condition
    # - if all nodes are above threshold, terminate the simulation
    terminate_cond(u, t, int) = minimum(u) >= int.p[5]

    # callback to terminate the simulation
    terminate_cb = DiscreteCallback(
        terminate_cond, integrator -> terminate!(integrator); 
        save_positions=(true, false)  # only save the data just before termination
    )
end

### core method of performing calculations ###
function sim_given_params!(vec_csvrow, csvpath, K, r, N, D, ξ, trajectories)
    # record K value in output vector
    vec_csvrow[begin] = K
    # ---Define the base problem---
    # update model parameters
    p = init_params(r, K, N, sqrt(2D), ξ)
    # sufficiently long tspan: hard-coded
    tspan = (0.0, 10000.0)
    # prepare initial condition
    u0 = zeros(N)
    # define the base problem
    prob = SDEProblem(f_global!, g!, u0, tspan, p)
    # ---Define the ensemble problem---
    # define the ensemble problem
    ensemble_prob = EnsembleProblem(
        prob, prob_func=prob_func, output_func=output_func
    )
    # ---Solve the ensemble problem with Callbacks---
    @everywhere begin
        # Callback to monitor escapes
        escape_cb = VectorContinuousCallback(
            escape_cond, escape_affect!, nothing, $N;
            # only save the data just before an escape
            save_positions=(true, false)
        )
        # set of Callbacks
        cb_set = CallbackSet(escape_cb, terminate_cb)
    end
    # solve the problem
    sim = solve(ensemble_prob, SOSRA(), EnsembleDistributed(); 
        abstol=5e-4, reltol=1e-3, maxiters=1e11,
        save_everystep=false, callback=cb_set, trajectories=trajectories
    )
    # check if all nodes escaped & record the results
    for (i, out) in enumerate(sim)
        if out[1] == ReturnCode.T(2)  # ReturnCode.Terminated = 2
            vec_csvrow[i+1] = out[2]
        else
            @warn "Calculation was not properly terminated!" K=K N=N sample=i retcode=out[1]
            vec_csvrow[i+1] = -1  # meaningless value
        end
    end
    open(csvpath, "a") do io
        writedlm(io, permutedims(vec_csvrow)::Matrix{T} where T <: Real, ',')
    end
end

### main part: measure escape times ###
# prepare an output csv file
const csvpath_base = joinpath(outdir, starttime_str * "_aet_")

vec_csvrow = Vector{Float64}(undef, 1 + num_sample)
for (id_N, N) in enumerate(N_vec)
    @info "$(Dates.format(now(tz"Asia/Tokyo"), "Y-mm-dd HH:MM:SS.s")) started N=$N"
    flush(io)
    # prepare an output csv file
    csvpath = csvpath_base * "n$(N).csv"
    # write header
    open(csvpath, "w") do io
        writedlm(io, hcat(["K"], permutedims(["$s" for s in seeds])), ',')
    end
    # measure escape times
    for (id_K, K) in enumerate(K_vec)
        @info "$(Dates.format(now(tz"Asia/Tokyo"), "Y-mm-dd HH:MM:SS.s")) - started K=$(round(K, digits=3))"
        flush(io)
        sim_given_params!(vec_csvrow, csvpath, K, r, N, D, ξ, num_sample)
        flush(io)
    end
end

# close the logging file
@info "$(Dates.format(now(tz"Asia/Tokyo"), "Y-mm-dd HH:MM:SS.s")) completed"
close(io)

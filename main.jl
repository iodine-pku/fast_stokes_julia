# Reproduce the results
# 1. Install the required packages
#using Pkg
#Pkg.add("Plots")
#Pkg.add("LinearAlgebra")
#Pkg.add("BenchmarkTools")

using LinearAlgebra

include("multigrid.jl")
include("uzawa.jl")
include("data.jl")
include("inexact_uzawa.jl")



#  Problem 1

println("Reproducing the results for Problem 1 ... ")

nlist = [64, 128, 256, 512, 1024, 2048]
errlist = zeros(length(nlist))
iterlist = zeros(length(nlist))
ct = 1
for n in nlist
    global ct 
    println("Begin computation for n = $n ... ")
    pde_er, u, v, p, iters = vcycle_stokes(n,3,3,1e-8,10,false)
    errlist[ct] = pde_er
    iterlist[ct] = iters
    println("n = $n, iters = $iters, error = $pde_er. Done.")
    ct += 1
end

println("Problem 1 Iteration counts: $iterlist")
println("Problem 1 finished...")



# Problem 2

println("Reproducing the results for Problem 2 ... ")
nlist = [64, 128, 256, 512]
errlist = zeros(length(nlist))
iterlist = zeros(length(nlist))
ct = 1
for n in nlist
    global ct 
    pde_er,iters,x,xx,xxx = uzawa(n)
    errlist[ct] = pde_er
    iterlist[ct] = iters
    println("n = ", n, ", error = ", pde_er)
    ct += 1
end

println("Problem 2 Iteration counts: $iterlist")
println("Problem 2 finished...")

# Problem 3
println("Reproducing the results for Problem 3 ... ")
nlist = [64, 128, 256, 512,1024,2048]
errlist = zeros(length(nlist))
iterlist = zeros(length(nlist))
ct = 1
v1 =2
v2 =2 
tol1 = 1e-8
tol2 = 1e-5
uzawa_iter = 4
pcg_iter = 5
vcycle_iter =2 


for n in nlist
    global ct 
    pde_er,iters,x,xx,xxx = inexact_uzawa(n, v1, v2, tol2, tol1, uzawa_iter, pcg_iter, vcycle_iter,1.0,true)
    errlist[ct] = pde_er
    iterlist[ct] = iters
    println("n = ", n, ", error = ", pde_er)
    ct += 1
end

println(errlist)
println(iterlist)
println("Problem 3 finished...")
println("All results reproduced. Done.")


include("neuralBL/EMSwingup.jl")

moe     = BSON.load("./neuralBL/savedWeights/for_paper_friction12_setdistancetraining_3NN_psi_5-4-4-3-3-1-thetak-5-10-10-4-4-1_elu copy.bson")[:param]
x0      = [0.0f0, 3.1415f0, 0.0f0, -0.1f0]
ψ, θk   = unstackParams(moe)
X,_     = integrate(x0, ψ, θk; totalTimeStep = totalTimeStep)

sleep(1)
animate(X)

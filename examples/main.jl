
include("bouncing_ball.jl")

G_THRESHOLD = 0.001
sys = BouncingBall(Float64)
lcp = Lcp(sys, Float64)

X, t, Λn, Λt, gn = fulltimestep(lcp, Float64)
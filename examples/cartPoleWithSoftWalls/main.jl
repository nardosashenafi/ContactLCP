using ContactLCP

include("dynamics.jl")
# include("../../src/lcp.jl")
# include("../../src/solver.jl")

# window = Window()
# vis = Visualizer()
# open(vis, window)

sys  = CartPoleWithSoftWalls()
lcp  = ContactLCP.Lcp(Float32, sys)
x0   = initialState(0.0f0)
X, t, Λn, Λt = ContactLCP.fulltimestep(lcp, x0, Float32[]; Δt = 0.001f0, 
                            totalTimeStep = 5000, expert=true)
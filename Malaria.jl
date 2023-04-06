using ModelingToolkit
using CairoMakie
using Latexify
using DifferentialEquations


@parameters t ϕ Z ph qv h ψh ψv μh λ β b δ θ
@variables Hs(t) Hp(t) Hr(t) Vs(t) Vp(t)


# Système d'équations différentielles
D = Differential(t)
Malaria_equations = [

    D(Hs) ~ ϕ + δ*Hr - Hs*(λ-Vp-θ),
    D(Hp) ~ λ*Vp*Hs - Hp*(ψh+μh+β),
    D(Hr) ~ β*Hp + θ*Hs + Hr*(ψh+δ),
    D(Vp) ~ λ*Hp*Vs + Z - Vs*ψv,
    D(Vs) ~ Z-λ*Vp*Hs
]

param = [ϕ => 0.05, 
        Z => 1,
        ph => 0.001,
        qv => 0.05,
        h => 0.01,
        ψh => 0.08,
        ψv => 0.08,
        μh => 0.05,
        λ => 0.075,
        β => 0.02,
        b => 0.015,
        δ => 0.025,
        θ => 0.02]
# initial conditions

u0 = [
    Hs => 5000,
    Hp => 1000,
    Hr => 500,
    Vs => 10000,
    Vp => 500]

# Durée
duree = (0.0, 100)

# Solve the system of differential equations
#Latexify.latexify(Malaria_equations) |> render

@named Malaria_system = ODESystem(Malaria_equations)

prob = ODEProblem(Malaria_system, u0, duree, param, jac=true)
sol = solve(prob, saveat=0.0:0.5:100, verbose=true)

#plt = plot(sol.t, Float32.(sol[end]), vars=[Hs, Hp, Hr, Vs, Vp], legend=:right, xlabel="Time (years)", ylabel="Population size")
#display(plt)

# Figure
fig = Figure(; resolution=(1000,500))
timecourse = Axis(fig[1,1]; xlabel="Temps", ylabel="Population")

# Add lines for each variable to the plot
lines!(timecourse, sol[1, :], label="Hs", color=:black)
lines!(timecourse, sol[2, :], label="Hp", color=:red)
lines!(timecourse, sol[3, :], label="Hr", color=:green)
lines!(timecourse, sol[4, :], label="Vs", color=:blue)
lines!(timecourse, sol[5, :], label="Vp", color=:orange)
xlims!(timecourse, (0., 100.))
ylims!(timecourse, (0., 10000.))

current_figure()

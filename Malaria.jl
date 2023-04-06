using ModelingToolkit
using CairoMakie
using Latexify
using DifferentialEquations

# Paramètres et variables

@parameters t ϕ Z ψh ψv μh λ β δ θ
# t : temps
# ϕ : Taux d'entrée pour les humains (naissance et immigration)
# Z : Taux d'entrée pour les vecteurs
# ψ : Taux de sortie des individus (mortalité naturelle et émmigration) 
# μ : Taux de mortalité du au parasite
# λ : Taux de transmission du parasite
# β : Taux de rétablissement des humains
# δ : Taux de perte d'immunité
# θ : Taux de vaccination
@variables Hs(t) Hp(t) Hr(t) Vs(t) Vp(t)
# H : Population d'humains
# V : Population de vecteurs
# s : Individus susceptibles
# p : Individus parasités
# r : Individus immunisés

# Système d'équations différentielles
D = Differential(t)

Malaria_equations = [
          # entrées              #sorties
    D(Hs) ~ ϕ +  δ*Hr         - Hs*(λ*Vp+θ+ψh),
    D(Hp) ~      λ*Vp*Hs      - Hp*(ψh+μh+β),
    D(Hr) ~      β*Hp + θ*Hs  - Hr*(ψh+δ),
    D(Vs) ~ Z                 - λ*Hp*Vs - Vs*ψv,
    D(Vp) ~      λ*Hp*Vs      - Vp*ψv

]

param = [ϕ => 150,    
        Z => 0.1,     
        ψh => 0.001,  
        ψv => 0.001,  
        μh => 0.001,  
        λ => 0.0001,  
        β => 0.02,    
        δ => 0.05,    
        θ => 0.001]   


# Conditions initiales

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

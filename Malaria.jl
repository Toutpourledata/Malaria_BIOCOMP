# Chargement des packages
using ModelingToolkit
using CairoMakie
using Latexify
using DifferentialEquations

# Paramètres et variables

@parameters t ϕ Z ph qv h ψh ψv μh λ β b δ θ
@variables Hs(t) Hp(t) Hr(t) Vs(t) Vp(t)


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
#valeurs des paramètres

param = [
        #paramètres NON multiplicatifs
        ϕ => 2,  #Entrants dans Hs (naissance et immigration)
        Z => 600,      #Entrants (vecteur)

        #paramètres multiplicatifs
        ψh => 0.001,  #Taux de sortie des humains
        ψv => 0.01,  #Taux de sortie des vecteurs
        μh => 0.001,  #Taux de mortalité dû au parasite
        λ => 0.0001,  #Taux de transmission
        β => 0.8,   #Taux de guérison
        δ => 0.5,  #Taux de perte d'immunité
        θ => 0.01]   #Taux de vaccination

# conditions initiales
u0 = [
    Hs => 5000,   # humains susceptibles (black)
    Hp => 15,   # humains parasités (red)
    Hr => 3,    # humains rétablis (green)
    Vs => 100000,  # vecteur susceptibles (orange)
    Vp => 100]    # vecteur parasités (blue)
# durée
duree = (0.0, 100)

# résoudre le système d'équations différentielles
@named Malaria_system = ODESystem(Malaria_equations)

prob = ODEProblem(Malaria_system, u0, duree, param, jac=true)
sol = solve(prob, saveat=0.0:0.5:100, verbose=true)

# graphique
fig = Figure(; resolution=(1000, 1000), layout=(2, 1))

# humains
timecourse_human = Axis(fig[1, 1]; xlabel="Temps", ylabel="Population humaine")
lines!(timecourse_human, sol[1, :]+sol[2, :]+sol[3, :], label="total_H", color=:yellow)
lines!(timecourse_human, sol[1, :], label="Hs", color=:black)
lines!(timecourse_human, sol[2, :], label="Hp", color=:red)
lines!(timecourse_human, sol[3, :], label="Hr", color=:green)
xlims!(timecourse_human, (0., 100.))
ylims!(timecourse_human, (0., 10000.))

# vecteurs 
timecourse_vector = Axis(fig[2, 1]; xlabel="Temps", ylabel="Population vecteur")
lines!(timecourse_vector, sol[4, :]+sol[5, :], label="total_V", color=:purple)
lines!(timecourse_vector, sol[4, :], label="Vs", color=:blue)
lines!(timecourse_vector, sol[5, :], label="Vp", color=:orange)
xlims!(timecourse_vector, (0., 100.))
ylims!(timecourse_vector, (0., 120000.))

# Display
fig

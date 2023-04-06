# Chargement des packages
using ModelingToolkit
using CairoMakie
using Latexify
using DifferentialEquations

# Ajout des paramètres et variables
@parameters t ϕ Z ph qv h ψh ψv μh λ β b δ θ
#
@variables Hs(t) Hp(t) Hr(t) Vs(t) Vp(t)
# H : Humains
# V : Vecteurs
# s : Individus susceptibles
# p : Individus parasités
# r : Individus immunisés

# système d'équations différentielles 
D = Differential(t)
Malaria_equations = [

    D(Hs) ~ ϕ + δ*Hr - Hs*(λ-Vp-θ),
    D(Hp) ~ λ*Vp*Hs - Hp*(ψh+μh+β),
    D(Hr) ~ β*Hp + θ*Hs + Hr*(ψh+δ),
    D(Vp) ~ λ*Hp*Vs + Z - Vs*ψv,
    D(Vs) ~ Z-λ*Vp*Hs
]
#valeurs des paramètres
param = [ϕ => 0.05,  #Entrants dans Hs (naissance et immigration)
        Z => 1,      #Entrants (vecteur)
        ph => 0.001, #Probabilité qu'un immigrant soit parasité?
        qv => 0.05,  #
        h => 0.01,   #
        ψh => 0.08,  #Taux de sortie des humains
        ψv => 0.08,  #Taux de sortie des vecteurs
        μh => 0.05,  #Taux de mortalité dû au parasite
        λ => 0.075,  #Taux de transmission
        β => 0.02,   #Taux de guérison
        b => 0.015,  #
        δ => 0.025,  #Taux de perte d'immunité
        θ => 0.02]   #Taux de vaccination

# conditions initiales
u0 = [
    Hs => 5000,   # humains susceptibles
    Hp => 1000,   # humains parasités
    Hr => 500,    # humains rétablis
    Vs => 10000,  # vecteur susceptibles
    Vp => 500]    # vecteur parasités

# durée
duree = (0.0, 100)

# résoudre le système d'équations différentielles
@named Malaria_system = ODESystem(Malaria_equations)

prob = ODEProblem(Malaria_system, u0, duree, param, jac=true)
sol = solve(prob, saveat=0.0:0.5:100, verbose=true)

# graphique
fig = Figure(; resolution=(1000,500))
timecourse = Axis(fig[1,1]; xlabel="Temps", ylabel="Population")

# on ajoute une ligne pour chaque variable
lines!(timecourse, sol[1, :], label="Hs", color=:black)
lines!(timecourse, sol[2, :], label="Hp", color=:red)
lines!(timecourse, sol[3, :], label="Hr", color=:green)
lines!(timecourse, sol[4, :], label="Vs", color=:blue)
lines!(timecourse, sol[5, :], label="Vp", color=:orange)
xlims!(timecourse, (0., 100.))
ylims!(timecourse, (0., 10000.))

#générer le graphique
current_figure()

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
    D(Vs) ~ Z *(1+(cos(2*pi*t/365)))      - λ*Hp*Vs - Vs*ψv,
    D(Vp) ~      λ*Hp*Vs      - Vp*ψv

]
#valeurs des paramètres

param = [ϕ => 150,  #Entrants dans Hs (naissance et immigration)
        Z => 0.1,      #Entrants (vecteur)
        ψh => 0.01,  #Taux de sortie des humains
        ψv => 0.01,  #Taux de sortie des vecteurs
        μh => 0.01,  #Taux de mortalité dû au parasite
        λ => 0.001,  #Taux de transmission
        β => 0.02,   #Taux de guérison
        δ => 0.05,  #Taux de perte d'immunité
        θ => 0.1
        ]#Taux de vaccination  

# conditions initiales
u0 = [
    Hs => 5000,   # humains susceptibles (black)
    Hp => 1000,   # humains parasités (red)
    Hr => 500,    # humains rétablis (green)
    Vs => 5000,  # vecteur susceptibles (orange)
    Vp => 50]    # vecteur parasités (blue)
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

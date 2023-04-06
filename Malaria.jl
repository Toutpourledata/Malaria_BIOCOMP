# Chargement des packages
using ModelingToolkit
using CairoMakie
using DifferentialEquations

@parameters t Zh Zv ψh ψv μ λ β δ θ
# t : Temps (jours)
# Z : Taux d'entrée dans la population (naissance et immigration)
# ψ : Taux de sortie de la population (mortalité naturelle et émigration)
# μ : Taux de mortalité des humains du au parasite
# λ : Taux de transmission
# β : Taux de rétablissement procurant une immunité temporaire
# δ : Taux de perte de l'immunité
# θ : Taux de vaccination
@variables Hs(t) Hp(t) Hr(t) Vs(t) Vp(t)
# H : Population d'humains
# V : Population de vecteurs
# s : Individus susceptibles
# p : Individus parasités
# r : Individus rétablis

# Système d'équations différentielles
D = Differential(t)

Malaria_equations = [
          # entrées              #sorties
    D(Hs) ~ Zh +  δ*Hr         - Hs*(λ*Vp+θ+ψh),
    D(Hp) ~      λ*Vp*Hs      - Hp*(ψh+μh+β),
    D(Hr) ~      β*Hp + θ*Hs  - Hr*(ψh+δ),
    D(Vs) ~ Zv                 - λ*Hp*Vs - Vs*ψv,
    D(Vp) ~      λ*Hp*Vs      - Vp*ψv

]
#valeurs des paramètres

param = [Zh => 150,  #Entrants dans Hs (naissance et immigration)
        Zv => 0.1,      #Entrants (vecteur)
        ψh => 0.001,  #Taux de sortie des humains
        ψv => 0.001,  #Taux de sortie des vecteurs
        μh => 0.001,  #Taux de mortalité dû au parasite
        λ => 0.0001,  #Taux de transmission
        β => 0.02,   #Taux de guérison
        δ => 0.05,  #Taux de perte d'immunité
        θ => 0.001]   #Taux de vaccination

# conditions initiales
u0 = [
    Hs => 5000,   # humains susceptibles (black)
    Hp => 1000,   # humains parasités (red)
    Hr => 500,    # humains rétablis (green)
    Vs => 6040,  # vecteur susceptibles (orange)
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

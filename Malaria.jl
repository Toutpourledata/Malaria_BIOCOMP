# Chargement des packages
using ModelingToolkit
using CairoMakie
using DifferentialEquations

@parameters t Zh Zv ψh ψv μ λ β δ θ
# t : Temps (jours)
# Z : Taux d'entrée dans la population (naissance et immigration)
# ψ : Taux de sortie de la population (mortalité naturelle et émigration)
# μ : Taux de mortalité des humains du au parasite
# λ : Taux de transmission du parasite
# β : Taux de rétablissement procurant une immunité temporaire
# δ : Taux de perte de l'immunité
# θ : Taux de vaccination

# Valeurs des paramètre
param = [Zh => 0.2,  
        Zv => 50,   
        ψh => 0.001, 
        ψv => 0.01, 
        μ => 0.01, 
        λ => 0.0001, 
        β => 0.2,   
        δ => 0.2,  
        θ => 0.1] 

@variables Hs(t) Hp(t) Hr(t) Vs(t) Vp(t)
# H : Population d'humains
# V : Population de vecteurs
# s : Individus susceptibles
# p : Individus parasités
# r : Individus rétablis

# Conditions initiales
u0 = [
    Hs => 5500,  
    Hp => 5,  
    Hr => 1000,   
    Vs => 6500,  
    Vp => 50]  

#Création du modèle
D = Differential(t)
Malaria_equations = [
          # entrées              #sorties
    D(Hs) ~ Zh +  δ*Hr         - Hs*(λ*Vp+θ+ψh),
    D(Hp) ~      λ*Vp*Hs      - Hp*(ψh+μ+β),
    D(Hr) ~      β*Hp + θ*Hs  - Hr*(ψh+δ),
    D(Vs) ~ Zv                 - λ*Hp*Vs - Vs*ψv,
    D(Vp) ~      λ*Hp*Vs      - Vp*ψv
]  

# Durée de la simulation
duree = (0.0, 350)

# Résolution du modèle
@named Malaria_system = ODESystem(Malaria_equations)
prob = ODEProblem(Malaria_system, u0, duree, param, jac=true)
sol = solve(prob, saveat=0.0:0.5:350, verbose=true)

# Création d'un canva pour le graphique
fig = Figure(; resolution=(1000,500))

# Identification des axes
timecourse = Axis(fig[1,1]; title = "Modèle de la propagation de la Malaria", xlabel="Temps", ylabel="Individus dans la population")
xlims!(timecourse, (0., 350.))
ylims!(timecourse, (0., 10000.))

# Ajout des courbes de tendances
lines!(timecourse, sol[1, :], label="Hs", color=:black)
lines!(timecourse, sol[2, :], label="Hp", color=:red)
lines!(timecourse, sol[3, :], label="Hr", color=:green)
lines!(timecourse, sol[4, :], label="Vs", color=:blue)
lines!(timecourse, sol[5, :], label="Vp", color=:orange)

# Ajout d'une légende
axislegend()

# Visualisation du graphique
current_figure()

#
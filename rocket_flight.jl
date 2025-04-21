
using DifferentialEquations
using Roots
using Plots

#============#
# Parameters #
#============#

thrust = 5900.0  # Newtons
burn_rate = 3.0  # L/s
fuel_density = 0.98  # kg/L
drag = 2.0
rocket_mass = 400.0  # kg
fuel_volume = 100.0  # L
g = 9.8

#=========#
# Phase 1 #
#=========#

m0 = rocket_mass + fuel_volume*fuel_density
Δm = burn_rate*fuel_density

m(t) =  m0 - Δm*t

f1(v,p,t) = -g + (thrust - drag*v)/m(t)

v1_init = 0.0

t_fuel_spent = fuel_volume/burn_rate
tspan1 = (0.0, t_fuel_spent) 

vprob1 = ODEProblem(f1, v1_init, tspan1)
vsol1 = solve(vprob1)

h1(y,p,t) = vsol1(t)
y1_init = 0.0
yprob1 = ODEProblem(h1, y1_init, tspan1)
ysol1 = solve(yprob1)

#=========#
# Phase 2 #
#=========#

f2(v,p,t) = -g - drag*v/rocket_mass

max_velocity = vsol1(t_fuel_spent)
@show max_velocity

v2_init = vsol1(t_fuel_spent)

t_final = 70.0
tspan2 = (t_fuel_spent, t_final)

vprob2 = ODEProblem(f2, v2_init, tspan2)
vsol2 = solve(vprob2)

h2(y,p,t) = vsol2(t)
y2_init = ysol1(t_fuel_spent)
yprob2 = ODEProblem(h2, y2_init, tspan2)
ysol2 = solve(yprob2)

#========================#
# Compute maximum height #
#========================#

t_max_height = find_zero(vsol2, tspan2)
max_height = ysol2(t_max_height)
@show t_max_height
@show max_height

#==================================#
# Compute impact time and velocity #
#==================================#

t_impact = find_zero(ysol2, tspan2)
v_impact = vsol2(t_impact)
@show t_impact
@show v_impact

#=======#
# Plots #
#=======#

v1(t) = thrust/drag*(1 - (m(t)/m0)^(drag/Δm)) +
          g*m0/(Δm - drag)*(m(t)/m0 - (m(t)/m0)^(drag/Δm))

v2(t) = exp(drag/rocket_mass*(t_fuel_spent - t))*
          (max_velocity + g*rocket_mass/drag) - g*rocket_mass/drag

p1 = plot(legend=:none)
ylabel!("velocity (m/s)")
hline!([0], c=:black)
plot!(vsol1, c=1)
plot!(vsol2, c=1)
#plot!(v1, 0.0, t_fuel_spent, c=:black)
#plot!(v2, t_fuel_spent, t_final, c=:black)
scatter!([0.0], [v1_init], c=1)
scatter!([t_fuel_spent], [max_velocity], c=1)
scatter!([t_max_height], [0.0], c=1)
scatter!([t_impact], [v_impact], c=1)
ylims!(v_impact*1.1, max_velocity*1.1)
xlabel!("")

k = drag
mr = rocket_mass
t̄ = t_fuel_spent

y1(t) = -thrust/drag*m(t)/Δm*(1 - Δm/(k+Δm)*(m(t)/m0)^(k/Δm)) +
          g*m0/(k-Δm)*m(t)/Δm*(m(t)/2m0 - Δm/(k+Δm)*(m(t)/m0)^(k/Δm)) +
          m0*(thrust - 1/2*g*m0)/Δm/(k + Δm)

y2(t) = mr/k*(v2_init + g*mr/k)*(1 - exp(k/mr*(t̄ - t))) +
          g*mr/k*(t̄ - t) + y2_init

p2 = plot(legend=:none)
ylabel!("height (m)")
hline!([0], c=:black)
plot!(ysol1, c=2)
plot!(ysol2, c=2)
#plot!(y1, 0.0, t_fuel_spent, c=:black)
#plot!(y2, t_fuel_spent, t_final, c=:black)
scatter!([0.0], [y1_init], c=2)
scatter!([t_fuel_spent], [y2_init], c=2)
scatter!([t_max_height], [max_height], c=2)
scatter!([t_impact], [0.0], c=2)
ylims!(-200.0, max_height*1.1)
xlabel!("time (s)")

plt = plot(p1, p2, layout=(2,1), size=(800,600))
xlims!(0.0, t_impact*1.01)

savefig(plt, "rocket_flight.png")

display(plt)
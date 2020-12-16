using DrWatson
@quickactivate "PECUZAL_Julia"

using DelayEmbeddings
using DynamicalSystemsBase
using Test
using DelimitedFiles
using Random
using PyPlot
pygui(true)

## We here focus on the dependence of the L-statistic on its parameter Tw and
# exemplary look at the Lorenz63 and the Lorenz96 system, as well as on random
# uniformly distributed noise.

## Generate Data:
# Lorenz63

# The time series have been computed in the script `/test/timeseries/produce_timeseries.jl`
s = readdlm("./test/timeseries/lorenz_pecora_uni_x.csv")
s1 = vec(s[1:5000]) # input timeseries = x component of lorenz
w1 = estimate_delay(vec(s), "mi_min")    # 1st minimum of mutual inf

# Lorenz96
idx = 487
N = 8 # number of oscillators
Fs = 3.5:0.002:5 # parameter spectrum
F = Fs[idx]
dt = 0.1 # sampling time
total = 5000  # time series length

# Parameters analysis:
ε = 0.05  # recurrence threshold
dmax = 10   # maximum dimension for traditional tde
lmin = 2   # minimum line length for RQA
trials = 80 # trials for MCDTS
taus = 0:100 # possible delays

#pick one time series
t_idx = 2

# init Lorenz96
u0 = [0.590; 0.766; 0.566; 0.460; 0.794; 0.854; 0.200; 0.298]
lo96 = Systems.lorenz96(N, u0; F = F)

#set_parameter!(lo96, 1, F)
data = trajectory(lo96, total*dt; dt = dt, Ttr = 2500 * dt)
s2 = data[1:5000,t_idx]
w2 = estimate_delay(vec(s2), "mi_min")    # 1st minimum of mutual inf

# Random numbers
Random.seed!(1234)
s3 = rand(1000)
w3 = estimate_delay(vec(s3), "mi_min")    # 1st minimum of mutual inf

# normalize data
D1 = regularize(Dataset(s1))
D2 = regularize(Dataset(s2))
D3 = regularize(Dataset(s3))

# generate embedding candidate
Y_t1 = DelayEmbeddings.hcat_lagged_values(D1, vec(Matrix(D1)), w1)
Y_t2 = DelayEmbeddings.hcat_lagged_values(D2, vec(Matrix(D2)), w2)
Y_t3 = DelayEmbeddings.hcat_lagged_values(D3, vec(Matrix(D3)), w3)

Tw = 100    # maximum Tw
L11 = zeros(Tw)
L12 = zeros(Tw)
L21 = zeros(Tw)
L22 = zeros(Tw)
L31 = zeros(Tw)
L32 = zeros(Tw)
for T = 1:Tw
    L11[T] = uzal_cost(D1; w = w1, samplesize = 1, Tw = T)
    L12[T] = uzal_cost(Y_t1; w = w1, samplesize = 1, Tw = T)
    L21[T] = uzal_cost(D2; w = w2, samplesize = 1, Tw = T)
    L22[T] = uzal_cost(Y_t2; w = w2, samplesize = 1, Tw = T)
    L31[T] = uzal_cost(D3; w = w3, samplesize = 1, Tw = T)
    L32[T] = uzal_cost(Y_t3; w = w3, samplesize = 1, Tw = T)
end

# Plotting

panelnames = ["A" "B" "C"]

fig = figure(figsize=[10,5])

lwg = 3         # linewidth of the graph
lwa = 2         # linewidth of the axis
fsp = 18        # Fontsize of panelnames
fsa = 12        # Fontsize of the axis
fsl = 10        # Fontsize of the legendentries
fst = 20        # Fontsize of title
axislabelsize = 12 # axislabelsize
ticklabelsize = 10  # labelsize of ticks
ms = 30         # markersize of maxima
ms_style = "*"  # markerstyle

subplots_adjust(
#left = 0.07,
#right = 0.96,
#bottom = 0.05,
#top = 0.96,
wspace = 0.35,
#hspace = 0.55)
)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
        "font.size" => fsa,
        "font.weight" => "bold",
        "axes.labelweight" => "bold",
        "axes.titleweight" => "bold",
        "axes.labelsize" => axislabelsize,
        "axes.linewidth" => lwa,
        "xtick.labelsize" => ticklabelsize,
        "ytick.labelsize" => ticklabelsize,
        "legend.fontsize" => fsl)
        #"text.usetex" => true)
merge!(rcParams, font0)

ax = fig.add_subplot(1, 3, 1)
p1 = plot(1:Tw, L11, label="time series", linewidth=lwg)
color1 = p1[1].get_color()
scatter(1:Tw, L11, label="")
p2 = plot(1:Tw, L12, label="2-dim embedding", linewidth=lwg)
color2 = p2[1].get_color()
scatter(1:Tw, L12, label="")
legend()
ax.text(-0.2, 1.1, panelnames[1], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("Lorenz 63")
xlabel(L"T_M")
ylabel("L")
grid()

ax = fig.add_subplot(1, 3, 2)
p1 = plot(1:Tw, L21, label="time series", linewidth=lwg)
color1 = p1[1].get_color()
scatter(1:Tw, L21, label="")
p2 = plot(1:Tw, L22, label="2-dim embedding", linewidth=lwg)
color2 = p2[1].get_color()
scatter(1:Tw, L22, label="")
#legend()
ax.text(-0.2, 1.1, panelnames[2], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("Lorenz 96")
xlabel(L"T_M")
grid()

ax = fig.add_subplot(1, 3, 3)
p1 = plot(1:Tw, L31, label="time series", linewidth=lwg)
color1 = p1[1].get_color()
scatter(1:Tw, L31, label="")
p2 = plot(1:Tw, L32, label="2-dim embedding", linewidth=lwg)
color2 = p2[1].get_color()
scatter(1:Tw, L32, label="")
#legend()
ax.text(-0.2, 1.1, panelnames[3], transform=ax.transAxes,
     fontsize=fsp, fontweight="bold", va="top")
title("Random numbers")
xlabel(L"T_M")
grid()

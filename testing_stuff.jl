using PyPlot
pygui(true)

x = tr[1:5000,1]
y = tr[1:5000,2]
z = tr[1:5000,3]

Y1 = regularize(DelayEmbeddings.hcat_lagged_values(Dataset(x), x, 18))
Y2 = regularize(DelayEmbeddings.hcat_lagged_values(Dataset(x), y, 0))
Y3 = regularize(DelayEmbeddings.hcat_lagged_values(Dataset(x), z, 0))
Y4 = regularize(DelayEmbeddings.hcat_lagged_values(Dataset(y), z, 0))

Tw_max = 100
Lx = zeros(Tw_max)
Ly = zeros(Tw_max)
Lz = zeros(Tw_max)
LY1 = zeros(Tw_max)
LY2 = zeros(Tw_max)
LY3 = zeros(Tw_max)
LY4 = zeros(Tw_max)

for Tw = 1:Tw_max
    Lx[Tw] = uzal_cost(Dataset(x); Tw = Tw, w = 17, samplesize=1)
    Ly[Tw] = uzal_cost(Dataset(y); Tw = Tw, w = 17, samplesize=1)
    Lz[Tw] = uzal_cost(Dataset(z); Tw = Tw, w = 17, samplesize=1)

    LY1[Tw] = uzal_cost(Y1; Tw = Tw, w = 17, samplesize=1)
    LY2[Tw] = uzal_cost(Y2; Tw = Tw, w = 17, samplesize=1)
    LY3[Tw] = uzal_cost(Y3; Tw = Tw, w = 17, samplesize=1)
    LY4[Tw] = uzal_cost(Y4; Tw = Tw, w = 17, samplesize=1)
end

Tws = 17:Tw_max
figure(figsize=(20,10))
plot(Tws, Lx[Tws], label = "Lx", linewidth = 2.5)
plot(Tws, Ly[Tws], label = "Ly", linewidth = 2.5)
plot(Tws, Lz[Tws], label = "Lz", linewidth = 2.5)
plot(Tws, LY1[Tws], label = "LXX18", linestyle="dotted", linewidth = 2.5)
plot(Tws, LY2[Tws], label = "LXY", linestyle="solid", linewidth = 2.5)
plot(Tws, LY3[Tws], label = "LXZ", linestyle="dashed", linewidth = 2.5)
plot(Tws, LY4[Tws], label = "LYZ", linestyle="dashdot", linewidth = 2.5)
legend()
grid()
xlabel("Tw")
ylabel("L")

Tws = 1:Tw_max
figure(figsize=(20,10))
subplot(2,1,1)
plot(Tws, Lx[Tws], label = "Lx", linewidth = 2.5)
plot(Tws, LY1[Tws], label = "LXX18", linestyle="dotted", linewidth = 2.5)
plot(Tws, LY2[Tws], label = "LXY", linestyle="solid", linewidth = 2.5)
plot(Tws, LY3[Tws], label = "LXZ", linestyle="dashed", linewidth = 2.5)
legend()
grid()
xlabel("Tw")
ylabel("L")
subplot(2,1,2)
plot(Tws, LY1[Tws].-Lx[Tws], label = "LXX18-DIFF", linewidth = 2.5)
plot(Tws, LY2[Tws].-Lx[Tws], label = "LXY-DIFF", linestyle="dotted", linewidth = 2.5)
plot(Tws, LY3[Tws].-Lx[Tws], label = "LXZ-DIFF", linestyle="solid", linewidth = 2.5)
legend()
grid()
xlabel("Tw")
ylabel("L-DIFF")


Tws = 1:Tw_max
figure(figsize=(20,10))
subplot(2,1,1)
plot(Tws, Ly[Tws], label = "Ly", linewidth = 2.5)
plot(Tws, LY4[Tws], label = "LYZ", linestyle="dotted", linewidth = 2.5)
legend()
grid()
xlabel("Tw")
ylabel("L")
subplot(2,1,2)
plot(Tws, LY4[Tws].-Ly[Tws], label = "LYZ-DIFF", linewidth = 2.5)
legend()
grid()
xlabel("Tw")
ylabel("L-DIFF")


L_sing = zeros(Tw_max)
L_mult = zeros(Tw_max)
for Tw = 1:Tw_max
    L_mult[Tw] = uzal_cost(regularize(Y); Tw = Tw, w = w, samplesize=1)
    L_sing[Tw] = uzal_cost(regularize(Y_s); Tw = Tw, w = w, samplesize=1)
end

figure()
plot(L_sing, label="single")
plot(L_mult, label="multi")
legend()
grid()

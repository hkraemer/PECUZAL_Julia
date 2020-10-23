
PyPlot.figure()
PyPlot.plot(data1[1:2:1000], marker=".", linestyle="none")

Y = embed(data1[1:20000], 3, 60)

figure()
plot3D(Y[:,1], Y[:,2], Y[:,3])

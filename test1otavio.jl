using PyCall

plt = pyimport("matplotlib.pyplot")

x = 3 .* 0:2pi/50:2pi
y = 2

ll = 32863.97979
strin = pyeval("f'{$ll:.5e}'")

plt.plot(sin.(x), label="$y")

plt.legend()

plt.show()

import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

fig, ax = plt.subplots()

GPa = 1e9 # Define Giga Pascal
K = 16.0*GPa # Set bulk modulus
G = 7.0*GPa # Set shear modulus
nu = (3*K-2*G) / (2*(3*K+G)) # Calculate Poissons ratio
rho = 2150.0 # Define density in kg/m^3
M = K + 4/3*G # Calculate compressional modulus
Vp = math.sqrt(M/rho) # Calculate P-wave velocity
Vs = math.sqrt(G/rho) # Calculate S-wave velocity
f = 30.0 # Set frequency
T = 1 / f # Calculate the period
omega = 2 * math.pi * f # Calculate the angular frequency
lmbd = Vp * T # Calculate the wave-length
k = 2 * math.pi / lmbd # Calculate the wave-number

# First figure

x = 1000
T = np.arange(0, 0.1, 0.001)
u = np.exp((k * x - omega * T) * 1j)
# plt.plot(T, u.real)

# plt.show()

# second figure

X = np.arange(-300, 300, 1)
t = 0
u = np.exp((k * X - omega * t) * 1j)
# plt.plot(X, u.real)

# plt.show()

# animation!!!

line, = ax.plot(X, u)

def animate(t):
    line.set_ydata(np.exp((k * X - omega * t) * 1j))
    return line,

def init():
    line.set_ydata(np.ma.array(X, mask=True))
    return line,

ani = animation.FuncAnimation(
    fig, animate, np.arange(1, 100, 0.001),
    init_func=init,
    interval=25,
    blit=True)

plt.title('title!!!')
plt.xlabel('distance (m)')
plt.ylabel('amplitude')
plt.show()

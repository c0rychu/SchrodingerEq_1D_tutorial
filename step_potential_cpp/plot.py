#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.rc('savefig', dpi=300)


data = np.genfromtxt('./step_potential.txt', delimiter=',', dtype=complex)
x = data[0, 1:].astype(float)
V = data[1, 1:].astype(float)
sol_t = data[2:, 0].astype(float)
sol_y = data[2:, 1:]

dx = x[1]-x[0]

# Make a plot of psi0 and V
fig = plt.figure(figsize=(15, 5))
plt.plot(x, V*0.01, "k--", label=r"$V(x) (x0.01)")
plt.plot(x, np.abs(sol_y[0, :])**2, "r", label=r"$\vert\psi(t=0,x)\vert^2$")
plt.plot(x, np.real(sol_y[0, :]), "g", label=r"$Re\{\psi(t=0,x)\}$")
plt.legend(loc=1, fontsize=8, fancybox=False)
fig.savefig('step_initial@2x.png')

print("Total Probability: ", np.sum(np.abs(sol_y[0, :])**2)*dx)


# Plotting
fig = plt.figure(figsize=(6, 4))
for i, t in enumerate(sol_t):
    plt.plot(x, np.abs(sol_y[i, :])**2)                  # Plot Wavefunctions
    print("Total Prob. in frame", i, "=", np.sum(np.abs(sol_y[i, :])**2)*dx)   # Print Total Probability (Should = 1)
plt.plot(x, V * 0.001, "k--", label=r"$V(x) (x0.001)")   # Plot Potential
plt.legend(loc=1, fontsize=8, fancybox=False)
fig.savefig('step@2x.png')


# Animation
fig = plt.figure(figsize=(8, 6))

ax1 = plt.subplot(2, 1, 1)
ax1.set_xlim(0, 10)
ax1.set_ylim(-1, 3)
title = ax1.set_title('')
line11, = ax1.plot([], [], "k--", label=r"$V(x)$ (x0.001)")
line12, = ax1.plot([], [], "b", label=r"$\vert \psi \vert^2$")
plt.legend(loc=1, fontsize=8, fancybox=False)

ax2 = plt.subplot(2, 1, 2)
ax2.set_xlim(0, 10)
ax2.set_ylim(-2, 2)
line21, = ax2.plot([], [], "k--", label=r"$V(x)$ (x0.001)")
line22, = ax2.plot([], [], "r", label=r"$Re\{ \psi \}$")
plt.legend(loc=1, fontsize=8, fancybox=False)


def init():
    line11.set_data(x, V * 0.001)
    line21.set_data(x, V * 0.001)
    return line11, line21


def animate(i):
    line12.set_data(x, np.abs(sol_y[i, :])**2)
    line22.set_data(x, np.real(sol_y[i, :]))
    title.set_text('Time = {0:1.3f}'.format(sol_t[i]))
    return line12, line22


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(sol_t), interval=200, blit=True)


# Save the animation into a short video
print("Generating mp4")
anim.save('step.mp4', fps=15, extra_args=['-vcodec', 'libx264'], dpi=150)
print("Generating GIF")
# anim.save('step@2x.gif', writer='pillow', fps=15)
anim.save('step@2x.gif', writer='imagemagick', fps=15, dpi=150)

import os, sys
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['axes.linewidth'] = 1.3
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)

# ---- Inputs ----------------------------
type_plot  = 'log'
program    = "run"

N          = int(sys.argv[1])

ti         = 1e-6
tf         = 1.
meas       = 1000

hthr       = 1.0
sthr       = 0.75

gamma_p    = 0.2
gamma_m    = 0.1

# ----------------------------------------


gnet       = gamma_p-gamma_m
if type_plot == 'log':
    gnet = -gnet

integrator = "Euler"
outf       = "output.txt"
# -----------------------


main_dir = os.getcwd()
exe = main_dir + "/" + program

options = ["--N", str(N), "--meas", str(meas), "--thr", str(hthr), "--sthr", str(sthr),
           "--integrator", str(integrator), "--ti", str(ti), "--tf", str(tf), "--file", outf,
           "--gamma_p", str(gamma_p), "--gamma_m", str(gamma_m)]
command = [exe] + options

process = subprocess.Popen(command)

plt.ion()
fig, ax = plt.subplots(figsize=(10,7))

line, = ax.plot([], [], 'b-', lw=3, label='Numerical')

trange = np.geomspace(ti,tf,100)
ax.plot(trange, -gnet*trange*N**2/4, c='k', ls='--', lw=1, label='Perturbative')

ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$\langle J_z\rangle$")
ax.set_xlim(ti, tf)
if type_plot == 'log':
    ax.set_ylim(1e-1, N)
else:
    ax.set_ylim(-N/1.8, 0)
ax.legend(frameon=False, fontsize=15)
ax.set_title(r'$N=%d$'%N, fontsize=18)
ax.set_xscale('log')
if type_plot == 'log':
    ax.set_yscale('log')



output_file = os.path.join(os.getcwd(), outf)

if os.path.exists(output_file):
    os.remove(output_file)

xdata, ydata = [], []

while process.poll() is None:
    if os.path.exists(output_file):
        try:
            with open(output_file, 'r') as f:
                lines = f.readlines()

            new_x, new_y = [], []
            for data_line in lines:
                if data_line.strip():
                    parts = data_line.strip().split()
                    if len(parts) >= 2:
                        try:
                            x, y = float(parts[0]), float(parts[1])
                            if type_plot == 'log':
                                y = np.abs(y)
                            new_x.append(x)
                            new_y.append(y)
                        except ValueError:
                            continue

            if new_x and new_y:
                xdata, ydata = new_x, new_y
                line.set_xdata(xdata)
                line.set_ydata(ydata)
                ax.relim()
                ax.autoscale_view()
                plt.draw()
                plt.pause(0.1)
        except Exception as e:
            print("Error while reading/plotting:", e)

    time.sleep(0.1)

plt.ioff()
plt.show()

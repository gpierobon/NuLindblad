import os, sys
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['axes.linewidth'] = 1.3
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)

# ---- Inputs -----------
program    = "run"
N          = int(sys.argv[1])
ti         = 0.001
tf         = 2.0
out        = 100
thr        = 1.0
integrator = "Euler"
# -----------------------

gnet = 0.01

main_dir = os.getcwd()
#main_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
exe = main_dir + "/" + program

options = ["--N", str(N), "--out", str(out), "--thr", str(thr),
           "--integrator", str(integrator), "--tf", str(tf)]
command = [exe] + options

process = subprocess.Popen(command)

plt.ion()
fig, ax = plt.subplots(figsize=(10,7))

line, = ax.plot([], [], 'b-', lw=3, label='Numerical')

trange = np.linspace(0,tf,100)
ax.plot(trange, -gnet*trange*N**2/4, c='k', ls='--', lw=1, label='Perturbative')

ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$\langle J_z\rangle$")
ax.set_xlim(ti, tf)
ax.set_ylim(-N/1.8, 0)
ax.legend(frameon=False, fontsize=15)
ax.set_title(r'$N=%d$'%N, fontsize=18)
#ax.set_xscale('log')



output_file = os.path.join(os.getcwd(), "output.txt")

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
                            new_x.append(x)
                            new_y.append(y)
                        except ValueError:
                            continue  # Skip malformed lines

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

# Final plot update after process ends
plt.ioff()
plt.show()

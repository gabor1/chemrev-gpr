import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


fig2e_data = np.load('fig2e_data.npz')
gp_stds = fig2e_data['gp_stds']
true_errors = fig2e_data['true_errors']
vac_dists = fig2e_data['vac_dists']

# Set colors.
init_val = 176/255
fin_val = 100/255

cdict = {'red':   [(0.0,  0.0, init_val),
                   (1.0,  0, 1.0)],

         'green': [(0.0,  0.0, 0.0),
                   (1.0, 0.0, 1.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (1.0,  init_val, 1.0)]}

colors = mpl.colors.LinearSegmentedColormap('grey-red', cdict)

plt.scatter(gp_stds, true_errors, c=vac_dists, marker='.', cmap=colors)

plt.plot([0, 0.2], [0, 0.4], 'k:', label='2 $\sigma$')
plt.xlabel('GP uncertainty (eV/$\AA$)')
plt.ylabel('True error (eV/$\AA$)')
cbar = plt.colorbar()
cbar.set_label('Distance to vacancy ($\AA$)')

plt.savefig("fig2e.pdf")

plt.show()

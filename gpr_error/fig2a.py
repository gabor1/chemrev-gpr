import numpy as np
import matplotlib.pyplot as plt

fig2a_data = np.load('fig2a_data.npz')

two_cuts = fig2a_data['two_cuts']
two_noise = fig2a_data['two_noise']
two_rmse = fig2a_data['two_rmse']

three_cuts = fig2a_data['three_cuts']
three_noise = fig2a_data['three_noise']
three_rmse = fig2a_data['three_rmse']

plt.plot(two_cuts, two_noise, 'g-', label='2b GP')
plt.plot(two_cuts, two_rmse, 'g:', label='2b true')
plt.plot(three_cuts, three_noise, 'b-', label='3b GP')
plt.plot(three_cuts, three_rmse, 'b:', label='3b true')

plt.xlabel('Cutoff ($\AA$)')
plt.ylabel('Error (meV/$\AA$)')
plt.legend()


plt.savefig("fig2a.pdf")
plt.show()
###
### This is a Python2.7 script, runs in the QUIP docker tagged :py2
###

from quippy import Descriptor
from ase import Atoms
import qlab as qp
import numpy as np
l_max = 12
n_max = 12
costs = np.zeros((l_max+1, n_max))
print costs
for (i, j), cost in np.ndenumerate(costs):
    string = "soap l_max=" + str(i) + " n_max=" + str(j+1) + " cutoff=4.2 atom_sigma=0.5 zeta=4.0 cutoff_transition_width=0.5"
    #print string
    single_atom = qp.Atoms(Atoms("C", positions=[[0, 0, 0]], cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]]))
    single_atom.cutoff = 4.2
    single_atom.calc_connect()
    desc = Descriptor(string)
    length = len(desc.calc(single_atom)["descriptor"][0])
    costs[i, j] = length-1

    i += 1
    j += 1
    
costs = costs.astype("int")
print costs

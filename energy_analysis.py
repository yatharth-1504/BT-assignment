import math as ma
import matplotlib.pyplot as plt
import numpy as np
from residues import residues_at_ph7
from residues import residues_at_ph2
from residues import residues_at_ph5
from residues import residues_at_ph14


# function to calculate total_energy and distance
# k_new (in a dilectric) = k/(dielectic_contant)
# dist_in_ang respresents distance in argstrong, it needs to be converted to S.I units
# assumption particles at distance greater than 15 angstroms are considered non-interactive
def cal_energy_and_dist(residues, pot_energy):
    k = 8.98755*10**9
    dielectric_const = 78
    k_new = k/dielectric_const
    total_energy = 0
    for i in range(len(residues) - 1):
        for j in range(i + 1, len(residues)):
            dist_in_ang = ma.sqrt((residues[i]["x"] - residues[j]["x"])**2 + (
                residues[i]["y"] - residues[j]["y"])**2 + (residues[i]["z"] - residues[j]["z"])**2)
            if dist_in_ang < 15:
                dist = dist_in_ang*10**-10
                pot_energy[i][j] = (pot_energy[j][i]) = k_new * \
                    residues[i]['charge'] * \
                    residues[j]['charge']/dist
                total_energy += pot_energy[i][j]
    return total_energy


print("Total no. of charged-residues_at_ph7 :", len(residues_at_ph7))
# pot_energy_at_ph7 represents matrix of pair wise interaction energy, neglecting the pot_enrgy between the same residues_at_ph7
pot_energy_at_ph7 = np.zeros((len(residues_at_ph7), len(residues_at_ph7)))
total_energy_at_ph7 = cal_energy_and_dist(
    residues_at_ph7, pot_energy_at_ph7)
print("TOTAL POTENTIAL ENERGY(at ph 7) : ", total_energy_at_ph7, "Joules")
print("Magnitude of Pairwise INTERACTION ENERGIES(at ph 7) :\n",
      np.abs(pot_energy_at_ph7))


print("Total no. of charged-residues_at_ph2 :", len(residues_at_ph2))
# pot_energy_at_ph2 represents matrix of pair wise interaction energy, neglecting the pot_enrgy between the same residues_at_ph2
pot_energy_at_ph2 = np.zeros((len(residues_at_ph2), len(residues_at_ph2)))
total_energy_at_ph2 = cal_energy_and_dist(
    residues_at_ph2, pot_energy_at_ph2)
print("TOTAL POTENTIAL ENERGY(at ph 2) : ", total_energy_at_ph2, "Joules")
print("Magnitude of Pairwise INTERACTION ENERGIES(at ph 2) :\n",
      np.abs(pot_energy_at_ph2))


print("Total no. of charged-residues_at_ph5 :", len(residues_at_ph5))
# pot_energy_at_ph5 represents matrix of pair wise interaction energy, neglecting the pot_enrgy between the same residues_at_ph5
pot_energy_at_ph5 = np.zeros((len(residues_at_ph5), len(residues_at_ph5)))
total_energy_at_ph5 = cal_energy_and_dist(
    residues_at_ph5, pot_energy_at_ph5)
print("TOTAL POTENTIAL ENERGY(at ph 5) : ", total_energy_at_ph5, "Joules")
print("Magnitude of Pairwise INTERACTION ENERGIES(at ph 2) :\n",
      np.abs(pot_energy_at_ph5))


print("Total no. of charged-residues_at_ph14 :", len(residues_at_ph14))
# pot_energy_at_ph14 represents matrix of pair wise interaction energy, neglecting the pot_enrgy between the same residues_at_ph14
pot_energy_at_ph14 = np.zeros((len(residues_at_ph14), len(residues_at_ph14)))
total_energy_at_ph14 = cal_energy_and_dist(
    residues_at_ph14, pot_energy_at_ph14)
print("TOTAL POTENTIAL ENERGY(at ph 14) : ", total_energy_at_ph14, "Joules")
print("Magnitude of Pairwise INTERACTION ENERGIES(at ph 14) :\n",
      np.abs(pot_energy_at_ph14))


# Ploting and Interaction Energy
energies = [total_energy_at_ph2, total_energy_at_ph5,
            total_energy_at_ph7, total_energy_at_ph14]
ph = [2, 5, 7, 14]
plt.plot(ph, energies)
plt.scatter(ph, energies)
plt.title("ph vs interaction energy graph")
plt.xlabel('----------------ph -------------->')
plt.ylabel('total interaction energies')
plt.colorbar()
plt.show()

# in the graph above it can be observed that the interaction energy obtained at ph 5 is minimum.
# hence, at ph 5 our proteien is most stable
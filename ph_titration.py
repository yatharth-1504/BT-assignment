# required import
from re import A
import numpy as np
import matplotlib.pyplot as plt

kw = 10**-14
y = []
conc = np.linspace(0.001, 1, 999)

# here we are using cubic polynomial for ph


def cal_ph_values(k):
    list = np.zeros(999)
    for i in range(len(conc)):
        coeff = [1, (k + conc[i]/(1 + conc[i])), conc[i]*k /
                 (1 + conc[i]) - k/(1 + conc[i]) + kw, k*kw]
        np.roots(coeff)
        coeff = sorted(x.real for x in coeff)
        list[i] = {coeff[0] > 0: coeff[0], coeff[1]
                   > 0: coeff[1]}.get(True, coeff[2])
    return list


# i have considered ka values of glutamic acid
k_values = {"ka1": 10**-4.25, "ka2": 10**-9.47}


plt.plot(conc, -np.log(cal_ph_values(k_values["ka1"])),
         label=f"curve for glu ka1", color="blue")
plt.plot(conc, -np.log(cal_ph_values(k_values["ka2"])),
         label=f"curve for glu ka2", color="red")
plt.xlabel("conc.")
plt.ylabel("ph")
plt.show()

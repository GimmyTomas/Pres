import numpy as np
import matplotlib.pyplot as plt

q = 0.001
def Pgw(R, eccentricity):
    return ((32/5) * q**2 / R**5) * (1+73*eccentricity**2/24+37*eccentricity**4/96) / (1-eccentricity**2)**(7/2)

alphas = []
Ns = []
R_list = []
Pres_list = []

inclination = 145 # in degrees, should be an integer as it's used as an index in ion211

with open("Pres-data/100_varying_alpha_inclin_"+str(inclination)+".txt", "r") as f:
    for line in f:
        parts = line.strip().split()
        if not parts:
            continue  # skip empty lines
        alpha = float(parts[0])
        N = int(parts[1])
        A = list(map(float, parts[2 : 2 + N]))
        B = list(map(float, parts[2 + N : 2 + 2 * N]))

        alphas.append(alpha)
        Ns.append(N)
        R_list.append(A)
        Pres_list.append(B)

ion211=np.loadtxt('211_ion_inclin.txt', unpack=True)

R_break = np.zeros(len(alphas))

# fig,ax=plt.subplots()

file = open("results/inclination_"+str(inclination)+".txt", "w")
for i in range(len(alphas)):
    Rion = ion211[0] * (alphas[i]/0.2)**(-2)
    PionPgw = ion211[2+4*inclination] * (alphas[i]/0.2)**(-5)

    Rres = np.array(R_list[i])
    PresPgw = np.array(Pres_list[i]) / Pgw(Rres,0)

    PionPgw_interp = np.interp(Rres, Rion, PionPgw)

    PtotPgw = PresPgw + PionPgw_interp
    j_min = np.argmin(PtotPgw)

    R_break[i] =  Rres[j_min]
    Mc_break =  0.01 * (- 1 / PtotPgw[j_min])

    file.write(str(alphas[i])+" "+str(R_break[i])+" "+str(Mc_break)+"\n")
    # print(alphas[i], R_break[i], Mc_break)

    # ax.plot(Rres, PresPgw + PionPgw_interp)

file.close()
# plt.show()
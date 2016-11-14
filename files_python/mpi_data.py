from pylab import *
from numpy import zeros


readfile = loadtxt('mpi_lattice40.txt')
readfile2 = loadtxt('mpi_lattice60.txt')
readfile3 = loadtxt('mpi_lattice100.txt')
readfile4 = loadtxt('mpi_lattice140.txt')

temp = readfile[:,1]

# divide each array with the number of spins
mean_energy = readfile[:,2]/(40*40)
mean_mag = readfile[:,3]/(40*40)
spec_heat = readfile[:,4]/(40*40)
mag_susc = readfile[:,5]/(40*40)

mean_energy2 = readfile2[:,2]/(60*60)
mean_mag2 = readfile2[:,3]/(60*60)
spec_heat2 = readfile2[:,4]/(60*60)
mag_susc2 = readfile2[:,5]/(60*60)

mean_energy3 = readfile3[:,2]/(100*100)
mean_mag3 = readfile3[:,3]/(100*100)
spec_heat3 = readfile3[:,4]/(100*100)
mag_susc3 = readfile3[:,5]/(100*100)

mean_energy4 = readfile4[:,2]/(140*140)
mean_mag4 = readfile4[:,3]/(140*140)
spec_heat4 = readfile4[:,4]/(140*140)
mag_susc4 = readfile4[:,5]/(140*140)


f, axarr = subplots(2,2)
suptitle(r'Mean energy, mean magnetic moment, specific heat and susceptibility against temperature', fontsize=18)

axarr[0, 0].plot(temp, mean_energy, label=r'$40 \times 40$')
axarr[0, 0].plot(temp, mean_energy2, '-.', label=r'$60 \times 60$')
axarr[0, 0].plot(temp, mean_energy3, '--', label=r'$100 \times 100$')
axarr[0, 0].plot(temp, mean_energy4, '-.', label=r'$140 \times 140$')
axarr[0, 0].set_xlabel(r'Temperature [kT/J]', fontsize=15)
axarr[0, 0].set_ylabel(r'Mean energy $<E>/N$', fontsize=15)
axarr[0, 0].legend(loc='best')

axarr[0, 1].plot(temp, mean_mag, label=r'$40 \times 40$')
axarr[0, 1].plot(temp, mean_mag2, '-.', label=r'$60 \times 60$')
axarr[0, 1].plot(temp, mean_mag3, '--', label=r'$100 \times 100$')
axarr[0, 1].plot(temp, mean_mag4, '-.', label=r'$140 \times 140$')
axarr[0, 1].set_xlabel(r'Temperature [kT/J]', fontsize=15)
axarr[0, 1].set_ylabel(r'Mean magnetic moment $<|M|>/N$', fontsize=15)
axarr[0, 1].legend(loc='best')

axarr[1, 0].plot(temp, spec_heat, label=r'$40 \times 40$')
axarr[1, 0].plot(temp, spec_heat2, '-.', label=r'$60 \times 60$')
axarr[1, 0].plot(temp, spec_heat3, '--', label=r'$100 \times 100$')
axarr[1, 0].plot(temp, spec_heat4, '-.', label=r'$140 \times 140$')
axarr[1, 0].set_xlabel(r'Temperature [kT/J]', fontsize=15)
axarr[1, 0].set_ylabel(r'Specific heat $C_V/N$', fontsize=15)
axarr[1, 0].legend(loc='best')

axarr[1, 1].plot(temp, mag_susc, label=r'$40 \times 40$')
axarr[1, 1].plot(temp, mag_susc2, '-.', label=r'$60 \times 60$')
axarr[1, 1].plot(temp, mag_susc3, '--', label=r'$100 \times 100$')
axarr[1, 1].plot(temp, mag_susc4, '-.', label=r'$140 \times 140$')
axarr[1, 1].set_xlabel(r'Temperature [kT/J]', fontsize=15)
axarr[1, 1].set_ylabel(r'Magnetic susceptibility $\chi /N$', fontsize=15)
axarr[1, 1].legend(loc='best')

show()	


# find maximum of peaks -> that's where the phase transition is

Tc_Cv = [ temp[where(spec_heat == max(spec_heat))], temp[where(spec_heat2 == max(spec_heat2))], temp[where(spec_heat3 == max(spec_heat3))], temp[where(spec_heat4 == max(spec_heat4))] ]
Tc_Chi = [ temp[where(mag_susc == max(mag_susc))], temp[where(mag_susc2 == max(mag_susc2))], temp[where(mag_susc3 == max(mag_susc3))], temp[where(mag_susc4 == max(mag_susc4))] ]

max1 = sum(Tc_Cv)/len(Tc_Cv)
max2 = sum(Tc_Chi)/len(Tc_Chi)

print Tc_Cv
print Tc_Chi

Tc_numerical = ( max1 + max2 ) / 2

print 'Numerical critical temperature: T_c = %g' % (Tc_numerical)

'''
Numerical critical temperature: T_c = 2.29375
'''
# Metropolis algorithm not very efficient close to the critical temperature, 
# as seen by the susceptibility values










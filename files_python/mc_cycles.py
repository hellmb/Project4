
from pylab import *
from numpy import zeros, histogram
from scipy import stats


TwentyByTwenty = False
if TwentyByTwenty:
    L = 20

    read_file = loadtxt('mc_cycles20.txt')

    E = read_file[:,0]
    M = read_file[:,1]
    mc_cycles = read_file[:,2]

    f, axarr = subplots(2)

    axarr[0].plot(mc_cycles, E, label=r'T = 1.0 kT/J')
    #axarr[0].plot(mc_cycles, E, label=r'T = 2.4 kT/J')
    axarr[0].set_title(r'Variation of energy for a $20 \times 20$ lattice', fontsize=18)
    #axarr[0].set_title(r'Variation of energy for an ordered $20 \times 20$ lattice', fontsize=18)
    axarr[0].set_ylabel(r'Energy', fontsize=15)
    axarr[0].legend()

    axarr[1].plot(mc_cycles, M, label=r'T = 1.0 kT/J')
    #axarr[1].plot(mc_cycles, M, label=r'T = 2.4 kT/J')
    axarr[1].set_title(r'Variation in magnetization for a $20 \times 20$ lattice', fontsize=18)
    #axarr[1].set_title(r'Variation in magnetization for an ordered $20 \times 20$ lattice', fontsize=18)
    axarr[1].set_xlabel(r'Monte Carlo cycles', fontsize=15)
    axarr[1].set_ylabel(r'Magnetization', fontsize=15)
    axarr[1].legend()

    axarr[0].set_ylim([-900, -600])         # temperature = 1
    #axarr[1].set_ylim([-500, -300])        # random matrix, temperature = 1
    #axarr[1].set_ylim([350, 450])           # ordered matrix, temperature = 1

    show()

accepted_configs = False
if accepted_configs:
    read_file = loadtxt('accepted1.0.txt')
    read_file2 = loadtxt('accepted2.4.txt')
    read_file3 = loadtxt('accepted5.0.txt')
    read_file4 = loadtxt('accepted10.0.txt')
    acc_conf1 = read_file[:,0]
    acc_conf2 = read_file2[:,0]
    acc_conf3 = read_file3[:,0]
    acc_conf4 = read_file4[:,0]
    mc_cycles = read_file2[:,1]

    # find slope of each line
    slope1, intercept, a, b, c_ = stats.linregress(mc_cycles, acc_conf1)
    slope2, intercept, a, b, c_ = stats.linregress(mc_cycles, acc_conf2)
    slope3, intercept, a, b, c_ = stats.linregress(mc_cycles, acc_conf3)
    slope4, intercept, a, b, c_ = stats.linregress(mc_cycles, acc_conf4)

    plot(mc_cycles, acc_conf1, label=r'T = 1.0 [kT/J] $\rightarrow$ slope = %g' % slope1)
    plot(mc_cycles, acc_conf2, label=r'T = 2.4 [kT/J] $\rightarrow$ slope = %g' % slope2)
    plot(mc_cycles, acc_conf3, label=r'T = 5.0 [kT/J]$\rightarrow$ slope = %g' % slope3)
    plot(mc_cycles, acc_conf4, label=r'T = 10.0 [kT/J]$\rightarrow$ slope = %g' % slope4)
    title(r'Number of accepted configurations', fontsize=20)
    xlabel(r'Monte Carlo cycles', fontsize=15)
    ylabel(r'Accepted configurations', fontsize=15)
    legend(loc='best')
    xlim([0, 2000])
    ylim([-0.5E5, 1.E6])
    show()

probability_dist = False
if probability_dist:
    read_file1 = loadtxt('probability1.0.txt')
    read_file2 = loadtxt('probability2.4.txt')

    energy_temp1 = read_file1[:,0]
    energy_temp2 = read_file2[:,0]

    mc_cycles = read_file1[:,1]
    last_mc = mc_cycles[-1] - 1000

    # find variance for each array of energies
    sum1 = sum(energy_temp1)
    sum2 = sum(energy_temp2)
    mean_energy1 = sum1/last_mc
    mean_energy2 = sum2/last_mc

    mean_sqr1 = sum(energy_temp1*energy_temp1) / last_mc
    mean_sqr2 = sum(energy_temp2*energy_temp2) / last_mc

    var1 = mean_sqr1 - (mean_energy1 * mean_energy1)
    var2 = mean_sqr2 - (mean_energy2 * mean_energy2)

    # ensure that the sum is 1
    sum_norm = [float(i)/sum(energy_temp1) for i in energy_temp1]
    sum_norm2 = [float(i)/sum(energy_temp2) for i in energy_temp2]
    print sum(sum_norm), sum(sum_norm2)

    print mean_energy1, mean_energy2


    # make subplots
    f, axarr = subplots(2)
    
    axarr[0].hist(energy_temp1, normed=True, label=r'$\sigma_E^2$ = %g' % var1)
    axarr[0].set_title(r'Probability distribution of energies at T = 1.0 [kT/J]', fontsize=15)
    axarr[0].set_ylabel(r'Probability', fontsize=15)
    axarr[0].legend()
    
    axarr[1].hist(energy_temp2, 100, normed=True, label=r'$\sigma_E^2$ = %g' % var2)
    axarr[1].set_title(r'Probability distribution of energies at T = 2.4 [kT/J]', fontsize=15)
    axarr[1].set_xlabel(r'Energy', fontsize=15)
    axarr[1].set_ylabel(r'Probability', fontsize=15)
    axarr[1].legend()
    show()   



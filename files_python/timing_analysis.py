from pylab import *

# values sampled from openmpi_project4.pro -> main.cpp
processors = [1, 2, 3, 4, 8, 12, 16]
time = [233.76, 124.254, 100.704, 81.3027, 80.9062, 80.6782, 79.1862]

plot(processors, time)
title(r'Timing analysis of a $40 \times 40$ lattice', fontsize=18)
xlabel(r'Number of processors', fontsize=15)
ylabel(r'Time [s]', fontsize=15)
show()
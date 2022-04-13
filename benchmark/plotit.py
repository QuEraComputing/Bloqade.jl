# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:12:59 2022

@author: jwurtz
"""

from matplotlib.pyplot import *
from numpy import *

rcParams.update({'font.size': 22})
close('all')

fil = open('pulser_compare.dat','r')
dat = [line.split(',') for line in fil.readlines()]
fil.close()
pulser_data = {int(a[0]):float(a[1]) for a in dat}

fil = open('quspin_compare.dat','r')
dat = [line.split(',') for line in fil.readlines()]
fil.close()
quspin_data = {int(a[0]):float(a[1]) for a in dat}

fil = open('bloqade_compare.dat','r')
dat = [line.split('\t') for line in fil.readlines()]
fil.close()
bloqade_data = {int(float(a[0])):float(a[1]) for a in dat}


fig = figure(figsize=(8,6))
ax  = subplot(1,1,1)

semilogy(list(pulser_data.keys()),list(pulser_data.values()),'o-',label='Pulser')
semilogy(list(quspin_data.keys()),list(quspin_data.values()),'o-',label='Quspin')
semilogy(list(bloqade_data.keys()),list(bloqade_data.values()),'o-',label='Bloqade')

xlabel('Number of qubits')
ylabel('Time (sec)')
legend()
ax.set_xticks(range(4,19,2))
axis([4,18,axis()[2],axis()[3]])
tight_layout()
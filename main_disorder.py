from mpi4py import MPI

import numpy as np
import pickle
import module_disorder as Maj

# In[3]:

comm = MPI.COMM_WORLD;
rank = comm.Get_rank();

file = open('variables_83.p', 'rb')
NS_dict = pickle.load(file)
file.close()

#q = 0;
# ======== G as a function of Vz ==============
#voltageMin = -0.3; voltageMax = 0.3; voltageNumber = 601;
#voltageRange = np.linspace(voltageMin, voltageMax, voltageNumber);

#VzStep = 0.005; NS_dict['Vz'] = rank*VzStep;

#if NS_dict['VZC']=='yes':
#        Delta_0 = NS_dict['Delta_0'];
#        NS_dict['Delta_0'] = Delta_0*np.sqrt(1 - (NS_dict['Vz']/NS_dict['Vzc'])**2);
#        print(NS_dict['Delta_0']);    

#gFile = open('G_L_rank'+ str(rank)+'.txt','w');
#for voltage in voltageRange:
#    NS_dict['voltage']=voltage;
#    gFile.write( str(Maj.conductance(NS_dict)) + ',' );
#gFile.write('\n');
#gFile.close();

# ======== G as a function of Tunneling Barrier ==============
#voltageMin = -0.15; voltageMax = 0.15; voltageNumber = 3001;
#voltageRange = np.linspace(voltageMin, voltageMax, voltageNumber);

#TB_Step = 25; NS_dict['Ebarrier'] = rank*TB_Step;

#NS_dict['Vz'] = 0.99;
#if NS_dict['VZC']=='yes':
#        Delta_0 = NS_dict['Delta_0'];
#        NS_dict['Delta_0'] = Delta_0*np.sqrt(1 - (NS_dict['Vz']/NS_dict['Vzc'])**2);
#        print(NS_dict['Delta_0']);    

#gFile = open('G_TB_Vz=' + str(NS_dict['Vz']) + '_rank'+ str(rank)+'.txt','w');
#for voltage in voltageRange:
#    NS_dict['voltage']=voltage;
#    gFile.write( str(Maj.conductance(NS_dict)) + ',' );
#gFile.write('\n');
#gFile.close();

# ======== TV as a function of Vz ==============
VzMin = 0; VzMax = 1.2; VzNumber = 241;
VzRange = np.linspace(VzMin, VzMax, VzNumber);

tvFile = open('TV_ugly_14.txt','w');
for voltage in VzRange:
    NS_dict['Vz']=voltage;
    temp = Maj.TV(NS_dict);
    if abs(temp.imag) < 10**(-5):
            tvFile.write( str(temp.real) + ',');
tvFile.write('\n');
tvFile.close();

# In[ ]:

# ======== G as a function of mu ==============

#voltageMin = -0.3; voltageMax = 0.3; voltageNumber = 601;
#voltageRange = np.linspace(voltageMin, voltageMax, voltageNumber);

#mu0=0; muStep = 0.01; NS_dict['mu'] = mu0 + rank*muStep;

#gFile = open('G_Vz'+str(NS_dict['Vz'])+'_L='+str(NS_dict['wireLength'])+'_rank'+ str(rank)+'.txt','w');
#for voltage in voltageRange:
#    NS_dict['voltage']=voltage;
#    gFile.write( str(Maj.conductance(NS_dict)) + ',' );
#gFile.write('\n');
#gFile.close();







import numpy as np
import pickle

wireLength = 500;
mu_std = 1.0;
V_imp = np.random.normal(0, mu_std, wireLength);
Delta_std = 0.06;
Delta_imp = np.random.normal(0, Delta_std, wireLength);
g_std = 0.8;
g = np.random.normal(1, g_std, wireLength);

NS_dict = {'alpha':3.0,'Delta_0':0.2,'wireLength':wireLength,'t':25.0, 'mu_lead':25.0, 'Nbarrier':2, 'Ebarrier':10.0, 'QD':'no', 'QD2':'no', 'VD1':1.3, 'VD2':2.3, 'dotLength1':25, 'dotLength2':15, 'inHM_0':'no', 'inHM_SE':'no', 'Vmax':1.2, 'sigma': 40, 'disorder_mu_0':'no', 'disorder_mu_SE':'no', 'V_imp': V_imp, 'disorder_Delta_0':'no', 'disorder_Delta_SE':'no', 'Delta_imp': Delta_imp, 'disorder_g_0':'no', 'disorder_g_SE':'no', 'g': g,'SE':'yes', 'VZC':'yes','Vzc':2.5,'Vz':2.0, 'voltage':0.0, 'varymu':'no','mu':0.5,'lamd':1.0,'gamma':0.00001};

file = open('variables.p','wb')
pickle.dump(NS_dict, file)
file.close()
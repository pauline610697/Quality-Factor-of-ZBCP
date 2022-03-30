import numpy as np
import matplotlib.pyplot as plt
import pickle

# In[1]: Load Data
Vz = 0.8;
Eb = 10;
F1 = np.array([11.60,3.121,2.736,1.624,1.095,0]); # epsilon = 0.05
F2 = np.array([23.08,4.678,4.098,2.304,1.419,1.141]); # epsilon = 0.10
F3 = np.array([49.33,7.545,6.737,3.610,2.058,1.556]); # epsilon = 0.20
T = np.array([0,0.0001, 0.0005,0.001,0.002,0.003]);

# In[2]: Save Data
#file = open('F_series_Ugly_19.p','wb')
#pickle.dump(Vz, file)
#pickle.dump(Eb, file)
#pickle.dump(F1, file)
#pickle.dump(F2, file)
#pickle.dump(F3, file)
#pickle.dump(T, file)
#file.close()

# In[2]: Load Data
file = open('F_series_Ugly_18a.p', 'rb');
Vz = pickle.load(file);
Eb = pickle.load(file);
F1 = pickle.load(file);
F2 = pickle.load(file);
F3 = pickle.load(file);
T = pickle.load(file);
file.close()
T = T*(10**4);

# In[3]: Plot F versus T
plt.figure(dpi = 400)
#plt.plot(T,F1, marker = 'o')
plt.plot(T,F1, marker = 'o', label = "$\epsilon=0.05$")
plt.plot(T,F2, marker = 'o', label = "$\epsilon=0.10$")
plt.plot(T,F3, marker = 'o',color="green", label = "$\epsilon=0.20$")
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('$T$ (mK)', fontsize=16)
plt.ylabel('$F$', fontsize=16)
plt.legend(fontsize=14)
plt.title('$V_z=$'+str(Vz)+' meV.', fontsize=16)
#plt.title('$E_{barrier}=$'+str(Eb)+' meV.', fontsize=16)
plt.xlim(2,30)
#plt.ylim(-0.2,2.13)
#plt.xticks([0,0.0005,0.001, 0.0015,0.002])
#plt.savefig('/Users/laiyihua/Google Drive/UMD/Research/Majorana/merit of Quantization/FvsT_ugly_19.png')
#plt.savefig('/Users/laiyihua/Google Drive/UMD/Research/Majorana/merit of Quantization/Fcom_bad_10.png')
plt.show()
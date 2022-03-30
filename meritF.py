import numpy as np
import matplotlib.pyplot as plt

Vz = 0.8;
V0 = 0;
V1 = -0.25;
V2 = 0.25;
epsilon = 0.05;

y0 = np.loadtxt('G_TG_Vz=' + str(Vz) + '_V=' + str(V0) + '.txt', delimiter=',',dtype='str')
y1 = np.loadtxt('G_TG_Vz=' + str(Vz) + '_V=' + str(V1) + '.txt', delimiter=',',dtype='str')
y2 = np.loadtxt('G_TG_Vz=' + str(Vz) + '_V=' + str(V2) + '.txt', delimiter=',',dtype='str')

TG_Min = 0; TG_Max = 200; TG_Number = 2001;
TG_Range = np.linspace(TG_Min, TG_Max, TG_Number);

G_z = [];
G_N = [];
G_N1 = [];
G_N2 = [];
temp = [];
for i in range(TG_Number):
        G_z.append(float(y0[i])/2); # Change the conductance to be in unit of (2e^2/h)
        G_N.append((float(y1[i]) + float(y2[i]))/4) # Average and Change the conductance to be in unit of (2e^2/h)
        G_N1.append(float(y1[i])/2);
        G_N2.append(float(y2[i])/2);
        temp.append(np.absolute(float(y0[i])/2 - 1) - epsilon);

sign_change_index = []
for idx in range(0, len(temp)-1):
        # Checking for successive opposite index
        if (temp[idx] > 0 and temp[idx+1] < 0) or (temp[idx] < 0 and temp[idx+1] > 0):
                sign_change_index.append(idx)

GN_min = np.minimum(G_N[sign_change_index[0]],G_N[sign_change_index[-1]])
GN_max = np.maximum(G_N[sign_change_index[0]],G_N[sign_change_index[-1]])

F = GN_max/GN_min;

# In[1]:
# ======== G as a function of Tunneling Barrier (Linecuts together)==============      
#plt.figure(dpi = 600)
#plt.plot(TG_Range, G_z, label = "$V_{bias}=$"+str(V0))
#plt.plot(TG_Range, G_N1, label = "$V_{bias}=$"+str(V1))
#plt.plot(TG_Range, G_N2, label = "$V_{bias}=$"+str(V2))
#plt.legend()
#plt.ylim(0,1.1)
#plt.xlim(200,0)
#plt.xlabel('$E_{barrier}(meV)$')
#plt.ylabel('$G(2e^2/h)$')
#plt.title('$V_z=$'+str(Vz)+' meV, $T=0$.')
#plt.savefig('/Users/laiyihua/Google Drive/UMD/Research/Majorana/merit of Quantization/TG_good_2_Vz='+str(Vz)+'_linecuts.png')
#plt.show()

# In[2]:
# ======== G as a function of Tunneling Barrier (average G_N)==============        
plt.figure(dpi = 600)
plt.plot(G_N,G_z)
plt.ylim(0,1.1)
plt.xlabel('$G_N(2e^2/h)$')
plt.ylabel('$G_z(2e^2/h)$')
plt.title('$V_z=$'+str(Vz)+' meV, $T=0$.')
plt.savefig('/Users/laiyihua/Google Drive/UMD/Research/Majorana/merit of Quantization/TG_ugly_11_Vz='+str(Vz)+'_GzGN.png')
plt.show()

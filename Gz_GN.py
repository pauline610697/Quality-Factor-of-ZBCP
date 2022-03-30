import numpy as np
import matplotlib.pyplot as plt

# In[1]: Parameters
Ebarrier = 10.0;
Vz = 0.8;
T_vec = [0.0001, 0.001,0.003,0.005]; # meV

Vmin = -0.3; Vmax = 0.3; Vpoints = 601;
Vstep = (Vmax - Vmin)/ (Vpoints - 1);
V0 = 0; V1 = -0.25; V2 = 0.25;
V_num0 = int((V0 - Vmin)/ Vstep); V_num1= int((V1 - Vmin)/ Vstep); V_num2 = int((V2 - Vmin)/ Vstep);

TG_Number = 1001;

# In[2]: Load data_Ebarrier and Plot Gz as a function of G_N
plt.figure(dpi = 400)
G_z = []; G_n = [];
temp1 = [];
for j in range(TG_Number): # T=0 case
        y2 = np.loadtxt('G_set=33_TB_Vz='+str(Vz)+'_rank'+ str(j)+'.txt', delimiter=',',dtype='str');        
        G_z.append(float(y2[V_num0])/2);
        G_n.append((float(y2[V_num1]) + float(y2[V_num2]))/4);
        print(j);
plt.plot(G_n,G_z, label = "$T=0$ mK.")

for T in T_vec:
        G_z = []; G_n = [];
        temp1 = [];
        for j in range(TG_Number): 
                y2 = np.loadtxt('G_TB_T='+str(T)+'_Vz='+str(Vz)+'_rank'+ str(j)+'.txt', delimiter=',',dtype='str');
                G_z.append(float(y2[V_num0])/2);
                G_n.append((float(y2[V_num1]) + float(y2[V_num2]))/4);               
                print(j);
        plt.plot(G_n,G_z, label = "$T=$"+str(T*10**4)+" mK.")

plt.legend(fontsize=14)                
#plt.ylim(0,1.2)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('$G_N$ $(2e^2/h)$',fontsize=16)
plt.ylabel('$G_z$ $(2e^2/h)$',fontsize=16)
plt.title('$V_z=$'+str(Vz)+' meV.', fontsize=16)
#plt.savefig('/Users/laiyihua/Google Drive/UMD/Research/Majorana/merit of Quantization/TG_ugly_19_Vz='+str(Vz)+'_GzGN_Tdep.png')
plt.show()
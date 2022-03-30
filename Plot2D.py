import numpy as np
import matplotlib.pyplot as plt

Vz = 0.8;
rankNumber = 241; # 400
VzStep = 0.005; # 0.005 or 0.1 or 0.2
Vpoints = 601;
colorbarL = 0;
colorbarU = 2.0;
Gmatrix = []

y = np.linspace(-0.3,0.3,Vpoints)
x = np.linspace(0,VzStep*(rankNumber-1),rankNumber)

rows = rankNumber;
cols = Vpoints;
for i in range(rows):
    row = []
    y1 = np.loadtxt('G_set=33_rank'+ str(i)+'.txt', delimiter=',',dtype='str')
    #y1 = np.loadtxt('G_L_rank'+ str(i)+'.txt', delimiter=',',dtype='str')
    #y1 = np.loadtxt('G_L_T=0.001_rank'+ str(i)+'.txt', delimiter=',',dtype='str')
    #y1 = np.loadtxt('G_set=33_TB_Vz='+str(Vz)+'_rank'+ str(i)+'.txt', delimiter=',',dtype='str')
    #y1 = np.loadtxt('G_TB_T=0.005_Vz='+str(Vz)+'_rank'+ str(i)+'.txt', delimiter=',',dtype='str')
    for j in range(cols):
        row.append(float(y1[j]))
    Gmatrix.append(row)

# Transpose the G matrix
result = [[0]*rows]*cols;
result = [[Gmatrix[j][i] for j in range(len(Gmatrix))] for i in range(len(Gmatrix[0]))]

# In[B]: Plot
#plt.figure(dpi = 600)
fig, ax = plt.subplots(dpi=400)
im = ax.pcolor(x,y,result, cmap='RdBu_r')
#im = ax.pcolor(x,y,result,vmin=0,vmax=2, cmap='RdBu_r')
fig.colorbar(im,ax=ax)
#plt.pcolor(x,y,result)
im.set_clim(colorbarL, colorbarU) # Set the axis scale

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.xlabel('$E_{barrier}$ (meV)', fontsize=14)
plt.xlabel('$V_z$ (meV)', fontsize=14)
plt.ylabel('$V$ (mV)', fontsize=14)
#plt.xlim(200,0)
#plt.ylim(-0.25,0.25)
plt.axvline(x=Vz, color='yellow', linestyle='--')
#plt.title('$G$ $(e^2/h)$, $V_z=$'+str(Vz)+' meV.', fontsize=14)
plt.title('$G$ $(e^2/h)$', fontsize=14)
#plt.savefig('/Users/laiyihua/Google Drive/UMD/Research/Majorana/merit of Quantization/zero-bias_ugly_19.png')
#plt.savefig('/Users/laiyihua/Google Drive/UMD/Research/Majorana/merit of Quantization/TG_ugly_19_Vz='+str(Vz)+'.png')
plt.show()
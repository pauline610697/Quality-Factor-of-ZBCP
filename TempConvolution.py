import numpy as np

rankNumber = 1001;
Vpoints = 6001;
voltageMin = -0.3; voltageMax = 0.3;
T = 0.005; # Unit: meV
Vz = 1.5;

E = np.linspace(voltageMin,voltageMax,Vpoints);
dE = E[1]-E[0];
f = 1/(np.exp(E/T)+1);
#dfdE = np.gradient(f,dE);
dfdE = -1/(4*T*(np.cosh(E/(2*T)))**2);

for i in range(rankNumber):
        #y1 = np.loadtxt('G_set=0_rank'+ str(i)+'.txt', delimiter=',',dtype='str'); # Read G from files. dataType=String
        y1 = np.loadtxt('G_TB_Vz=' + str(Vz) + '_rank'+ str(i)+'.txt', delimiter=',',dtype='str'); # Read G from files. dataType=String
        
        y1 = ','.join(y1); # Convert String arrays into one single string
        y1=np.fromstring(y1,dtype=np.float, sep=','); # Convert sting into Float arrays
        G = -np.convolve(y1,dfdE, 'same')*dE;
        #gFile = open('G_L_T='+str(T)+'_rank'+str(i)+'.txt','w');    
        gFile = open('G_TB_T='+str(T)+'_Vz=' + str(Vz) + '_rank'+ str(i)+'.txt','w');    
        for j in range(Vpoints):
                gFile.write( str(G[j]) + ',');
        gFile.write('\n');
        gFile.close();
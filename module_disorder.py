import numpy as np
import kwant
from math import pi
from cmath import sqrt
import numpy.linalg as LA

# =========== global variables ====================
s0 = np.array([[1.0, 0.0], [0.0, 1.0]]); sx = np.array([[0.0, 1.0], [1.0, 0.0]]); 
sy = np.array([[0.0, -1j], [1j, 0.0]]); sz = np.array([[1.0, 0.0], [0.0, -1.0]]);

t0 = np.array([[1.0, 0.0], [0.0, 1.0]]); tx = np.array([[0.0, 1.0], [1.0, 0.0]]);
ty = np.array([[0.0, -1j], [1j, 0.0]]); tz = np.array([[1.0, 0.0], [0.0, -1.0]]);

tzs0 = np.kron(tz,s0); t0s0 = np.kron(t0,s0); t0sx = np.kron(t0,sx);
txs0 = np.kron(tx,s0); tzsy = np.kron(tz,sy); t0sz = np.kron(t0,sz);


# In[327]:

def NSjunction(args_dict):
    t = args_dict['t'];
    alpha = args_dict['alpha'];
    Vz = args_dict['Vz'];
    Delta_0 = args_dict['Delta_0'];
    mu = args_dict['mu'];
    mu_lead = args_dict['mu_lead'];
    wireLength = args_dict['wireLength'];
    Nbarrier = args_dict['Nbarrier'];
    Ebarrier = args_dict['Ebarrier'];
    gamma = args_dict['gamma'];
    lamd=args_dict['lamd'];
    voltage=args_dict['voltage'];
    
    #========= set-up of NS-junction =========
    junction = kwant.Builder(); a=1; # Lattice constant
    lat = kwant.lattice.chain(a);     

    for x in range(wireLength):
        junction[ lat(x) ] = (2*t - mu)*tzs0 + Vz*t0sx + Delta_0*txs0 - 1j*gamma*t0s0;
    
    if args_dict['varymu']=='yes':
        mu1=args_dict['mu1']
        mu2=args_dict['mu2']
        uncoverLength=args_dict['uncoverLength']
        mus=np.linspace(mu1,mu2,wireLength)
        for x in range(wireLength):
            junction[ lat(x) ] = (2*t -mus[x])*tzs0 + Vz*t0sx + Delta_0*txs0 - 1j*gamma*t0s0
        for x in range(uncoverLength):
            junction[ lat(x) ] = (2*t -mus[x])*tzs0 + Vz*t0sx - 1j*gamma*t0s0

    if args_dict['SE'] == 'yes':
        SelfE=lamd*np.sign(voltage-Delta_0)*(voltage*t0s0+Delta_0*txs0)/sqrt(Delta_0**2-voltage**2+1e-9j);
        for x in range(wireLength):
            junction[ lat(x) ] = (2*t - mu)*tzs0 + Vz*t0sx + SelfE - 1j*gamma*t0s0;
	
    if args_dict['QD'] == 'yes':
        dotLength1 = args_dict['dotLength1'];
        VD1 = args_dict['VD1'];
        for x in range(0,dotLength1):
           junction[ lat(x) ] = (2*t - mu + VD1*np.exp(-(x/dotLength1)**2))*tzs0 + Vz*t0sx - 1j*gamma*t0s0;
	
    if args_dict['QD2'] == 'yes': #Quantum Dot away from the lead.
            dotLength2 = args_dict['dotLength2'];
            VD2 = args_dict['VD2'];
            for x in range(0,dotLength2):
                    junction[ lat(wireLength - x) ] = (2*t - mu + VD2*np.exp(-(x/dotLength2)**2))*tzs0 + Vz*t0sx - 1j*gamma*t0s0;

    if args_dict['inHM_0'] == 'yes':
        sigma = args_dict['sigma'];
        Vmax = args_dict['Vmax'];
        for x in range(wireLength):
           junction[ lat(x) ] = (2*t - mu + Vmax*np.exp(-0.5*(x/sigma)**2))*tzs0 + Vz*t0sx + Delta_0*txs0 - 1j*gamma*t0s0;
    
    if args_dict['inHM_SE'] == 'yes':
        sigma = args_dict['sigma'];
        Vmax = args_dict['Vmax'];
        SelfE=lamd*np.sign(voltage-Delta_0)*(voltage*t0s0+Delta_0*txs0)/sqrt(Delta_0**2-voltage**2+1e-9j);
        for x in range(wireLength):
           junction[ lat(x) ] = (2*t - mu + Vmax*np.exp(-0.5*(x/sigma)**2))*tzs0 + Vz*t0sx + SelfE - 1j*gamma*t0s0;

    if args_dict['disorder_mu_0'] == 'yes':
        V_imp = args_dict['V_imp'];
        for x in range(wireLength):
           junction[ lat(x) ] = (2*t - mu + V_imp[x])*tzs0 + Vz*t0sx + Delta_0*txs0 - 1j*gamma*t0s0;
    
    if args_dict['disorder_mu_SE'] == 'yes':
        V_imp = args_dict['V_imp'];
        SelfE=lamd*np.sign(voltage-Delta_0)*(voltage*t0s0+Delta_0*txs0)/sqrt(Delta_0**2-voltage**2+1e-9j);
        for x in range(wireLength):
           junction[ lat(x) ] = (2*t - mu + V_imp[x])*tzs0 + Vz*t0sx + SelfE- 1j*gamma*t0s0;
           
    if args_dict['disorder_Delta_0'] == 'yes':
        Delta_imp = args_dict['Delta_imp'];
        for x in range(wireLength):
           junction[ lat(x) ] = (2*t - mu)*tzs0 + Vz*t0sx + (Delta_0 + Delta_imp[x])*txs0 - 1j*gamma*t0s0;    

    if args_dict['disorder_Delta_SE'] == 'yes':
        Delta_imp = args_dict['Delta_imp'];
        for x in range(wireLength):
           SelfE=lamd*np.sign(voltage- (Delta_0 + Delta_imp[x]))*(voltage*t0s0+(Delta_0 + Delta_imp[x])*txs0)/sqrt((Delta_0 + Delta_imp[x])**2-voltage**2+1e-9j);     
           junction[ lat(x) ] = (2*t - mu)*tzs0 + Vz*t0sx + SelfE - 1j*gamma*t0s0;               

    if args_dict['disorder_g_0'] == 'yes':
        g = args_dict['g'];
        for x in range(wireLength):
           junction[ lat(x) ] = (2*t - mu)*tzs0 + g[x]*Vz*t0sx + Delta_0*txs0 - 1j*gamma*t0s0;

    if args_dict['disorder_g_SE'] == 'yes':
        g = args_dict['g'];
        SelfE=lamd*np.sign(voltage-Delta_0)*(voltage*t0s0+Delta_0*txs0)/sqrt(Delta_0**2-voltage**2+1e-9j);
        for x in range(wireLength):
           junction[ lat(x) ] = (2*t - mu)*tzs0 + g[x]*Vz*t0sx + SelfE - 1j*gamma*t0s0;
                                             
    for x in range(Nbarrier):
        junction[ lat(x) ] = (2*t - mu + Ebarrier)*tzs0 + Vz*t0sx - 1j*gamma*t0s0;
    
    for x in range( 1, wireLength ):
        junction[ lat(x-1), lat(x) ] = -t*tzs0 - 1j*alpha*tzsy;
    
    symLeft = kwant.TranslationalSymmetry([-a]);
    lead = kwant.Builder(symLeft);
    #lead = kwant.Builder(symLeft, particle_hole=ty);
    lead[ lat(0) ] = (2*t - mu_lead)*tzs0 + Vz*t0sx;
    lead[ lat(0), lat(1) ] = -t*tzs0 - 1j*alpha*tzsy;
    junction.attach_lead(lead);
    junction = junction.finalized();
    
    return junction;

def conductance(args_dict):
    voltage = args_dict['voltage'];
    junction = NSjunction(args_dict);
    S_matrix = kwant.smatrix(junction, voltage, check_hermiticity=False);
    R = S_matrix.submatrix(0,0); 
    
    G = 2.0;
    for (i,j) in [(0,0),(0,1),(1,0),(1,1)]:
        G = G - abs(R[i,j])**2 + abs(R[2+i,j])**2; ##

    return G;

def TV(args_dict):
    args_dict['voltage'] = 0.0; 
    junction = NSjunction(args_dict);
    
    S_matrix = kwant.smatrix(junction, args_dict['voltage'], check_hermiticity=False);
    R = S_matrix.submatrix(0,0); # "0" for The first lead index
    tv0 = LA.det(R);
    basis_wf = S_matrix.lead_info[0].wave_functions; # 'lead_info[i].wave_functions' contains the wavefunctions of the propagating modes in lead "i"
    #This is a 2D array, where the first axis corresponds to the orbitals in a unit cell of the lead, and the second axis is the mode number.
    normalize_dict = {0:0,1:0,2:3,3:3,4:0,5:0,6:3,7:3};
    phase_dict = {};
    
    for n in range(8):
        m = normalize_dict[n];
        phase_dict[n] = (-1)**m*basis_wf[m,n]/abs(basis_wf[m,n]); ##
        #phase_dict[n] = basis_wf[m,n]/abs(basis_wf[m,n]);
        
    tv = tv0*np.conjugate(phase_dict[0]*phase_dict[1]*phase_dict[2]*phase_dict[3])*phase_dict[4]*phase_dict[5]*phase_dict[6]*phase_dict[7];
    
    return tv
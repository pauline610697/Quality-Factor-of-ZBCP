function rho=dos_disorderMuH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V_imp,V,s)
H = disorderMuH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V_imp,V+1i*s);
temp = (V + 1i*s)*speye(4*N_tot);
G = inv(temp - H);
rho = -imag(trace(G))./pi;
end
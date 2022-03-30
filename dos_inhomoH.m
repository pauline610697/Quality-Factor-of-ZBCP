function rho=dos_inhomoH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,sigma,V,s)
H = inhomoH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,sigma,V+1i*s);
temp = (V + 1i*s)*speye(4*N_tot);
G = inv(temp - H);
rho = -imag(trace(G))./pi;
end
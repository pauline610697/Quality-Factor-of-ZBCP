function rho=dos_goodH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V,s)
H = goodH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V+1i*s);
temp = (V + 1i*s)*speye(4*N_tot);
G = inv(temp - H);
rho = -imag(trace(G))./pi;
end
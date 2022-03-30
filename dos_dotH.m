function rho=dos_dotH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,N_dot,V,s)
H = dotH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,N_dot,V+1i*s);
temp = (V + 1i*s)*speye(4*N_tot);
G = inv(temp - H);
rho = -imag(trace(G))./pi;
end
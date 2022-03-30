function H = dotH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,N_dot,omega)
%% Spin & Particle-hole spaces (Pauli matrices)
s0 = [1 0; 0 1]; sx = [0 1; 1 0]; sy = [0 -1i; 1i 0]; sz = [1 0; 0 -1]; % Spin space
t0 = [1 0; 0 1]; tx = [0 1; 1 0]; ty = [0 -1i; 1i 0]; tz = [1 0; 0 -1]; % particle-hole space

tzs0 = kron(tz,s0); tzsy = kron(tz,sy); t0s0 = kron(t0,s0); txs0 = kron(tx,s0); t0sx = kron(t0,sx);

band11sm = spdiags([ones(N_tot,1) ones(N_tot,1)],[-1,1],N_tot,N_tot);
band1m1sm = spdiags([ones(N_tot,1) -ones(N_tot,1)],[-1,1],N_tot,N_tot);
eyesm = speye(N_tot);

%% Extra Potentials
barrier = spdiags(cat(1,ones(N_barrier,1),zeros(N_tot-N_barrier,1)),0,N_tot,N_tot);

dot_vec1 = (N_barrier+1):N_dot;
dotV1 = VD*exp(-(dot_vec1./N_dot).^2)';
dot1 = spdiags(cat(1,zeros(N_barrier,1),dotV1,zeros(N_tot-N_dot,1)),0,N_tot,N_tot);
NOT_dot = spdiags(cat(1,zeros(N_dot,1),ones(N_tot-N_dot,1)),0,N_tot,N_tot); % QD on the left.

%%
H = zeros(4*N_tot,4*N_tot);
SelfE = -lambda.*((omega)*t0s0 + Delta*txs0)./sqrt(Delta.^2 -(omega).^2);
H = kron(eyesm, (2*t - mu)*tzs0 + Vz*t0sx) + kron(barrier,E_barrier*tzs0) + kron(dot1,tzs0) + kron(NOT_dot, SelfE); % Diagonal Part --- Nanowire
%H = kron(eyesm, (2*t - mu)*tzs0 + Vz*t0sx) + kron(barrier,E_barrier*tzs0) + kron(dot1,tzs0) + kron(NOT_dot, Delta*txs0);

H = H + kron(band11sm, -t*tzs0) + kron(band1m1sm, 1i*alpha*tzsy); % Off-diagonal Part
end
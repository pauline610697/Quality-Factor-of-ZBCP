%% Version 1: This code is to calculate the wave function inside the nanowire for a specific Vz.
% 
clear;
tic;
%% Parameters Setting
% Note that the length scale is in unit of lattice constant, which is 10nm.

% Basic Parameters (nanowire + lead)
t = 25; %unit: meV
Delta_0 = 0.2; %unit: meV
wireLength = 500; %unit: 10nm
alpha = 2.5; %unit: meV/m
E_barrier = 10; %unit: meV
N_barrier = 1; %unit: 10nm
lambda = 0.2; %unit: meV
Vc = 1.2; %unit: meV
mu = 1.0; %unit: meV
Vz = 0.8;

% Quantum Dot Parameters
VD = 0.3; %unit:meV
sigma = 40;  %unit: 10nm
dotLength = 15; %unit: 10nm
s = 1e-2;
Type = "disorderMu"; % 'disorderMu', 'inhomo', 'dot', 'good'

if Type=="disorderMu"
    fileID = fopen('mu_imp.txt','r');
    V_imp = fscanf(fileID,'%f');
end

%% Construct the Vectors
N_tot = wireLength;

VzMin = 0; VzMax = 1.2; VzNumber = 241;
VzStep = (VzMax - VzMin)./(VzNumber - 1);
VzRange = linspace(VzMin,VzMax,VzNumber);

Vmin = -0.3; Vmax = 0.3; Vnumber = 601;
Vstep = (Vmax - Vmin)./(Vnumber - 1);
Vrange = linspace(Vmin,Vmax,Vnumber);

rho = zeros(VzNumber,Vnumber);
dosmap = cell(1,VzNumber);

%% Calculate the total DOS
parfor k = 1:VzNumber
    Vz = VzMin + (k-1).*VzStep;
    disp(k);
    Delta = Delta_0.*sqrt(1 - (Vz./Vc).^2).*(Vz<Vc);
    %Delta = Delta_0;

    if Type=="disorderMu"
        rho(k,:) = arrayfun(@(V) dos_disorderMuH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V_imp,V,s), Vrange);
    elseif Type=="inhomo"
        rho(k,:) = arrayfun(@(V) dos_inhomoH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,sigma,V,s), Vrange);
    elseif Type=="dot"
        rho(k,:) = arrayfun(@(V) dos_dotH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,dotLength,V,s), Vrange);
    elseif Type=="good"
        rho(k,:) = arrayfun(@(V) dos_goodH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V,s), Vrange);
    end
    
    [~,loc] = findpeaks(rho(k,:));
    Vloc = Vrange(loc);
    dosmap{k} = Vloc;        
end

%save Energy_data_1.mat
%% Plot Energy Spectrum
figure()
for k = 1:VzNumber
    Vz = VzMin + (k-1).*VzStep;
    scatter(ones(1,length(dosmap{k}))*Vz,dosmap{k},'b','.');
    hold on
end

box on
hold off

%title('$t=25$ meV, $\Delta_0=0.2$ meV, $V_c=1.2$ meV, $\alpha=2.5$ meV, $\mu=1.0$ meV, $\lambda=0.2$ meV, $L=3.0$ um.','interpreter','latex','FontSize',16)
axis([VzMin VzMax Vmin Vmax])
TQPT = sqrt(lambda^2 + mu^2);
line([TQPT,TQPT],[Vmin,Vmax],'Color','black','LineStyle','--','LineWidth',0.75)
%line([Vc,Vc],[Vmin,Vmax],'Color','black','LineStyle','--','LineWidth',0.75)
xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
ylabel('$$E$$ (meV)','interpreter','latex','FontSize',16)
set(gca,'FontSize',16)
xticks([1 2 3 4 5])
%% Wave-function
Delta = Delta_0.*sqrt(1 - (Vz./Vc).^2).*(Vz<Vc);

if Type=="disorderMu"
    Rho = arrayfun(@(V) dos_disorderMuH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V_imp,V,s), Vrange);
elseif Type=="inhomo"
    Rho = arrayfun(@(V) dos_inhomoH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,sigma,V,s), Vrange);
elseif Type=="dot"
    Rho = arrayfun(@(V) dos_dotH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,dotLength,V,s), Vrange);
elseif Type=="good"
    Rho = arrayfun(@(V) dos_goodH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V,s), Vrange);
end

[~,loc] = findpeaks(Rho);
Vloc = Vrange(loc);

if mod(length(Vloc),2)==0
    firstE = Vloc(1+length(Vloc)./2);
else
    firstE = Vloc((length(Vloc)+1)./2);
end

if Type=="disorderMu"
    H1 = disorderMuH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V_imp,firstE);
    H2 = disorderMuH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,V_imp,-firstE);
elseif Type=="inhomo"
    H1 = inhomoH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,sigma,firstE);
    H2 = inhomoH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,sigma,-firstE);
elseif Type=="dot"
    H1 = dotH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,dotLength,firstE);
    H2 = dotH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,VD,dotLength,-firstE);
elseif Type=="good"
    H1 = goodH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,firstE);
    H2 = goodH(t,alpha,mu,Delta,Vz,lambda,N_tot,E_barrier,N_barrier,-firstE);
end

[WF1,D1] = eigs(H1,6,'SM');
D1 = diag(D1);
WF1 = WF1(:,D1>=0);
D1 = D1(D1>=0);
[~,index1] = min(abs(D1-firstE));
Min1 = D1(index1);
WF_up = WF1(:,index1);

[WF2,D2] = eigs(H2,6,'SM');
D2 = diag(D2);
WF2 = WF2(:,D2<=0);
D2 = D2(D2<=0);
[~,index2] = min(abs(D2+firstE));
Min2 = D2(index2);
WF_down = WF2(:,index2);

WF_1A = (WF_up + WF_down)./sqrt(2);
WF_1B = -1i.*(WF_up - WF_down)./sqrt(2);

P_1A = zeros(1,N_tot);
P_1B = zeros(1,N_tot);

for n = 1:N_tot
    P_1A(n) = (norm(WF_1A(4*(n-1)+1:4*(n-1)+4))).^2;
    P_1B(n) = (norm(WF_1B(4*(n-1)+1:4*(n-1)+4))).^2;
end

% normalization
P1A = P_1A/(sum(P_1A));
P1B = P_1B/(sum(P_1B));

P1A = [0 P1A 0];
P1B = [0 P1B 0];
%% Plot wave-function
figure()
xRange = 0:(N_tot+1);
%xRange = 0:F:(wireLength+1).*F;
fill(xRange,P1A,'b','EdgeColor','b','FaceAlpha',0.2,'EdgeAlpha',1);
hold on;
fill(xRange,P1B,'r','EdgeColor','r','FaceAlpha',0.2,'EdgeAlpha',1);
xlim([0 wireLength])
title('Lowest-lying wave-function: $V_z=0.8$ meV','interpreter','latex','FontSize',22)
%legend({'1st lowest state, Majorana A','1st lowest state, Majorana B'},'FontSize',14)
%legend({'1st lowest state, \phi_\epsilon + \phi_{-\epsilon}','1st lowest state, -i(\phi_\epsilon - \phi_{-\epsilon})'},'FontSize',14)
%legend({'Majorana A','Majorana B'},'FontSize',14)
legend({'\phi_\epsilon + \phi_{-\epsilon}','i(\phi_\epsilon - \phi_{-\epsilon})'},'FontSize',20)
xlabel('$$x$$ (10nm)','interpreter','latex','FontSize',24)
ylabel('$$|\Psi|^2$$','interpreter','latex','FontSize',24)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontName','Times','fontsize',22)
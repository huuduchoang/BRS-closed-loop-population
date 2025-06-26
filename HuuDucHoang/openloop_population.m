function z = openloop_population(t,u,params)

%     params.N           - number of neurons
%     params.El       - N×1 vector of E_L draws
%     params.gnap        - N×1 vector of g_NaP draws
%     params.gsyn        - N×N matrix of g_syn draws (zero diagonal)
%     params.phi         - N×1 vector of chemosensory φ draws
%     params.thetaO2     - N×1 vector of θ_O2 draws
%     params.sigmaO2     - N×1 vector of σ_O2 draws
%     params.gtonic      - value of tonic conductance

%% Unpack State variables

N = params.N;
v = u(1 : N);
n = u(N+1 : 2*N);
h = u(2*N+1 : 3*N);
s = u(3*N+1 : 4*N);
alpha = u(4*N+1);
vollung = u(4*N+2);
PO2lung = u(4*N+3);
PO2blood = u(4*N+4);

%% Draw Heterogeneous Parameters
El = params.Eleak;
gnap = params.gnap;

% Fully-connected synaptic matrix with zero diagonal
gsyn = params.gsyn;

%% CPG 

% capacitance
C = 21;  

% maximal conductances
gna=28; gk=11.2; gl=2.8;

% reversal potentials
Ena=50; Ek=-85; Esyn=0;

% persistent sodium
theta_mp = -40;  sigma_mp = -6; 
theta_h = -48; sigma_h = 6; taumax_h = 10000;

mp_inf = 1 ./ (1+exp((v-theta_mp)/sigma_mp));
h_inf = 1 ./ (1+exp((v-theta_h)/sigma_h));
tau_h = taumax_h ./ cosh((v-theta_h)/(2*sigma_h));

Inap = gnap .* mp_inf .* h .* (v-Ena);

% transient sodium
theta_m = -34; sigma_m = -5;

m_inf = 1 ./ (1+exp((v-theta_m)/sigma_m));

Ina = gna .* (m_inf.^3) .* (1-n) .* (v-Ena);

% potassium
theta_n = -29; sigma_n = -4; taumax_n = 10;

Ik = gk .* (n.^4) .* (v-Ek);

n_inf = 1 ./ (1+exp((v-theta_n)/sigma_n));
tau_n = taumax_n ./ cosh((v-theta_n)/(2*sigma_n));

% leak
Il = gl .* (v-El);

% synaptic gating
tau_s = 5;
k_r = 1;
thetas = -10; % half-activation
sigmas = -5;
s_inf = 1./ (1 + exp((v - thetas)/sigmas));
g_in = gsyn * s;
Isyn = g_in .* (v - Esyn);

% tonic
Itonic = params.gtonic .* (v-Esyn);

%% Motor pool

r = 0.001; Tmax = 1; VT = 2; Kp = 5;

NT = 1 ./ (1+exp(-(v-VT)/Kp));
Tpop = Tmax * mean(NT);

%% Lung volume

E1 = 0.0025; E2 = 0.4; Vol0 = 2;

dvolrhs=max(0,-E1*(vollung-Vol0)+E2*alpha);

%% Lung oxygen

PO2ext = (760-47)*.21;  R = 62.364; Temp = 310; 

taulb = 500;

%% Blood oxygen

M = 8e-6;
Hb = 150; volblood = 5; eta = Hb*1.36; gamma = volblood/22400; betaO2 = 0.03;

c = 2.5; K = 26;
SaO2 = (PO2blood^c)/(PO2blood^c+K^c);
CaO2 = eta*SaO2+betaO2*PO2blood;
partial = (c*PO2blood^(c-1))*(1/(PO2blood^c+K^c)-(PO2blood^c)/((PO2blood^c+K^c)^2));

Jlb=(1/taulb)*(PO2lung-PO2blood)*(vollung/(R*Temp));
Jbt=M*CaO2*gamma;


%% Differential equations

z = zeros(size(u));
z(1 : N) = (-Inap-Ina-Ik-Il-Itonic-Isyn)./C;
z(N+1 : 2*N) = (n_inf-n)./tau_n;
z(2*N+1 : 3*N) = (h_inf-h)./tau_h;
z(3*N+1 : 4*N) = ((1 - s) .* s_inf  -  k_r .* s) ./ tau_s;
z(4*N+1) = r*Tpop*(1-alpha)-r*alpha;
z(4*N+2) = -E1*(vollung-Vol0)+E2*alpha;
z(4*N+3) = (1/vollung)*(PO2ext-PO2lung)*dvolrhs-Jlb*(R*Temp/vollung);
z(4*N+4) = (Jlb-Jbt)/(gamma*(betaO2+eta*partial));

z=z(:);





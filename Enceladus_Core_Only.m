%% Poroviscoelastic Enceladus model, core-only 
% Example to compute the proviscoelastic response of Enceladus. 
% We consider the model of Liao et al. 2020: https://doi.org/10.1029/2019JE006209. In such case the boundary conditions are given by Eq. (A16)
%%
close all
clear 
clc
addpath(genpath(pwd))
set(0,'defaulttextInterpreter','latex') 
%% PHYSICAL PARAMETERS 
R0=[0 191.1e3;]; %layering [m]
rhof0=1E3; %density of the fluid phase kg/m^3
mu0=1E9; %shear modulus [GPa]
Ks0=10E9; %bulk modulus solid [GPa]
Kf0=2.2E9; %bulkd modulus fluid [GPa]
k_perm0=1E-8; %fluid permeability [m^2]
etaf0=1.9E-3; % viscosity of the mantle's liquid phase
liquid=0; %is the layer liquid 
poro=0.2; %porosity
alpha=0.95;
etas0=1e16; % viscosity of the mantle's solid phase
omega0= 2*pi/(33*3600); %Encelauds
eccen=0.0047; % eccentricity Enceladus
l=2; %degree of the forcing
Gg0=6.6743E-11; %gravitational constant
% set average density of the core
rhoc=2422;
rhos0(1)=(rhoc-poro(1)*rhof0(1))/(1-poro(1));
rho0=(1-poro).*rhos0+poro.*(rhof0);

%% Numerical parameters 
nrr=400; %number of points used for the radial integration
resample=1;
%% Non-dimentionalize 
T_nd=2*pi/omega0; %time scale, we use Io's orbital period.
R=R0/R0(end);
rho=rho0/rho0(end);
rhos=rhos0/rho0(end);
rhof=rhof0/rho0(end);
Ks=Ks0/mu0(end);
Kf=Kf0/mu0(end);
mu=mu0/mu0(end);
etas=etas0/(mu0(end)*T_nd);
etaf=etaf0/(mu0(end)*T_nd);
k_perm=k_perm0/(R0(end)^2);
global Gg
Gg=Gg0*rho0(end)^2*R0(end)^2/mu0(end);
omega=omega0*T_nd; 
%% Compute radial functions 
% compute tidal response for a moon with an internal ocean layer in terms of the y functions defined in Appendix A. Note that variables are non-dimensional 
% compute the imposed strain 
gs=0.125;
G=6.6743e-11;
g=0.125;
Mu_star=1/(1/mu0+1/(1j*omega0*etas0));
epsilon_0=1;
% The maximum displacement imposed at the pole is given by (see Section 3.1):
epsilon_pole=abs((9/(4*pi))*eccen*(omega0^2/(Gg0*rho0))*5/3*(3/2)/(1+19*Mu_star./(2*rho0*gs*R0(end))));
% obtain y(r)functions. Note in this case the following options are used: 
    % (1) self_gravity: self gravity is turned off 
    % (2) tidal_fluid: The fluid is not directly affected by the tidal force, it is massless 
    % (3) pressure_BC: At the surface of the core, the boundary conditions \sigma_rr+p=0, is used
    % (4) strain_BC: A strain is prescribed as boundary condition
[y, r]=tidal(l,R,rho,rhof,mu,Ks,etas,alpha,poro,k_perm,etaf,Kf,liquid,omega,'self_gravity',0,'tidal_fluid',0, 'radial_points',nrr,'pressure_BC',0,'strain_BC',epsilon_0,'gravity_on',0,'resample',resample);
%% Obtain the different relevant fields, strain, stress, etc. 
%The different fields can be obtained as indicated in Appendix (D).
[colat,lon,rr,displacements_20,flux_20,strain_20,stress_20,gravpot_20,p_fluid_20,C_fluid_20]=build_solution(y,r,R,rhof,rho,Ks,Kf,mu,etas,etaf,k_perm,liquid,omega,alpha,poro,2,0,'tidal_fluid',0,'lat_points',300);
%(2) normalize so that the strain at the poles is given by epsilon_pole
norm=epsilon_pole/abs(strain_20(1,1,end,1));
strain=norm*strain_20;
stress=norm*stress_20;
p_fluid=norm*p_fluid_20;
C_fluid=norm*C_fluid_20;
flux=norm*flux_20;
displacements=norm*displacements_20;
gravpot=norm*gravpot_20;
%% Compute energy dissipation 
constantE=mu0/T_nd; %used to get energy in dimensional units 
[energy_solidV, energy_solid_pore,energy_fluidV,energy_solid_surface, energy_solid_pore_surface,energy_fluid_surface,energy_solid_total_out, energy_fluid_total_out]=compute_energy(strain,stress,flux,p_fluid,C_fluid,omega,etaf,k_perm,rr,colat,lon,poro);
energy_solid=constantE*R0(end)^3*energy_solid_total_out;
energy_fluid=constantE*R0(end)^3*energy_fluid_total_out;
energy_total=energy_solid+energy_fluid;
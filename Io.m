%% Viscoelastic Io Model
% Example to compute the tidal response of a viscoelastic Io. 
%The interior model of Io corresponds to models A and B of Steinke et al. 2020. https://doi.org/10.1016/j.icarus.2019.05.001. Im(k2) should give 0.015 
close all
clear 
clc
addpath(genpath(pwd))
set(0,'defaulttextInterpreter','latex') 

%% Example 2: Viscoelastic Io Model 
R0=[0 965e3 1591e3 1791e3 1821e3]; %layering [m]
rhos0=[5.150e3 3.244E3 3.244E3 3.244E3]; %density of the solid phase kg/m^3
rho0=[5.150e3 3.244E3 3.244E3 3.244E3]; %density of the solid phase kg/m^3
rhof0=[0 0 0 0]; %no fluid
Ks0=[10E40 10E20 10E20 10E20]; %bulk modulus solid [GPa], assumed incompressible 
Kf0=[2.2E9 2.2E9 2.2E9 2.2E9]; %bulkd modulus fluid [GPa]
k_perm0=[1E-7 1E-6 1E-6 1E-7]; %fluid permeability 
alpha=[ 0 0 0 0]; %Biot's coefficient 
poro=[ 0 0 0 0]; %porosity 
etaf0=[1 1 1 1]; % viscosity of the mantle's liquid phase
liquid=[1 0 0 0]; % Liquid core 
%Interior parameters for the two models, models A and B 
mu0_B=[1E9 2E9 9E7 6.5E10]; %shear modulus [GPa]
etas0_B=[10E15 8E14 3.5E12 1E23]; % viscosity of the mantle's solid phase
%astenosphere only heating
mu0_A=[1E9 6E10 7.8E5 6.5E10]; %shear modulus [GPa]
etas0_A=[10E15 1E20 1E11 1E23]; % viscosity of the mantle's solid phase
mu0=mu0_B;
etas0=etas0_B;

l=2; %degree of the forcing
Gg0=6.6743E-11; %gravitational constant
%Orbital parameters Io 
eccen=4.1e-03;% eccentricity Io
omega0=4.1086E-05; %Io
%% Numerical variables
nrr=400; %number of points used for the radial integration
resample=4; % number of points skipped to obtain the different fields
lat_points=80; %lat_points: number of points used in the latitide,longitude grid
%% Non-dimentionalize
T_nd=2*pi/omega0;
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
[y, r]=tidal(l,R,rho,rhof,mu,Ks,etas,alpha,poro,k_perm,etaf,Kf,liquid,omega,'radial_points',nrr,'resample',resample);
k2_T=y(5,end,end)-1;

%% Obtain the different relevant fields, strain, stress, etc. 
%The different fields can be obtained as indicated in Appendix (D).
% Degree 2 order 0
[colat,lon,rr,displacements_20,flux_20,strain_20,stress_20,gravpot_20,p_fluid_20,C_fluid_20,Delta_poro_20,div_20, divq_20]=build_solution(y,r,R,rhof,rhos,Ks,Kf,mu,etas,etaf,k_perm,liquid,omega,alpha,poro,2,0,'lat_points',lat_points);
%degree 2, order 2
[colat,lon,rr,displacements_22c,flux_22c,strain_22c,stress_22c,gravpot_22c,p_fluid_22c,C_fluid_22c,Delta_poro_22c,div_22c,divq_22c]=build_solution(y,r,R,rhof,rhos,Ks,Kf,mu,etas,etaf,k_perm,liquid,omega,alpha,poro,2,2,'lat_points',lat_points);
%degree 2, order -2
[colat,lon,rr,displacements_22s,flux_22s,strain_22s,stress_22s,gravpot_22s,p_fluid_22s,C_fluid_22s,Delta_poro_22s,div_22s,divq_22s]=build_solution(y,r,R,rhof,rhos,Ks,Kf,mu,etas,etaf,k_perm,liquid,omega,alpha,poro,2,-2,'lat_points',lat_points);
%Add different components of the tidal potential Eq. (D6)
l=2;m=0;
norm20=sqrt((2*l+1)*factorial(l-m)/factorial(l+m)/(4*pi));
l=2;m=2;
norm22=sqrt(2*(2*l+1)*factorial(l-m)/factorial(l+m)/(4*pi));

strain_ND=3/2/norm20*strain_20-3/4/norm22*strain_22c+1i/norm22*strain_22s;
stress_ND=3/2/norm20*stress_20-3/4/norm22*stress_22c+1i/norm22*stress_22s;
p_fluid_ND=3/2/norm20*p_fluid_20-3/4/norm22*p_fluid_22c+1i/norm22*p_fluid_22s;
C_fluid_ND=3/2/norm20*C_fluid_20-3/4/norm22*C_fluid_22c+1i/norm22*C_fluid_22s;
Deltaporo_ND=3/2/norm20*Delta_poro_20-3/4/norm22*Delta_poro_22c+1i/norm22*Delta_poro_22s;
div_ND=3/2/norm20*div_20-3/4/norm22*div_22c+1i/norm22*div_22s;
divq_ND=3/2/norm20*divq_20-3/4/norm22*divq_22c+1i/norm22*divq_22s;
flux_ND=3/2/norm20*flux_20-3/4/norm22*flux_22c+1i/norm22*flux_22s;
displacements_ND=3/2/norm20*displacements_20-3/4/norm22*displacements_22c+1i/norm22*displacements_22s;
gravpot_ND=3/2/norm20*gravpot_20-3/4/norm22*gravpot_22c+1i/norm22*gravpot_22s;

% write in dimensional units
constant_f=(omega0*R0(end))^2*rho0(end)/mu0(end)*eccen; % Constant used to go from unit a tidal foricing to the actual magnitude of the tidal forcing

strain=constant_f*strain_ND;
stress=mu0(end)*constant_f*stress_ND;
p_fluid=mu0(end)*constant_f*p_fluid_ND;
flux=R0(end)/T_nd*constant_f*flux_ND;
displacements=R0(end)*constant_f*displacements_ND;
gravpot=mu0(end)/rho0(end)*constant_f*gravpot_ND;
C_fluid=constant_f*C_fluid_ND;
Deltaporo=constant_f*Deltaporo_ND;
div=constant_f*div_ND;
divq=constant_f*divq_ND;

%% Compute energy dissipation 
[energy_solid, energy_solid_pore,energy_fluid,energy_solid_surface, energy_solid_pore_surface,energy_fluid_surface,energy_solid_total, energy_fluid_total]=compute_energy(strain_ND,stress_ND,flux_ND,p_fluid_ND,C_fluid_ND,omega,etaf,k_perm,rr,colat,lon,poro);
E1=mu0(end)/T_nd*constant_f^2*R0(end)^3*energy_solid_total*1e-9;
E2=-21/2*(omega0*R0(end))^5/Gg0*eccen^2*imag(y(5,end,end))*1e-9;
disp(['Energy computed using direct Integraation'   num2str(E1, '%10.5e')  ' GW'])
% we can also compute it using the Love number in case there are not solid layers
disp(['Energy Computed using Love Number GW ' num2str(E2)])
disp(['Difference ' num2str((E1-E2)/(E1)*100) '%'])

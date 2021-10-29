%% FUNCTION: build_solution
% Description: Given the solution vector y(r), compute the solution in a colat-lon-r grid for a given order and degree.  
% Author: M. Rovira-Navarro 

%% INPUT 
%y: solution vector y(1:8,1:end,layer)
%r: radial points correspoding to y(1:8,:,layer)
%rhof: fluid density
%rhos: solid density
%Ks: bulk modulus of the solid 
%Kf: bulk modulus of the fluid
%mu: shear modulus of the solid
%etas: viscosity of the solid
%etaf: viscosity of the fluid
%liquid: 1 if layer is liquid
%k_perm: permeability 
%w_moon: forcing frequency of the moon
%alpha: Biot constant 
%poro: porosity
%l: degree of the tidal forcing
%m: order of the tidal forcing: 
    %0,-2 or 2.

%All physical parameters are given in non-dimensional units

% Optional input variables 
    % tidal_fluid: (1) if tides affect the fluid, (0) if not. 
    % lat_points: number of points used in the latitide,longitude grid. Default 70



%% OUTPUT 
% colat: colatitude where solution is given 
% lon: longitudes where solution is given 
% rr: list of radial points where solution is returned
% displacements: displcement vector  
                % displacements(icolat,ilon,ir,1): radial component of the displacement
                % displacements(icolat,ilon,ir,2): latitudinal component of the
                % displacement
                % displacements(icolat,ilon,ir,3): longitudinal component of the
                % diplacement
% flux: flux vector
                % flux(icolat,ilon,ir,1): radial component of the flux
                % flux(icolat,ilon,ir,2): latitudinal component of the flux
                % flux(icolat,ilon,ir,3): longitudinal component of the flux
% stress:
                %stress(icolat,ilon,ir,1)=\sigma_r_r;
                %stress(icolat,ilon,ir,2)=\sigma_theta_theta;
                %stress(icolat,ilon,ir,3)=\sigma_phi_phi;
                %stress(icolat,ilon,ir,4)=\sigma_r_theta;
                %stress(icolat,ilon,ir,5)=\sigma_r_phi;
                %stress(icolat,ilon,ir,6)=\sigma_theta_phi; 
% strain: strain tensor 
                %strain(icolat,ilon,ir,1)=\epsilon_r_r;
                %strain(icolat,ilon,ir,2)=\epsilon_theta_theta;
                %strain(icolat,ilon,ir,3)=\epsilon_phi_phi;
                %strain(icolat,ilon,ir,4)=\epsilon_r_theta;
                %strain(icolat,ilon,ir,5)=\epsilon_r_phi;
                %strain(icolat,ilon,ir,6)=\epsilon_theta_phi;                 
% gravpot: perturbing gravitational potential
                %gravpot(icolat,ilon,ir,1)
% p_fluid: fluid pore pressure
                %p_fluid(icolat,ilon,ir,1)
% C_fluid: variation of fluid content: 
                %C_fluid(icolat,ilon,ir,1)
% varargout can be used to get some extra outputs
                %varargout{1}: porosity change 
                
                
%% START OF THE FUNCTION
function [colat,lon,rr,displacements,flux,strain,stress,gravpot,p_fluid,C_fluid,varargout]=build_solution(y,r,R,rhof,rhos,Ks,Kf,mu,etas,etaf,k_perm,liquid,w_moon,alpha,poro,l,m,varargin)
%% READ OPTIONAL INPUTS
lat_points=70;
tidal_fluid=1;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'lat_points')
        lat_points=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'tidal_fluid')
        tidal_fluid=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end
tidal_on=tidal_fluid;
%% Compute some physical constants
Nlayers=length(R)-1;
nrr=length(r(1,:,1));
for i=1:length(R)-1
    % density 
    rho=(1-poro(i))*rhos(i)+poro(i)*rhof(i); 
    % drained bulk modulus
    if alpha(i)==1 %incompressible solid but matrix can still deform. (Cheng Section 4.8.1)
        Kd(i)=Ks(i);        %the matrix can rearrange and deform.  
        Ks(i)=Kd(i)*1e5;    %the solid matrix is incompressible 
        Ku(i)=Kd(i)+Kf(i)/poro(i);
    else
        Kd(i)=(1-alpha(i))*Ks(i);
        Ku(i)=Kd(i)+Kf(i)*Ks(i)*alpha(i)^2/(poro(i)*Ks(i)+(alpha(i)-poro(i))*Kf(i));
    end
    % undrained bulk modulus
    if poro(i)==0 | alpha(i)==0
        Ku(i)=0;
    end
    if liquid(i)==1
    muC(i)=1;
    lambda(i)=1;
    else
    muC(i)=mu(i)*1i*etas(i)*w_moon/(mu(i)+1i*etas(i)*w_moon);
    % lame parameter
    lambda(i)=Kd(i)-2/3*muC(i);
    end
end
%% OBTAIN SPHERICAL HARMONIC FUNCTIONS 
% spherical harmonic functions
colat=linspace(0,pi,lat_points);
lon=linspace(0,2*pi,lat_points);
Ylm=zeros(length(colat),length(lon));
Ylm_theta=zeros(length(colat),length(lon));
Ylm_theta_theta=zeros(length(colat),length(lon));
Ylm_phi=zeros(length(colat),length(lon));
Zlm=zeros(length(colat),length(lon));
Xlm=zeros(length(colat),length(lon));
if m==0
    norm=sqrt((2*l+1)*factorial(l-m)/factorial(l+m)/(4*pi));
    for i=1:length(colat)
        for j=1:length(lon)
            Ylm(i,j)=norm*(3*cos(colat(i))^2-1)/2*cos(m*lon(j));
            Ylm_theta(i,j)=-norm*3*cos(colat(i))*sin(colat(i))*cos(m*lon(j));
            Ylm_theta_theta(i,j)=norm*3*(sin(colat(i))^2-cos(colat(i))^2)*cos(m*lon(j));
            Ylm_phi(i,j)=0;
            Zlm(i,j)=0;
            Xlm(i,j)=2*Ylm_theta_theta(i,j)+l*(l+1)*Ylm(i,j);
        end
    end
elseif m==2
    norm=sqrt(2*(2*l+1)*factorial(l-m)/factorial(l+m)/(4*pi));
    for i=1:length(colat)
            for j=1:length(lon)
                Ylm(i,j)=norm*3*sin(colat(i))^2*cos(m*lon(j));
                Ylm_theta(i,j)=norm*6*cos(colat(i))*sin(colat(i))*cos(m*lon(j));
                Ylm_theta_theta(i,j)=norm*6*(cos(colat(i))^2-sin(colat(i))^2)*cos(m*lon(j));
                Ylm_phi(i,j)=-norm*6*sin(colat(i))*sin(m*lon(j));
                Zlm(i,j)=-norm*12*cos(colat(i))*sin(m*lon(j));
                Xlm(i,j)=2*Ylm_theta_theta(i,j)+l*(l+1)*Ylm(i,j);
            end
    end
else
    norm=sqrt(2*(2*l+1)*factorial(l-abs(m))/factorial(l+abs(m))/(4*pi));
    for i=1:length(colat)
        for j=1:length(lon)
            Ylm(i,j)=norm*3*sin(colat(i))^2*sin(abs(m)*lon(j));
            Ylm_theta(i,j)=norm*6*cos(colat(i))*sin(colat(i))*sin(abs(m)*lon(j));
            Ylm_theta_theta(i,j)=norm*6*(cos(colat(i))^2-sin(colat(i))^2)*sin(abs(m)*lon(j));
            Ylm_phi(i,j)=norm*6*sin(colat(i))*cos(abs(m)*lon(j));
            Zlm(i,j)=norm*12*cos(colat(i))*cos(abs(m)*lon(j));
            Xlm(i,j)=2*Ylm_theta_theta(i,j)+l*(l+1)*Ylm(i,j);
        end
    end    
end
    
%% BUILD THE SOLUTION PER LAYER
strain_lay=zeros(length(colat),length(lon),length(r(1,:,1)),6,Nlayers);
stress_lay=zeros(length(colat),length(lon),length(r(1,:,1)),6,Nlayers);
displacements_lay=zeros(length(colat),length(lon),length(r(1,:,1)),3,Nlayers);
flux_lay=zeros(length(colat),length(lon),length(r(1,:,1)),3,Nlayers);
C_fluid_lay=zeros(length(colat),length(lon),length(r(1,:,1)),1,Nlayers);
p_fluid_lay=zeros(length(colat),length(lon),length(r(1,:,1)),1,Nlayers);
gravpot_lay=zeros(length(colat),length(lon),length(r(1,:,1)),1,Nlayers);
Deltaporo_lay=zeros(length(colat),length(lon),length(r(1,:,1)),1,Nlayers);
dlm_lay=zeros(length(colat),length(lon),length(r(1,:,1)),1,Nlayers);
divq_lay=zeros(length(colat),length(lon),length(r(1,:,1)),1,Nlayers);
for ilay=1:Nlayers
    if liquid(ilay)==1 && ilay==1
    rr(1,1,:)=r(:,:,ilay);
    else
    ulm(1,1,:)=y(1,1:end,ilay);
    vlm(1,1,:)=y(2,1:end,ilay);
    y3(1,1,:)=y(3,1:end,ilay);
    y4(1,1,:)=y(4,1:end,ilay);
    y5(1,1,:)=y(5,1:end,ilay);
    y6(1,1,:)=y(6,1:end,ilay);
    y7(1,1,:)=y(7,1:end,ilay);
    y8(1,1,:)=y(8,1:end,ilay);
    rr(1,1,:)=r(:,:,ilay);
    ulm_r=y3/(lambda(ilay)+2*muC(ilay))+lambda(ilay)*l*(l+1)*vlm./((2*muC(ilay)+lambda(ilay))*rr)-2*lambda(ilay)*ulm./((2*muC(ilay)+lambda(ilay))*rr)+alpha(ilay)*y7/(lambda(ilay)+2*muC(ilay));
    vlm_r=y4/muC(ilay)+vlm./rr-ulm./rr;
    slm=vlm_r+(ulm-vlm)./rr;
    dlm=2*ulm./rr+ulm_r-l*(l+1)*vlm./rr;
    %compute dr(y8)
    prr1=1i*2*alpha(ilay)*w_moon*etaf(ilay)./(k_perm(ilay)*rr)*(1-lambda(ilay)/(2*mu(ilay)+lambda(ilay))).*ulm;
    prr2=1i*w_moon*alpha(ilay)*l*(l+1)*etaf(ilay)./(k_perm(ilay)*rr)*(-1+lambda(ilay)/(lambda(ilay)+2*mu(ilay))).*vlm;
    prr3=1i*w_moon*alpha(ilay)*etaf(ilay)./(k_perm(ilay)*(2*mu(ilay)+lambda(ilay))).*y3;
    prr5=tidal_fluid*l*(l+1)*rhof(ilay)./rr(ilay).^2.*y5;
    prr7=l*(l+1)./rr.^2.*y7+1i*w_moon*alpha(ilay)^2*etaf(ilay)/k_perm(ilay)*(1/(2*mu(ilay)+lambda(ilay))+1/(Ku(ilay)-Kd(ilay))).*y7;
    prr8=-2./rr.*y8;
    prr=prr1+prr2+prr3+prr5+prr7+prr8;
    divqml=-k_perm(ilay)/etaf(ilay)*(prr+2./rr.*y8-l*(l+1)./rr.^2.*(y7+rhof(ilay)*y5));    
    if alpha(ilay)==0
    Clm=zeros((size(dlm)));
    else
    Clm=alpha(ilay)*dlm+alpha(ilay)^2*y7/(Ku(ilay)-Kd(ilay));
    end
    %strain
    e_r_r=ulm_r.*Ylm;
    e_theta_theta=(ulm-l*(l+1)/2*vlm)./rr.*Ylm+vlm./rr.*Xlm(:,:)/2;
    e_phi_phi=(ulm-l*(l+1)/2*vlm)./rr.*Ylm-vlm./rr.*Xlm/2;
    e_r_theta=slm.*Ylm_theta/2;
    e_r_phi=slm.*Ylm_phi/2;
    e_theta_phi=vlm/2./rr.*Zlm;
    strain_lay(:,:,:,1,ilay)=e_r_r;
    strain_lay(:,:,:,2,ilay)=e_theta_theta;
    strain_lay(:,:,:,3,ilay)=e_phi_phi;
    strain_lay(:,:,:,4,ilay)=e_r_theta;
    strain_lay(:,:,:,5,ilay)=e_r_phi;
    strain_lay(:,:,:,6,ilay)=e_theta_phi;
    %fluid pressure
    p_fluid_lay(:,:,:,ilay)=y7.*Ylm;
    %varu=iation in fluid content
    C_fluid_lay(:,:,:,ilay)=Clm.*Ylm;
    %grav potential
    gravpot_lay(:,:,:,ilay)=y5.*Ylm;
    %flux and displacements
    if poro(ilay)>0
    flux_lay(:,:,:,1,ilay)=-k_perm(ilay)/(etaf(ilay))*y8.*Ylm;
    flux_lay(:,:,:,2,ilay)=-k_perm(ilay)/(etaf(ilay))*(y7+tidal_on*rhof(ilay)*y5)./rr.*Ylm_theta;
    flux_lay(:,:,:,3,ilay)=-k_perm(ilay)/(etaf(ilay))*(y7+tidal_on*rhof(ilay)*y5)./rr.*Ylm_phi;
    % compute the variable porosity
    Deltaporo_lay(:,:,:,1,ilay)=(1-poro(ilay))*(1/Kd(ilay)-1/((1-poro(ilay))*Ks(ilay)))*(Kd(ilay)*dlm-(alpha(ilay)-1)*y7).*Ylm;
    divq_lay(:,:,:,1,ilay)=divqml.*Ylm;
    end
    displacements_lay(:,:,:,1,ilay)=ulm.*Ylm;
    displacements_lay(:,:,:,2,ilay)=vlm.*Ylm_theta;
    displacements_lay(:,:,:,3,ilay)=vlm.*Ylm_phi;
    %divergence
    dlm_lay(:,:,:,1,ilay)=dlm.*Ylm;
    %stress
    s_r_r=lambda(ilay)*dlm.*Ylm+2*muC(ilay)*e_r_r-alpha(ilay)*p_fluid_lay(:,:,:,ilay);
    s_theta_theta=lambda(ilay)*dlm.*Ylm+2*muC(ilay)*e_theta_theta-alpha(ilay)*p_fluid_lay(:,:,:,ilay);
    s_phi_phi=lambda(ilay)*dlm.*Ylm+2*muC(ilay)*e_phi_phi-alpha(ilay)*p_fluid_lay(:,:,:,ilay);
    s_r_theta=2*muC(ilay)*e_r_theta;
    s_r_phi=2*muC(ilay)*e_r_phi;
    s_theta_phi=2*muC(ilay)*e_theta_phi;
    stress_lay(:,:,:,1,ilay)=s_r_r;
    stress_lay(:,:,:,2,ilay)=s_theta_theta;
    stress_lay(:,:,:,3,ilay)=s_phi_phi;
    stress_lay(:,:,:,4,ilay)=s_r_theta;
    stress_lay(:,:,:,5,ilay)=s_r_phi;
    stress_lay(:,:,:,6,ilay)=s_theta_phi;
    end
end
%% PUT THE SOLUTIONS TOGETHER
strain=zeros(length(colat),length(lon),length(r)*Nlayers,6);
stress=zeros(length(colat),length(lon),length(r)*Nlayers,6);
div=zeros(length(colat),length(lon),length(r)*Nlayers,1);
displacements=zeros(length(colat),length(lon),length(r)*Nlayers,3);
flux=zeros(length(colat),length(lon),length(r)*Nlayers,3);
C_fluid=zeros(length(colat),length(lon),length(r)*Nlayers);
p_fluid=zeros(length(colat),length(lon),length(r)*Nlayers);
gravpot=zeros(length(colat),length(lon),length(r)*Nlayers);
Deltaporo=zeros(length(colat),length(lon),length(r)*Nlayers);
divq=zeros(length(colat),length(lon),length(r)*Nlayers);
rr=zeros(1,length(r)*Nlayers);
i_start=1;
%this can be quite slow...
if Nlayers==1
    rr=r;
    strain=strain_lay;
    stress=stress_lay;
    displacements=displacements_lay;
    flux=flux_lay;
    C_fluid= C_fluid_lay;
    div= dlm_lay;
    divq= divq_lay;
    p_fluid= p_fluid_lay;
    gravpot= gravpot_lay;
    Deltaporo= Deltaporo_lay;
else
    for ilay=1:Nlayers
        rr(1,i_start:i_start+nrr-1)=r(1,:,ilay);
        strain(:,:,i_start:i_start+nrr-1,1)=strain_lay(:,:,:,1,ilay);
        strain(:,:,i_start:i_start+nrr-1,2)=strain_lay(:,:,:,2,ilay);
        strain(:,:,i_start:i_start+nrr-1,3)=strain_lay(:,:,:,3,ilay);
        strain(:,:,i_start:i_start+nrr-1,4)=strain_lay(:,:,:,4,ilay);
        strain(:,:,i_start:i_start+nrr-1,5)=strain_lay(:,:,:,5,ilay);
        strain(:,:,i_start:i_start+nrr-1,6)=strain_lay(:,:,:,6,ilay);
        stress(:,:,i_start:i_start+nrr-1,1)=stress_lay(:,:,:,1,ilay);
        stress(:,:,i_start:i_start+nrr-1,2)=stress_lay(:,:,:,2,ilay);
        stress(:,:,i_start:i_start+nrr-1,3)=stress_lay(:,:,:,3,ilay);
        stress(:,:,i_start:i_start+nrr-1,4)=stress_lay(:,:,:,4,ilay);
        stress(:,:,i_start:i_start+nrr-1,5)=stress_lay(:,:,:,5,ilay);
        stress(:,:,i_start:i_start+nrr-1,6)=stress_lay(:,:,:,6,ilay);
        displacements(:,:,i_start:i_start+nrr-1,1)=displacements_lay(:,:,:,1,ilay);
        displacements(:,:,i_start:i_start+nrr-1,2)=displacements_lay(:,:,:,2,ilay);
        displacements(:,:,i_start:i_start+nrr-1,3)=displacements_lay(:,:,:,3,ilay);
        flux(:,:,i_start:i_start+nrr-1,1)=flux_lay(:,:,:,1,ilay);
        flux(:,:,i_start:i_start+nrr-1,2)=flux_lay(:,:,:,2,ilay);
        flux(:,:,i_start:i_start+nrr-1,3)=flux_lay(:,:,:,3,ilay);
        C_fluid(:,:,i_start:i_start+nrr-1)= C_fluid_lay(:,:,:,ilay);
        div(:,:,i_start:i_start+nrr-1)= dlm_lay(:,:,:,ilay);
        divq(:,:,i_start:i_start+nrr-1)= divq_lay(:,:,:,ilay);
        p_fluid(:,:,i_start:i_start+nrr-1)= p_fluid_lay(:,:,:,ilay);
        gravpot(:,:,i_start:i_start+nrr-1)= gravpot_lay(:,:,:,ilay);
        Deltaporo(:,:,i_start:i_start+nrr-1)= Deltaporo_lay(:,:,:,ilay);
        i_start=i_start+nrr;
    end
end
varargout{1}=Deltaporo;
varargout{2}=div;
varargout{3}=divq;
end
%% FUNCTION: tidal
% Description: Propagate the solution from the center to the surface  
% Author: M. Rovira-Navarro 
%% INPUT 
%l: spherical harmonic degree
%R: vector containing the upper boundary of each layer
%rho: vector containing the average density of each layer
%rhof: vector containing the fluid average density of each layer
%mu: vector containing the shear modulus of each layer
%Ks: vector containing the solid bulk modulus of each layer
%etas: vector containing the viscosity of the solid 
%alpha: vector conatining alpha
%poro: vector containing porosity for each layer 
%k_perm: vector containing permeability for each layer
%etaf: vector containing fluid viscosity for each layer
%Kf: vector containing the fluid bulk modulus of each layer
%liquid: vector stating if layer is fluid (1) or not (0). Here used to indicate the model has a fluid core
%omega: forcing frequency


% some optional input in varargin: 
    %radial_points: number of radial points
    %pressure_BC: use pressure boundary condition instead of constant flux
    %strain_BC: use strain
    %print_results: print some results on the screen solution at the surface, love numbers and y plots
    %self_gravity: (1) if self-gravity is used in the momentum equation, (0) if not
    %tidal_fluid: (1) if tidal potential affects the fluid, (0) if not
    %gravity_on: turn gravity on (1) and off (0). Default 1 

%% OUTPUT 
%y(1:8,1:nrr,1:Nlayers): solution vector for each layer
    % y1: normal displacement 
    % y2: tangential displacements
    % y3: normal stress
    % y4: tangential stress
    % y5: perturbing potential 
    % y6: potential stress
    % y7: pore pressure
    % y8: radial flux
%r: radial points where solution is obtained 
%% FUNCTION 
function [y,r,varargout]=tidal(l,R,rho,rhof,mu,Ks,etas,alpha,poro,k_perm,etaf,Kf,liquid,omega, varargin)
global Gg
%% SOME PHYSICAL AND NUMERICAL CONSTANTS
% PHYSICAL PARAMETERS
Rs=R(end); %surface radius 
glayer(1)=0;
Nlayers=length(R)-1;
layporo=NaN;
for i=1:length(R)-1
    % drained bulk modulus 
    if alpha(i)==1 %incompressible solid but matrix can still deform. (Cheng Section 4.8.1)
        Kd(i)=Ks(i);        %the matrix can rearrange and deform, the bulk modulus is already given by the user.   
        Ks(i)=Kd(i)*1e5;    %the solid matrix is incompressible 
        Ku(i)=Kd(i)+Kf(i)/poro(i);
    else
        if poro(i)==0 | alpha(i)==0
            Ku(i)=0;
            Kd(i)=(1-alpha(i))*Ks(i);
        else
            Kd(i)=(1-alpha(i))*Ks(i);
            Ku(i)=Kd(i)+Kf(i)*Ks(i)*alpha(i)^2/(poro(i)*Ks(i)+(alpha(i)-poro(i))*Kf(i));
        end
    end
    % Andrade
%     alphaA=0.3;
%     gammaA=0.8975;
%     betaA=mu(i)^(alphaA-1)*etas(i)^(-alphaA);
%     JJA=1/mu(i)+1/(etas(i)*omega*1i)+betaA*(omega*1i)^(-alphaA)*gammaA;
%     JJA=1/JJA;
%     muC(i)=JJA;     
    % Maxwell
    muC(i)=mu(i)*1i*etas(i)*omega/(mu(i)+1i*etas(i)*omega);
    % lame parameter
    lambda(i)=Kd(i)-2/3*muC(i);
    glayer(i+1)=glayer(i)*(R(i)/R(i+1))^2+4/3*pi*Gg*rho(i)*(R(i+1)^3-R(i)^3)/R(i+1)^2;
end

%% NUMERICAL CONSTANTS + USER FLAGS 
pressure=NaN;
nrr=1000; %number of points in the radial direction
print_results=0;
strain_BC=0;
shear_BC=0;
displacement_BC=0;
resample=0;
self_gravity=1;
tidal_fluid=1;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'radial_points')
        nrr=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'pressure_BC')
        pressure=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
     end
     if strcmpi(varargin{k},'strain_BC')
        strain_BC=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
     end
     if strcmpi(varargin{k},'shear_BC')
        shear_BC=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
     end
     if strcmpi(varargin{k},'self_gravity')
        self_gravity=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'tidal_fluid')
        tidal_fluid=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
     if strcmpi(varargin{k},'displacement_BC')
        displacement_BC=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
     end
     if strcmpi(varargin{k},'gravity_on')
        gravity_on=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
        rho=rho*gravity_on;
        rhof=rho*gravity_on;
    end
    if strcmpi(varargin{k},'print_results')
        print_results=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'resample')
        resample=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end
% constants used in the RK4 integrator
A2=1/5; A3=3/10; A4=3/5; A5=1; A6=7/8;
B21=1/5;
B31=3/40; B32=9/40;
B41=3/10; B42=-9/10; B43=6/5;
B51=-11/54; B52=5/2; B53=-70/27; B54=35/27;
B61=1631/55296; B62=175/512; B63=575/13824; B64=44275/110592; B65=253/4096; 
AC1=37/378; AC2=0; AC3=250/621; AC4=125/594; AC5=0; AC6=512/1771; 
DC1=2825/27648; DC2=0; DC3=18575/48384; DC4=13525/55296; DC5=277/14336; DC6=1/4;
%% INITIALIZE SOLUTION  
%%%%%%%%%%%%%%% 
% Initialize vectors
y1=zeros(8,nrr+1,Nlayers); %solution 1
y2=zeros(8,nrr+1,Nlayers); %solution 2
y3=zeros(8,nrr+1,Nlayers); %solution 3
y4=zeros(8,nrr+1,Nlayers); %solution 4

r=zeros(1,nrr+1,Nlayers);  %r
g=zeros(1,nrr+1,Nlayers); %g

%% BOUNDARY CONDITIONS  AT THE CENTER/CORE
%%%%%%%%%%%%%%% 
%vector of integration constnats
C1=zeros(4,1);
C1(1,1)=1;
C2=zeros(4,1);
C2(2,1)=1;
C3=zeros(4,1);
C3(3,1)=1;
C4=zeros(4,1);
C4(4,1)=0;
if liquid(1)==0 %no liquid core
    Icore=zeros(8,4);
    %normal stress
    Icore(3,1)=1;
    % tangential stress
    Icore(4,2)=1;
    % potential stress
    Icore(6,3)=1;
    % pressure
    Icore(7,4)=1;
    %at which layer does integration start 
    layer_start=1;
    %lag in the integration
    R(1)=R(2)*1e-3;
    r(1,1,layer_start)=R(1);
    g(1,1,layer_start)=0;
else %liquid core
    Ac  = 4/3*pi*Gg*rho(1);
    Icore=zeros(8,4);
    Icore(1,1)= -R(2)^(l-1)/Ac;
    Icore(5,1)= R(2)^l;
    Icore(6,1)= 2*(l-1)*R(2)^(l-1);
    Icore(2,2)= 1;
    Icore(1,3)= 1;
    Icore(3,3)= rho(1)*Ac*R(2);
    Icore(6,3)= 3*Ac;
    Icore(7,4)=1;
    layer_start=2;
    %starting of the integration;
    r(1,1,layer_start)=R(2);
    r(1,1:nrr+1,1)=linspace(R(1),R(2),nrr+1);
    g(1,1,layer_start)=glayer(2);
end

y1(:,1,layer_start)=Icore*C1;
y2(:,1,layer_start)=Icore*C2;
y3(:,1,layer_start)=Icore*C3;
y4(:,1,layer_start)=Icore*C4;



%% PROPAGATE SOLUTION THROUGH THE DIFFERENT LAYERS 
%%%%%%%%%%%%%%% 
for ilay=layer_start:Nlayers 
    % Starting vector
    if ilay>layer_start
        y1(1:6,1,ilay)=y1(1:6,ifin(ilay-1),ilay-1);
        y2(1:6,1,ilay)=y2(1:6,ifin(ilay-1),ilay-1);
        y3(1:6,1,ilay)=y3(1:6,ifin(ilay-1),ilay-1);
        y4(1:6,1,ilay)=y4(1:6,ifin(ilay-1),ilay-1);
        r(1,1,ilay)=r(1,ifin(ilay-1),ilay-1);
        g(1,1,ilay)=g(1,ifin(ilay-1),ilay-1);
    end
    % If the layer is porous introuduce new integration constant C4
    if poro(ilay)>0 
        %y4(:,1,ilay)=[0;0;0;0;0;0;1;0]; if radial velocity is 0 
        y4(:,1,ilay)=[0;0;0;0;0;0;1;0];
        layporo=ilay;
    end
    % Delta_r
    Delta_r=(R(ilay+1)-R(ilay))/nrr;
    i=1;
% propagate inside the layer
        while (r(1,i,ilay)<R(ilay+1)) 
        if (abs(r(1,i,ilay)+Delta_r-R(ilay+1))/R(ilay+1)<1e-8) 
                  Delta_r=R(ilay+1)-r(1,i,ilay);
        end
        r(1,i+1,ilay)=r(1,i,ilay)+Delta_r;
        g(1,i+1,ilay)=g(1,1,ilay)*(R(ilay)/r(1,i+1,ilay))^2+4/3*pi*Gg*rho(ilay)*(r(1,i+1,ilay)^3-R(ilay)^3)/r(1,i+1,ilay)^2;
        %glayer(i+1)=glayer(i)*(R(i)/R(i+1))^2+4/3*pi*Gg*rho(i)*(R(i+1)^3-R(i)^3)/R(i+1)^2
        rK=r(1,i,ilay);
        gK=g(1,i,ilay);
        %%%%%%%%% Step 1
        Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
        k11=Delta_r*(Ymatrix*y1(1:8,i,ilay));
        k21=Delta_r*(Ymatrix*y2(1:8,i,ilay));
        k31=Delta_r*(Ymatrix*y3(1:8,i,ilay));
        k41=Delta_r*(Ymatrix*y4(1:8,i,ilay));
        %%%%%%%%% Step 2
        rK=r(1,i,ilay)+A2*Delta_r;
        gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
        Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
        k12=Delta_r*Ymatrix*(y1(1:8,i,ilay)+B21*k11);
        k22=Delta_r*Ymatrix*(y2(1:8,i,ilay)+B21*k21);
        k32=Delta_r*Ymatrix*(y3(1:8,i,ilay)+B21*k31);
        k42=Delta_r*Ymatrix*(y4(1:8,i,ilay)+B21*k41);
        %%%%%%%%% Step 3
        rK=r(1,i,ilay)+A3*Delta_r;
        gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
        Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
        k13=Delta_r*Ymatrix*(y1(1:8,i,ilay)+B31*k11+B32*k12);
        k23=Delta_r*Ymatrix*(y2(1:8,i,ilay)+B31*k21+B32*k22);
        k33=Delta_r*Ymatrix*(y3(1:8,i,ilay)+B31*k31+B32*k32);
        k43=Delta_r*Ymatrix*(y4(1:8,i,ilay)+B31*k41+B32*k42);
        %%%%%%%%% Step 4
        rK=r(1,i,ilay)+A4*Delta_r;
        gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
        Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
        k14=Delta_r*Ymatrix*(y1(1:8,i,ilay)+B41*k11+B42*k12+B43*k13);
        k24=Delta_r*Ymatrix*(y2(1:8,i,ilay)+B41*k21+B42*k22+B43*k23);
        k34=Delta_r*Ymatrix*(y3(1:8,i,ilay)+B41*k31+B42*k32+B43*k33);
        k44=Delta_r*Ymatrix*(y4(1:8,i,ilay)+B41*k41+B42*k42+B43*k43);
        %%%%%%%%% Step 5
        rK=r(1,i,ilay)+A5*Delta_r;
        gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
        Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
        k15=Delta_r*Ymatrix*(y1(1:8,i,ilay)+B51*k11+B52*k12+B53*k13+B54*k14);
        k25=Delta_r*Ymatrix*(y2(1:8,i,ilay)+B51*k21+B52*k22+B53*k23+B54*k24);
        k35=Delta_r*Ymatrix*(y3(1:8,i,ilay)+B51*k31+B52*k32+B53*k33+B54*k34);
        k45=Delta_r*Ymatrix*(y4(1:8,i,ilay)+B51*k41+B52*k42+B53*k43+B54*k44);
        %%%%%%%%% Step 6
        rK=r(1,i,ilay)+A6*Delta_r;
        gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
        Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
        k16=Delta_r*Ymatrix*(y1(1:8,i,ilay)+B61*k11+B62*k12+B63*k13+B64*k14+B65*k15);
        k26=Delta_r*Ymatrix*(y2(1:8,i,ilay)+B61*k21+B62*k22+B63*k23+B64*k24+B65*k25);
        k36=Delta_r*Ymatrix*(y3(1:8,i,ilay)+B61*k31+B62*k32+B63*k33+B64*k34+B65*k35);
        k46=Delta_r*Ymatrix*(y4(1:8,i,ilay)+B61*k41+B62*k42+B63*k43+B64*k44+B65*k45);
        %%%%%%%%%% Get next step
        y1(1:8,i+1,ilay)=y1(1:8,i,ilay)+AC1*k11+AC2*k12+AC3*k13+AC4*k14+AC5*k15+AC6*k16;
        y2(1:8,i+1,ilay)=y2(1:8,i,ilay)+AC1*k21+AC2*k22+AC3*k23+AC4*k24+AC5*k25+AC6*k26;
        y3(1:8,i+1,ilay)=y3(1:8,i,ilay)+AC1*k31+AC2*k32+AC3*k33+AC4*k34+AC5*k35+AC6*k36;
        y4(1:8,i+1,ilay)=y4(1:8,i,ilay)+AC1*k41+AC2*k42+AC3*k43+AC4*k44+AC5*k45+AC6*k46;  
        i=i+1;    
        end
        ifin(ilay)=i;
end

%% BOUNDARY CONDITIONS
%%%%%%%%%%%%%%% 
% boundary conditions at the surface    
% row 2: tangent stress 
M(2,1)=y1(4,ifin(end),Nlayers);
M(2,2)=y2(4,ifin(end),Nlayers);
M(2,3)=y3(4,ifin(end),Nlayers);
M(2,4)=y4(4,ifin(end),Nlayers);

B=zeros(4,1);
if strain_BC~=0  %stain boundary condition
    %impose strain 
    B(1)=strain_BC;
    M(1,1)=y1(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y1(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y1(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y1(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    M(1,2)=y2(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y2(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y2(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y2(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    M(1,3)=y3(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y3(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y3(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y3(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    M(1,4)=y4(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y4(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y4(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y4(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    %impose displacemenent
%      B(1)=strain_BC/gK;
%      M(1,1)=y1(1,ifin);
%      M(1,2)=y2(1,ifin);
%      M(1,3)=y3(1,ifin);
%      M(1,4)=y4(1,ifin);
     % impose non-rotational displacement field
     if shear_BC==0
        M(2,1)=2*(y1(1,ifin(end),Nlayers)-y1(2,ifin(end),Nlayers))-y1(4,ifin(end),Nlayers)/(muC(Nlayers));
        M(2,2)=2*(y2(1,ifin(end),Nlayers)-y2(2,ifin(end),Nlayers))-y2(4,ifin(end),Nlayers)/(muC(Nlayers));
        M(2,3)=2*(y3(1,ifin(end),Nlayers)-y3(2,ifin(end),Nlayers))-y3(4,ifin(end),Nlayers)/(muC(Nlayers));
        M(2,4)=2*(y4(1,ifin(end),Nlayers)-y4(2,ifin(end),Nlayers))-y4(4,ifin(end),Nlayers)/(muC(Nlayers)); 
     else
        M(2,1)=y1(4,ifin(end),Nlayers);
        M(2,2)=y2(4,ifin(end),Nlayers);
        M(2,3)=y3(4,ifin(end),Nlayers);
        M(2,4)=y4(4,ifin(end),Nlayers);
     end
%      M(2,1)=2*(y1(1,ifin(end),Nlayers)-y1(2,ifin(end),Nlayers))-y1(4,ifin(end),Nlayers)/(muC(Nlayers));
%      M(2,2)=2*(y2(1,ifin(end),Nlayers)-y2(2,ifin(end),Nlayers))-y2(4,ifin(end),Nlayers)/(muC(Nlayers));
%      M(2,3)=2*(y3(1,ifin(end),Nlayers)-y3(2,ifin(end),Nlayers))-y3(4,ifin(end),Nlayers)/(muC(Nlayers));
%      M(2,4)=2*(y4(1,ifin(end),Nlayers)-y4(2,ifin(end),Nlayers))-y4(4,ifin(end),Nlayers)/(muC(Nlayers));     
     B(4)=0;
     M(4,1)=y1(6,ifin(end),Nlayers);
     M(4,2)=y2(6,ifin(end),Nlayers);
     M(4,3)=y3(6,ifin(end),Nlayers);
     M(4,4)=y4(6,ifin(end),Nlayers);     
elseif displacement_BC~=0  %displacement boundary condition 
    B(1)=displacement_BC;
    M(1,1)=y1(1,ifin(end),Nlayers);
    M(1,2)=y2(1,ifin(end),Nlayers);
    M(1,3)=y3(1,ifin(end),Nlayers);
    M(1,4)=y4(1,ifin(end),Nlayers);
    
    % potential stress 
    B(4)=0;
    M(4,1)=y1(6,ifin(end),Nlayers);
    M(4,2)=y2(6,ifin(end),Nlayers);
    M(4,3)=y3(6,ifin(end),Nlayers);
    M(4,4)=y4(6,ifin(end),Nlayers);
else %tidally forced problem     
    % row 1: Normal stress is 0 
    B(1)=0;
    M(1,1)=y1(3,ifin(end),Nlayers);
    M(1,2)=y2(3,ifin(end),Nlayers);
    M(1,3)=y3(3,ifin(end),Nlayers);
    M(1,4)=y4(3,ifin(end),Nlayers);
    % potential stress 
    B(4)=(2*l+1)/Rs;
    M(4,1)=y1(6,ifin(end),Nlayers);
    M(4,2)=y2(6,ifin(end),Nlayers);
    M(4,3)=y3(6,ifin(end),Nlayers);
    M(4,4)=y4(6,ifin(end),Nlayers);
end
% BOUNDARY CONDITIONS AT THE POROUS BOUNDARY

if (isnan(layporo)==1) %No porus layers
    M(3,3)=0;
    M(3,1)=0;
    M(3,2)=0;
    M(3,4)=1; 
else
    if isnan(pressure)==1  %no flux boundary condition
    B(3)=0;
    M(3,1)=y1(8,ifin(layporo),layporo);
    M(3,2)=y2(8,ifin(layporo),layporo);
    M(3,3)=y3(8,ifin(layporo),layporo);
    M(3,4)=y4(8,ifin(layporo),layporo);
    else  %\sigma_rr+p=0
    B(3)=pressure;
    M(3,1)=y1(7,ifin(layporo),layporo)+y1(3,ifin(layporo),layporo);
    M(3,2)=y2(7,ifin(layporo),layporo)+y2(3,ifin(layporo),layporo);
    M(3,3)=y3(7,ifin(layporo),layporo)+y3(3,ifin(layporo),layporo);
    M(3,4)=y4(7,ifin(layporo),layporo)+y4(3,ifin(layporo),layporo);
    end
end

%% INVERT MATRIX 
sol=M\B;
%sol=double(inv(sym(M))*B);
y=sol(1)*y1+sol(2)*y2+sol(3)*y3+sol(4)*y4;

%resample for other more expensive computations
if resample>0
        yres=y(:,1:resample:end,:);
        y=yres;
        rres=r(:,1:resample:end,:);
        r=rres;
end
% sol=M(1:3,1:3)\B(1:3);
% y=sol(1)*y1(1:8,1:ifin)+sol(2)*y2(1:8,1:ifin)+sol(3)*y3(1:8,1:ifin);

%% RETURN SOME OPTIONAL OUTPUTS
varargout{1}=y(3,:,Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y(2,:,Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y(1,:,Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y(7,:,Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));

%% PRINT SOME USEFUL INFORMATION 
if print_results==1
    disp("====================================")
    cprintf('*Blue' , 'Solution at the surface \n')
    for i=1:8
        disp(['y' num2str(i) ' '   num2str(y(i,ifin(Nlayers),Nlayers))])
    end
    cprintf('*Blue' , 'Love numbers \n')
    disp(['h2 '   num2str(g(1,ifin(Nlayers),Nlayers)*y(1,ifin(Nlayers),Nlayers), '%10.5e')])
    disp(['l2 '   num2str(g(1,ifin(Nlayers),Nlayers)*y(2,ifin(Nlayers),Nlayers), '%10.5e')])
    disp(['k2 '   num2str(-1-y(5,ifin(Nlayers),Nlayers), '%10.5e')])

    %% PLOTS FOR DEBUG
    figure
    title("y1")
    for i=1:8
        subplot(2,4,i)
        for j=layer_start:Nlayers
        plot(r(1,1:ifin(j),j),abs(y1(i,1:ifin(j),j)));
        hold on
        end
        set(gca, 'YScale', 'log')
    end
    figure
    title("y2")
    for i=1:8
    subplot(2,4,i)
        for j=layer_start:Nlayers
        plot(r(1,1:ifin(j),j),abs(y2(i,1:ifin(j),j)));
        hold on
        end
    set(gca, 'YScale', 'log')
    end
    figure
    for i=1:8
    subplot(2,4,i)
        for j=layer_start:Nlayers
        plot(r(1,1:ifin(j),j),abs(y3(i,1:ifin(j),j)));
        hold on 
        end
    set(gca, 'YScale', 'log')
    end
    figure
    title("y4")
    for i=1:8
    subplot(2,4,i)
        for j=layer_start:Nlayers
        plot(r(1,1:ifin(j),j),abs(y4(i,1:ifin(j),j)));
        hold on
        end
    set(gca, 'YScale', 'log')
    end

    figure
    title("y")
    for i=1:8
    subplot(2,4,i)
        for j=layer_start:Nlayers
        plot(r(1,1:ifin(j),j),real(y(i,1:ifin(j),j)),'LineWidth',2);
        hold on
        end
    end
    %plot the curl of the displacementes 
    figure
    for j=layer_start:Nlayers
    maxu=abs((y(1,1:ifin(j),j)));
    maxu=max(maxu(:));
    rot=2./r(1,1:ifin(j),j).*((y(1,1:ifin(j),j))-(y(2,1:ifin(j),j)))-(y(4,1:ifin(j),j))/muC(j);
    plot(r(1,1:ifin(j),j),abs(rot)/maxu,'LineWidth',2);
    hold on
    end
    xlabel('r')
    ylabel('$\nabla\times u/max(u)$')
    xlim([0.1 1])
    
    a_poro=find(poro>0);

    if isempty(a_poro)==0
    figure
    subplot(1,3,1)
    for j=layer_start:Nlayers   
    plot(r(1,1:ifin(j),j),real(y(7,1:ifin(j),j)+tidal_fluid*rhof(j)*y(5,1:ifin(j),j)),'LineWidth',2,'color','k');
    hold on
    plot(r(1,1:ifin(j),j),imag(y(7,1:ifin(j),j)+tidal_fluid*rhof(j)*y(5,1:ifin(j),j)),'LineWidth',2,'color','k','LineStyle','--');
    hold on
    end
    xlabel('r')
    ylabel('Driving force Darcy flow')
    a_poro=find(poro>0);
    xlim([r(1,1,a_poro) r(1,end,a_poro)])
   
    subplot(1,3,2)
    for j=layer_start:Nlayers   
    plot(r(1,1:ifin(j),j),real(y(8,1:ifin(j),j)),'LineWidth',2,'color','k');
    hold on
    plot(r(1,1:ifin(j),j),imag(y(8,1:ifin(j),j)),'LineWidth',2,'color','k','LineStyle','--');
    hold on
    end
    xlabel('r')
    ylabel('Radial flow')
    a_poro=find(poro>0);
    xlim([r(1,1,a_poro) r(1,end,a_poro)])
    
    subplot(1,3,3)
    for j=layer_start:Nlayers   
    plot(r(1,1:ifin(j),j),real(y(7,1:ifin(j),j)),'LineWidth',2,'color','k');
    hold on
    plot(r(1,1:ifin(j),j),imag(y(7,1:ifin(j),j)),'LineWidth',2,'color','k','LineStyle','--');
    hold on
    plot(r(1,1:ifin(j),j),real(rhof(j)*y(5,1:ifin(j),j)),'LineWidth',2,'color','b');
    hold on
    plot(r(1,1:ifin(j),j),imag(rhof(j)*y(5,1:ifin(j),j)),'LineWidth',2,'color','b','LineStyle','--');
    end
    xlabel('r')
    ylabel('pressure/potential')
    a_poro=find(poro>0);
    xlim([r(1,1,a_poro) r(1,end,a_poro)])
    end
    
end


end
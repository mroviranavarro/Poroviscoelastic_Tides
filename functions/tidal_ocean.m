%% FUNCTION: tidal_ocean
% Description: Propagate the solution from the center to the surface but with the possibility of a subsurface ocean 
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
%liquid: vector stating of layer is fluid (1) or not (0)
%omega: forcing frequency


%some optional input in varargin: 
    %radial_points: number of radial points
    %pressure_BC: use pressure boundary condition instead of constant flux
    %strain_BC: use strain
    %print_results: print some results on the screen, solution at the surface, love numbers and y plots
    % resample: To solve the PDE it is good to have  a lot of points in the  radial direction, but then to build the solution u need to deal with big matrices that slows things down. In such case the resample option is useful.
    %self_gravity: (1) if self-gravity is used in the momentum equation, (0) if not
    %tidal_fluid: (1) if tidal potential affects the fluid, (0) if not
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
function [y,r]=tidal_ocean(l,R,rho,rhof,mu,Ks,etas,alpha,poro,k_perm,etaf,Kf,liquid,omega,compe,tidal_fluid, varargin)
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
    muC(i)=mu(i)*1i*etas(i)*omega/(mu(i)+1i*etas(i)*omega);
    lambda(i)=Kd(i)-2/3*muC(i);
    glayer(i+1)=glayer(i)*(R(i)/R(i+1))^2+4/3*pi*Gg*rho(i)*(R(i+1)^3-R(i)^3)/R(i+1)^2;

end

%% NUMERICAL CONSTANTS + USER FLAGS 
pressure=NaN;
nrr=1000; %number of points in the radial direction
print_results=0;
strain_BC=0;
self_gravity=1;
tidal_fluid=1;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'radial_points')
        nrr=varargin{k+1}; 
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

%% FIND OCEAN LAYER
layer_ocean=find(liquid(2:end)>0)+1;
%% INITIALIZE SOLUTION  
%Before the ocean
y1=zeros(8,nrr+1,Nlayers); %solution 1
y2=zeros(8,nrr+1,Nlayers); %solution 2
y3=zeros(8,nrr+1,Nlayers); %solution 3
y4=zeros(8,nrr+1,Nlayers); %solution 4
%In the ocean 
y5=zeros(8,nrr+1,Nlayers); %solution 5, also z5 (tidal potential)
y6=zeros(8,nrr+1,Nlayers); %solution 6, also z6 (potential stress)
%After the ocean
y7=zeros(8,nrr+1,Nlayers); %solution 7
y8=zeros(8,nrr+1,Nlayers); %solution 8
y9=zeros(8,nrr+1,Nlayers); %solution 9
y10=zeros(8,nrr+1,Nlayers); %solution 10

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
    R(1)=R(2)*5e-4;
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
    %starting of the integration
    r(1,1,layer_start)=R(2);
    r(1,1:nrr+1,1)=linspace(R(1),R(2),nrr+1);
    g(1,1,layer_start)=glayer(2);
end

y1(:,1,layer_start)=Icore*C1;
y2(:,1,layer_start)=Icore*C2;
y3(:,1,layer_start)=Icore*C3;
y4(:,1,layer_start)=Icore*C4;

%% OCEAN-MANTLE INTERFACE BOUNDARY CONDITIONS
C5=zeros(2,1);
C5(1,1)=1;
C6=zeros(2,1);
C6(2,1)=1;
Ioceanmantle=zeros(8,2);
Ioceanmantle(5,1)=1;
Ioceanmantle(6,2)=1;

y5(:,1,layer_ocean)=Ioceanmantle*C5;
y6(:,1,layer_ocean)=Ioceanmantle*C6;

%% OCEAN-ICE INTERFACE BOUNDARY CONDITIONS
C7=zeros(4,1);
C7(1,1)=1;
C8=zeros(4,1);
C8(2,1)=1;
C9=zeros(4,1);
C9(3,1)=1;
C10=zeros(4,1);
C10(4,1)=1;

g_top=glayer(layer_ocean+1);
rho_ocean=rho(layer_ocean);

Ioceanice=zeros(8,4);
Ioceanice(1,1)=1;
Ioceanice(1,3)=-1/g_top;
Ioceanice(2,2)=1;
Ioceanice(3,1)=rho_ocean*g_top;
Ioceanice(5,3)=1;
Ioceanice(6,4)=1;

y7(:,1,layer_ocean+1)=Ioceanice*C7;
y8(:,1,layer_ocean+1)=Ioceanice*C8;
y9(:,1,layer_ocean+1)=Ioceanice*C9;
y10(:,1,layer_ocean+1)=Ioceanice*C10;

%% PROPAGATE SOLUTION THROUGH THE DIFFERENT LAYERS  
%%%%%%%%%%%%%%% 
for ilay=layer_start:Nlayers
    if ilay>layer_start
        g(1,1,ilay)=g(1,ifin(ilay-1),ilay-1);
        r(1,1,ilay)=r(1,ifin(ilay-1),ilay-1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LAYERS BEFORE THE OCEAN
    if  ilay<layer_ocean
        % INITIALIZE
        if ilay>layer_start
        y1(1:6,1,ilay)=y1(1:6,ifin(ilay-1),ilay-1);
        y2(1:6,1,ilay)=y2(1:6,ifin(ilay-1),ilay-1);
        y3(1:6,1,ilay)=y3(1:6,ifin(ilay-1),ilay-1);
        y4(1:6,1,ilay)=y4(1:6,ifin(ilay-1),ilay-1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE POROUS LAYER GOES HERE 
        if poro(ilay)>0 
            y4(:,1,ilay)=[0;0;0;0;0;0;1;0];
            layporo=ilay;
        end
        Delta_r=(R(ilay+1)-R(ilay))/nrr;
        i=1;
        %INTEGRATE THE SOLUTION UPWARDS 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCEAN STARTS HERE
    if ilay==layer_ocean
        Delta_r=(R(ilay+1)-R(ilay))/nrr;
        i=1;
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
            Ymatrix=propagation_matrix_ocean(l,rK,rho(ilay),gK);
            k11=Delta_r*(Ymatrix*y5(1:8,i,ilay));
            k21=Delta_r*(Ymatrix*y6(1:8,i,ilay));
            %%%%%%%%% Step 2
            rK=r(1,i,ilay)+A2*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_ocean(l,rK,rho(ilay),gK);
            k12=Delta_r*Ymatrix*(y5(1:8,i,ilay)+B21*k11);
            k22=Delta_r*Ymatrix*(y6(1:8,i,ilay)+B21*k21);
            %%%%%%%%% Step 3
            rK=r(1,i,ilay)+A3*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_ocean(l,rK,rho(ilay),gK);
            k13=Delta_r*Ymatrix*(y5(1:8,i,ilay)+B31*k11+B32*k12);
            k23=Delta_r*Ymatrix*(y6(1:8,i,ilay)+B31*k21+B32*k22);
            %%%%%%%%% Step 4
            rK=r(1,i,ilay)+A4*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_ocean(l,rK,rho(ilay),gK);
            k14=Delta_r*Ymatrix*(y5(1:8,i,ilay)+B41*k11+B42*k12+B43*k13);
            k24=Delta_r*Ymatrix*(y6(1:8,i,ilay)+B41*k21+B42*k22+B43*k23);
            %%%%%%%%% Step 5
            rK=r(1,i,ilay)+A5*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_ocean(l,rK,rho(ilay),gK);
            k15=Delta_r*Ymatrix*(y5(1:8,i,ilay)+B51*k11+B52*k12+B53*k13+B54*k14);
            k25=Delta_r*Ymatrix*(y6(1:8,i,ilay)+B51*k21+B52*k22+B53*k23+B54*k24);
            %%%%%%%%% Step 6
            rK=r(1,i,ilay)+A6*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_ocean(l,rK,rho(ilay),gK);
            k16=Delta_r*Ymatrix*(y5(1:8,i,ilay)+B61*k11+B62*k12+B63*k13+B64*k14+B65*k15);
            k26=Delta_r*Ymatrix*(y6(1:8,i,ilay)+B61*k21+B62*k22+B63*k23+B64*k24+B65*k25);
            %%%%%%%%%% Get next step
            y5(1:8,i+1,ilay)=y5(1:8,i,ilay)+AC1*k11+AC2*k12+AC3*k13+AC4*k14+AC5*k15+AC6*k16;
            y6(1:8,i+1,ilay)=y6(1:8,i,ilay)+AC1*k21+AC2*k22+AC3*k23+AC4*k24+AC5*k25+AC6*k26;   
            i=i+1;    
        end
        ifin(ilay)=i;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LAYERS ABOVE THE OCEAN 
    if ilay>layer_ocean
        if ilay>layer_ocean+1
        y7(1:6,1,ilay)=y7(1:6,ifin(ilay-1),ilay-1);
        y8(1:6,1,ilay)=y8(1:6,ifin(ilay-1),ilay-1);
        y9(1:6,1,ilay)=y9(1:6,ifin(ilay-1),ilay-1);
        y10(1:6,1,ilay)=y10(1:6,ifin(ilay-1),ilay-1);
        end
        %INTEGRATE
        Delta_r=(R(ilay+1)-R(ilay))/nrr;
        i=1;
        %INTEGRATE THE SOLUTION UPWARDS 
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
            k11=Delta_r*(Ymatrix*y7(1:8,i,ilay));
            k21=Delta_r*(Ymatrix*y8(1:8,i,ilay));
            k31=Delta_r*(Ymatrix*y9(1:8,i,ilay));
            k41=Delta_r*(Ymatrix*y10(1:8,i,ilay));
            %%%%%%%%% Step 2
            rK=r(1,i,ilay)+A2*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
            k12=Delta_r*Ymatrix*(y7(1:8,i,ilay)+B21*k11);
            k22=Delta_r*Ymatrix*(y8(1:8,i,ilay)+B21*k21);
            k32=Delta_r*Ymatrix*(y9(1:8,i,ilay)+B21*k31);
            k42=Delta_r*Ymatrix*(y10(1:8,i,ilay)+B21*k41);
            %%%%%%%%% Step 3
            rK=r(1,i,ilay)+A3*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
            k13=Delta_r*Ymatrix*(y7(1:8,i,ilay)+B31*k11+B32*k12);
            k23=Delta_r*Ymatrix*(y8(1:8,i,ilay)+B31*k21+B32*k22);
            k33=Delta_r*Ymatrix*(y9(1:8,i,ilay)+B31*k31+B32*k32);
            k43=Delta_r*Ymatrix*(y10(1:8,i,ilay)+B31*k41+B32*k42);
            %%%%%%%%% Step 4
            rK=r(1,i,ilay)+A4*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
            k14=Delta_r*Ymatrix*(y7(1:8,i,ilay)+B41*k11+B42*k12+B43*k13);
            k24=Delta_r*Ymatrix*(y8(1:8,i,ilay)+B41*k21+B42*k22+B43*k23);
            k34=Delta_r*Ymatrix*(y9(1:8,i,ilay)+B41*k31+B42*k32+B43*k33);
            k44=Delta_r*Ymatrix*(y10(1:8,i,ilay)+B41*k41+B42*k42+B43*k43);
            %%%%%%%%% Step 5
            rK=r(1,i,ilay)+A5*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
            k15=Delta_r*Ymatrix*(y7(1:8,i,ilay)+B51*k11+B52*k12+B53*k13+B54*k14);
            k25=Delta_r*Ymatrix*(y8(1:8,i,ilay)+B51*k21+B52*k22+B53*k23+B54*k24);
            k35=Delta_r*Ymatrix*(y9(1:8,i,ilay)+B51*k31+B52*k32+B53*k33+B54*k34);
            k45=Delta_r*Ymatrix*(y10(1:8,i,ilay)+B51*k41+B52*k42+B53*k43+B54*k44);
            %%%%%%%%% Step 6
            rK=r(1,i,ilay)+A6*Delta_r;
            gK=g(1,1,ilay)*(R(ilay)/rK)^2+4/3*pi*Gg*rho(ilay)*(rK^3-R(ilay)^3)/rK^2;
            Ymatrix=propagation_matrix_porosity(l,rK,rho(ilay),rhof(ilay),muC(ilay),lambda(ilay),Ks(ilay),Ku(ilay),Kd(ilay),Kf(ilay),alpha(ilay),poro(ilay),k_perm(ilay)/etaf(ilay),gK,omega,self_gravity,tidal_fluid);
            k16=Delta_r*Ymatrix*(y7(1:8,i,ilay)+B61*k11+B62*k12+B63*k13+B64*k14+B65*k15);
            k26=Delta_r*Ymatrix*(y8(1:8,i,ilay)+B61*k21+B62*k22+B63*k23+B64*k24+B65*k25);
            k36=Delta_r*Ymatrix*(y9(1:8,i,ilay)+B61*k31+B62*k32+B63*k33+B64*k34+B65*k35);
            k46=Delta_r*Ymatrix*(y10(1:8,i,ilay)+B61*k41+B62*k42+B63*k43+B64*k44+B65*k45);
            %%%%%%%%%% Get next step
            y7(1:8,i+1,ilay)=y7(1:8,i,ilay)+AC1*k11+AC2*k12+AC3*k13+AC4*k14+AC5*k15+AC6*k16;
            y8(1:8,i+1,ilay)=y8(1:8,i,ilay)+AC1*k21+AC2*k22+AC3*k23+AC4*k24+AC5*k25+AC6*k26;
            y9(1:8,i+1,ilay)=y9(1:8,i,ilay)+AC1*k31+AC2*k32+AC3*k33+AC4*k34+AC5*k35+AC6*k36;
            y10(1:8,i+1,ilay)=y10(1:8,i,ilay)+AC1*k41+AC2*k42+AC3*k43+AC4*k44+AC5*k45+AC6*k46;   
            i=i+1;    
        end
        ifin(ilay)=i;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF INTEGRATION        
end

%% BOUNDARY CONDITIONS
% Matrix with boundary conditions is assembled here. We need 10 BC to get
% C1,....C10.
    %x3 at the surface:                 ROWS 1,2,3
    %x2 at the ice-ocean interfacce     ROWS 4,5
    %x4 at the mantle-ocean interface   ROWS 6,7,8,9
    %x1 at the porous boundary          ROW 10
M=zeros(10,10);
B=zeros(10,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOUNDARY CONDITIONS AT THE SURFACE    
% ROW 2: tangent stress at the surface is 0
M(2,7)=y7(4,ifin(end),Nlayers);
M(2,8)=y8(4,ifin(end),Nlayers);
M(2,9)=y9(4,ifin(end),Nlayers);
M(2,10)=y10(4,ifin(end),Nlayers);

B=zeros(10,1);
if strain_BC==0 
    % ROW 1: Normal stress at the surface is 0
    B(1)=0;
    M(1,7)=y7(3,ifin(end),Nlayers);
    M(1,8)=y8(3,ifin(end),Nlayers);
    M(1,9)=y9(3,ifin(end),Nlayers);
    M(1,10)=y10(3,ifin(end),Nlayers);
    %ROW 3: potential stress at the surface  
    B(3)=(2*l+1)/Rs;
    M(3,7)=y7(6,ifin(end),Nlayers);
    M(3,8)=y8(6,ifin(end),Nlayers);
    M(3,9)=y9(6,ifin(end),Nlayers);
    M(3,10)=y10(6,ifin(end),Nlayers);
else %strain boundary condition 
    %ROW 1: STRAIN AT THE SURFACE  
    B(1)=strain_BC;
    M(1,7)=y7(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y7(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y7(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y7(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    M(1,8)=y8(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y8(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y8(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y8(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    M(1,9)=y9(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y9(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y9(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y9(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    M(1,10)=y10(3,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers))+lambda(Nlayers)*l*(l+1)*y10(2,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))-2*lambda(Nlayers)*y10(1,ifin(end),Nlayers)/((2*muC(Nlayers)+lambda(Nlayers)))+alpha(Nlayers)*y10(7,ifin(end),Nlayers)/(lambda(Nlayers)+2*muC(Nlayers));
    % ROW 2: NON-ROTATIONAL 
    M(2,7)=2*(y7(1,ifin(end),Nlayers)-y7(2,ifin(end),Nlayers))-y7(4,ifin(end),Nlayers)/(muC(Nlayers));
    M(2,8)=2*(y8(1,ifin(end),Nlayers)-y8(2,ifin(end),Nlayers))-y8(4,ifin(end),Nlayers)/(muC(Nlayers));
    M(2,9)=2*(y9(1,ifin(end),Nlayers)-y9(2,ifin(end),Nlayers))-y9(4,ifin(end),Nlayers)/(muC(Nlayers));
    M(2,10)=2*(y10(1,ifin(end),Nlayers)-y10(2,ifin(end),Nlayers))-y10(4,ifin(end),Nlayers)/(muC(Nlayers));     
    % ROW 3: Potential stress at the surface 
    B(3)=0;
    M(3,7)=y7(6,ifin(end),Nlayers);
    M(3,8)=y8(6,ifin(end),Nlayers);
    M(3,9)=y9(6,ifin(end),Nlayers);
    M(3,10)=y10(6,ifin(end),Nlayers);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOUNDARY CONDITIONS AT THE POROUS BOUNDARY
if (isnan(layporo)==1) %No porus layers
    M(10,1)=0;
    M(10,2)=0;
    M(10,3)=0;
    M(10,4)=1; 
else
    if isnan(pressure)==1  %no flux boundary condition
    %ROW 10: No flux boundary conditions 
    B(10)=0;
    M(10,1)=y1(8,ifin(layporo),layporo);
    M(10,2)=y2(8,ifin(layporo),layporo);
    M(10,3)=y3(8,ifin(layporo),layporo);
    M(10,4)=y4(8,ifin(layporo),layporo);
    else  %\sigma_rr+p=0
    B(10)=pressure;
    M(10,1)=y1(7,ifin(layporo),layporo)+y1(3,ifin(layporo),layporo);
    M(10,2)=y2(7,ifin(layporo),layporo)+y2(3,ifin(layporo),layporo);
    M(10,3)=y3(7,ifin(layporo),layporo)+y3(3,ifin(layporo),layporo);
    M(10,4)=y4(7,ifin(layporo),layporo)+y4(3,ifin(layporo),layporo);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%% BOUNDARY CONDITIONS AT THE ICE-OCEAN
% ROW 4: CONTINUITY OF GRAVITATIONAL POTENTIAL 
M(4,5)=y5(5,ifin(layer_ocean),layer_ocean);
M(4,6)=y6(5,ifin(layer_ocean),layer_ocean);
M(4,7)=-y7(5,1,layer_ocean+1);
M(4,8)=-y8(5,1,layer_ocean+1);
M(4,9)=-y9(5,1,layer_ocean+1);
M(4,10)=-y10(5,1,layer_ocean+1);
% ROW 5: CONTINUITY OF THE POTENTIAL STRESS (see equations 3.37-3.42
% Jara-Orue)
M(5,5)=y5(6,ifin(layer_ocean),layer_ocean)+4*pi*Gg*rho(layer_ocean)/glayer(layer_ocean+1)*y5(5,ifin(layer_ocean),layer_ocean);
M(5,6)=y6(6,ifin(layer_ocean),layer_ocean)+4*pi*Gg*rho(layer_ocean)/glayer(layer_ocean+1)*y6(5,ifin(layer_ocean),layer_ocean);
M(5,7)=-y7(6,1,layer_ocean+1)+4*pi*Gg*rho(layer_ocean)*y7(1,1,layer_ocean+1);
M(5,8)=-y8(6,1,layer_ocean+1)+4*pi*Gg*rho(layer_ocean)*y8(1,1,layer_ocean+1);
M(5,9)=-y9(6,1,layer_ocean+1)+4*pi*Gg*rho(layer_ocean)*y9(1,1,layer_ocean+1);
M(5,10)=-y10(6,1,layer_ocean+1)+4*pi*Gg*rho(layer_ocean)*y10(1,1,layer_ocean+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOUNDARY CONDITIONS AT THE OCEAN-MANTLE INTERFACE
% ROW 6: CONTINUITY OF THE GRAVITATIONAL POTENTIAL
M(6,1)=-y1(5,ifin(layer_ocean-1),layer_ocean-1);
M(6,2)=-y2(5,ifin(layer_ocean-1),layer_ocean-1);
M(6,3)=-y3(5,ifin(layer_ocean-1),layer_ocean-1);
M(6,4)=-y4(5,ifin(layer_ocean-1),layer_ocean-1);
M(6,5)=y5(5,1,layer_ocean);
M(6,6)=y6(5,1,layer_ocean);
% ROW 7: CONTINUITY OF THE POTENTIAL STRESS
M(7,1)=-y1(6,ifin(layer_ocean-1),layer_ocean-1)+4*pi*Gg*rho(layer_ocean)*y1(1,ifin(layer_ocean-1),layer_ocean-1);
M(7,2)=-y2(6,ifin(layer_ocean-1),layer_ocean-1)+4*pi*Gg*rho(layer_ocean)*y2(1,ifin(layer_ocean-1),layer_ocean-1);
M(7,3)=-y3(6,ifin(layer_ocean-1),layer_ocean-1)+4*pi*Gg*rho(layer_ocean)*y3(1,ifin(layer_ocean-1),layer_ocean-1);
M(7,4)=-y4(6,ifin(layer_ocean-1),layer_ocean-1)+4*pi*Gg*rho(layer_ocean)*y4(1,ifin(layer_ocean-1),layer_ocean-1);
M(7,5)=y5(6,1,layer_ocean)+4*pi*Gg*rho(layer_ocean)/glayer(layer_ocean)*y5(5,1,layer_ocean);
M(7,6)=y6(6,1,layer_ocean)+4*pi*Gg*rho(layer_ocean)/glayer(layer_ocean)*y6(5,1,layer_ocean);
% ROW 8: TANGENTIAL STRESS
M(8,1)=y1(4,ifin(layer_ocean-1),layer_ocean-1);
M(8,2)=y2(4,ifin(layer_ocean-1),layer_ocean-1);
M(8,3)=y3(4,ifin(layer_ocean-1),layer_ocean-1);
M(8,4)=y4(4,ifin(layer_ocean-1),layer_ocean-1);
% ROW 9: BOUYANCY CONDITION 
M(9,1)=y1(1,ifin(layer_ocean-1),layer_ocean-1)-y1(3,ifin(layer_ocean-1),layer_ocean-1)/(rho(layer_ocean)*glayer(layer_ocean))+y1(5,ifin(layer_ocean-1),layer_ocean-1)/(glayer(layer_ocean));
M(9,2)=y2(1,ifin(layer_ocean-1),layer_ocean-1)-y2(3,ifin(layer_ocean-1),layer_ocean-1)/(rho(layer_ocean)*glayer(layer_ocean))+y2(5,ifin(layer_ocean-1),layer_ocean-1)/(glayer(layer_ocean));
M(9,3)=y3(1,ifin(layer_ocean-1),layer_ocean-1)-y3(3,ifin(layer_ocean-1),layer_ocean-1)/(rho(layer_ocean)*glayer(layer_ocean))+y3(5,ifin(layer_ocean-1),layer_ocean-1)/(glayer(layer_ocean));
M(9,4)=y4(1,ifin(layer_ocean-1),layer_ocean-1)-y4(3,ifin(layer_ocean-1),layer_ocean-1)/(rho(layer_ocean)*glayer(layer_ocean))+y4(5,ifin(layer_ocean-1),layer_ocean-1)/(glayer(layer_ocean));
%% INVERT MATRIX 
sol=M\B;
%sol=pinv(M)*B;
%% COMPUTE THE SOLUTION
y=sol(1)*y1+sol(2)*y2+sol(3)*y3+sol(4)*y4+sol(5)*y5+sol(6)*y6+sol(7)*y7+sol(8)*y8+sol(9)*y9+sol(10)*y10;

% sol=M(1:3,1:3)\B(1:3);
% y=sol(1)*y1(1:8,1:ifin)+sol(2)*y2(1:8,1:ifin)+sol(3)*y3(1:8,1:ifin);
if resample>0
        yres=y(:,1:resample:end,:);
        y=yres;
        rres=r(:,1:resample:end,:);
        r=rres;
end

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
    %% CHECK BOUNDARY CONDITIONS ARE OK
    cprintf('*Blue' , 'Check that B.C are ok. All results should be 0\n')
    % SURFACE
    cprintf('*Blue' , 'Surface Boundary Conditions \n')
    disp(['Normal stress '   num2str(y(3,ifin(Nlayers),Nlayers), '%10.5e')])
    disp(['Tangential stress '   num2str(y(4,ifin(Nlayers),Nlayers), '%10.5e')])
    disp(['Potential stress '   num2str(y(6,ifin(Nlayers),Nlayers), '%10.5e')])
    % OCEAN-ICE
    cprintf('*Blue' , 'Ocean-ice\n')
    disp(['Tangential Stress '   num2str(y(4,1,layer_ocean+1), '%10.5e')])
    disp(['Buyoancy  '   num2str(y(1,1,layer_ocean+1)-y(3,1,layer_ocean+1)/(rho(layer_ocean)*glayer(layer_ocean+1))+y(5,1,layer_ocean+1)/glayer(layer_ocean+1), '%10.5e')])
    disp(['Continuity of grav potential '   num2str(y(5,1,layer_ocean+1)-y(5,ifin(layer_ocean),layer_ocean), '%10.5e')])
    disp(['Potential Stress '  num2str(y(6,ifin(layer_ocean),layer_ocean)-y(6,1,layer_ocean+1)+4*pi*Gg*rho(layer_ocean)*(y(1,1,layer_ocean+1)+y(5,ifin(layer_ocean),layer_ocean)/glayer(layer_ocean+1)), '%10.5e')])
    % OCEAN-MANTLE
    cprintf('*Blue' , 'Ocean-mantle\n')
    disp(['Tangential Stress '   num2str(y(4,ifin(layer_ocean-1),layer_ocean-1), '%10.5e')])
    disp(['Buyoancy  '   num2str(y(1,ifin(layer_ocean-1),layer_ocean-1)-y(3,ifin(layer_ocean-1),layer_ocean-1)/(rho(layer_ocean)*glayer(layer_ocean))+y(5,ifin(layer_ocean-1),layer_ocean-1)/glayer(layer_ocean), '%10.5e')])
    disp(['Continuity of grav potential '   num2str(y(5,1,layer_ocean)-y(5,ifin(layer_ocean-1),layer_ocean-1), '%10.5e')])
    disp(['Potential Stress '  num2str(y(6,1,layer_ocean)-y(6,ifin(layer_ocean-1),layer_ocean-1)+4*pi*Gg*rho(layer_ocean)*(y(1,ifin(layer_ocean-1),layer_ocean-1)+y(5,1,layer_ocean)/glayer(layer_ocean)), '%10.5e')])
    %% PLOTS FOR DEBUG
%     figure
%     title("y1")
%     for i=1:8
%         subplot(2,4,i)
%         for j=layer_start:Nlayers
%         plot(r(1,1:ifin(j),j),abs(y1(i,1:ifin(j),j)));
%         hold on
%         end
%         set(gca, 'YScale', 'log')
%     end
%     figure
%     title("y2")
%     for i=1:8
%     subplot(2,4,i)
%         for j=layer_start:Nlayers
%         plot(r(1,1:ifin(j),j),abs(y2(i,1:ifin(j),j)));
%         hold on
%         end
%     set(gca, 'YScale', 'log')
%     end
%     figure
%     for i=1:8
%     subplot(2,4,i)
%         for j=layer_start:Nlayers
%         plot(r(1,1:ifin(j),j),abs(y3(i,1:ifin(j),j)));
%         hold on 
%         end
%     set(gca, 'YScale', 'log')
%     end
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
    title("y4")
    for i=1:8
    subplot(2,4,i)
        for j=layer_start:Nlayers
        plot(r(1,1:ifin(j),j),abs(y3(i,1:ifin(j),j)),'color','r');
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
    
    figure
    for j=layer_start:Nlayers
        plot(r(1,1:ifin(j),j),real(y(7,1:ifin(j),j))+rhof(j)/rho(end)*real(y(5,1:ifin(j),j)),'LineWidth',2);
        hold on        
    end
    title("y7+rhoy5")
end


end
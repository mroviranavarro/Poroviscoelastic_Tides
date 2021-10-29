%% FUNCTION: propagation_matrix_porosity
% Description: Compute the propagation matrix for the poroviscoelastic normal mode problem  
% Author: M. Rovira-Navarro 
%% INPUT 
%lin: spherical harmonic degree
%rL: radial position
%rhoL: average density
%rhofL: fluid density
%muL: shear modulus
%lambdaL: Lame parameter
%KsL: solid bulk modulus
%KuL: undrained bulk modulus
%KdL: drained bulk modulus
%KfL: fluid bulk modulus
%alphaL: alpha
%poroL: porosity
%k_perm2L: permeability/viscosity
%gL: gravity 
%omega: forcing frequency
%self_gravity: (1) if self-gravity is used in the momentum equation, (0) if not
%tidal_fluid: (1) if tidal potential affects the fluid, (0) if not



%% OUTPUT
%Y(8,8): propagator matrix 
%% Start of the function
function Y=propagation_matrix_porosity(l,rL,rhoL,rhofL,muL,lambdaL,KsL,KuL,KdL,KfL,alphaL,poroL,k_perm2L,gL,omega,self_gravity,tidal_fluid)
global Gg
Y=zeros(8,8);
K=-4*pi*Gg*rhofL*k_perm2L/(omega*1i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(1,1)=-2*lambdaL/((2*muL+lambdaL)*rL);
Y(1,2)=l*(l+1)*lambdaL/((2*muL+lambdaL)*rL);
Y(1,3)=1/(2*muL+lambdaL);
Y(1,4)=0;
Y(1,5)=0;
Y(1,6)=0;
Y(1,7)=alphaL/(2*muL+lambdaL);
Y(1,8)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(2,1)=-1/rL;
Y(2,2)=1/rL;
Y(2,3)=0;
Y(2,4)=1/muL;
Y(2,5)=0;
Y(2,6)=0;
Y(2,7)=0;
Y(2,8)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(3,1)=4/rL*(muL/rL+2*muL*lambdaL/((2*muL+lambdaL)*rL)-rhoL*gL/2-lambdaL*gL*rhoL/(2*(2*muL+lambdaL))) ...
	+self_gravity*(rhofL*alphaL-rhoL)*2*gL*(-lambdaL/(2*muL+lambdaL)+1)/rL;
Y(3,2)=-l*(l+1)/rL*(-rhoL*gL*lambdaL/(2*muL+lambdaL)+2*muL/rL+4*muL*lambdaL/(rL*(2*muL+lambdaL))) ...
    +self_gravity*(rhofL*alphaL-rhoL)*gL*l*(l+1)*(lambdaL/(2*muL+lambdaL)-1)/rL;
Y(3,3)=1/(2*muL+lambdaL)*(-4*muL/rL+rhoL*gL) ...
    +self_gravity*(rhofL*alphaL-rhoL)*gL/(2*muL+lambdaL);
Y(3,4)=l*(l+1)/rL;
Y(3,5)=-rhoL*(l+1)/rL;
Y(3,6)=rhoL;
Y(3,7)=alphaL/(lambdaL+2*muL)*(-4*muL/rL+rhoL*gL) ...
    +self_gravity*(rhofL*alphaL-rhoL)*gL*alphaL/(2*muL+lambdaL)+self_gravity*gL*rhofL*(poroL/KfL+(alphaL-poroL)/(KsL));
Y(3,8)=-K*rhoL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(4,1)=-1/rL*((2+4*lambdaL/(2*muL+lambdaL))*muL/rL-rhoL*gL);
Y(4,2)=2*muL/rL^2*(l^2+l-1+lambdaL*l*(l+1)/(2*muL+lambdaL));
Y(4,3)=-1/rL*(1-2*muL/(2*muL+lambdaL));
Y(4,4)=-3/rL;
Y(4,5)=rhoL/rL;
Y(4,6)=0;
Y(4,7)=2*alphaL*muL/(rL*(lambdaL+2*muL));
Y(4,8)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(5,1)=-4*pi*Gg*rhoL;
Y(5,2)=0;
Y(5,3)=0;
Y(5,4)=0;
Y(5,5)=-(l+1)/rL;
Y(5,6)=1;
Y(5,7)=0;
Y(5,8)=-K;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(6,1)=-4*pi*Gg*rhoL/rL*(l-1+2*lambdaL/(2*muL+lambdaL)) ...
        +self_gravity*8*pi*Gg*(rhofL*alphaL-rhoL)*(-lambdaL/(2*muL+lambdaL)+1)/rL ...
        +K*1i*2*alphaL*omega/(k_perm2L*rL)*(1-lambdaL/(2*muL+lambdaL));
Y(6,2)=4*pi*Gg*rhoL*l*(l+1)*lambdaL/((2*muL+lambdaL)*rL) ...
        +self_gravity*4*pi*Gg*(rhofL*alphaL-rhoL)*l*(l+1)*(lambdaL/(2*muL+lambdaL)-1)/rL ...
        +K*1i*omega*alphaL*l*(l+1)/(k_perm2L*rL)*(-1+lambdaL/(lambdaL+2*muL));
Y(6,3)=4*pi*Gg*rhoL/(2*muL+lambdaL) ...
        +self_gravity*4*pi*Gg*(rhofL*alphaL-rhoL)/(2*muL+lambdaL) ...
        +K*1i*omega*alphaL/(k_perm2L*(2*muL+lambdaL));
Y(6,4)=0;
Y(6,5)=K*tidal_fluid*l*(l+1)*rhofL/rL^2;
Y(6,6)=(l-1)/rL;
Y(6,7)=4*pi*Gg*rhoL*alphaL/(lambdaL+2*muL) ...
        +self_gravity*4*pi*Gg*((rhofL*alphaL-rhoL)*alphaL/(2*muL+lambdaL)+rhofL*(poroL/KfL+(alphaL-poroL)/(KsL)))...
        +K*(l*(l+1)/rL^2+1i*omega/k_perm2L*(alphaL^2/(2*muL+lambdaL)+poroL/KfL+(alphaL-poroL)/(KsL)));
Y(6,8)=-K*(l-1)/rL-K*2/rL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(7,1)=tidal_fluid*4*pi*Gg*rhofL*rhoL;
Y(7,2)=0;
Y(7,3)=0;
Y(7,4)=0;
Y(7,5)=tidal_fluid*rhofL*(l+1)/rL; %old tidal_fluid*rhofL*l*(l+1)/rL
Y(7,6)=-tidal_fluid*rhofL;
Y(7,7)=-rhofL*gL/KfL;
Y(7,8)=1+K*tidal_fluid*rhofL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(8,1)=1i*2*alphaL*omega/(k_perm2L*rL)*(1-lambdaL/(2*muL+lambdaL));
Y(8,2)=1i*omega*alphaL*l*(l+1)/(k_perm2L*rL)*(-1+lambdaL/(lambdaL+2*muL));
Y(8,3)=1i*omega*alphaL/(k_perm2L*(2*muL+lambdaL));
Y(8,4)=0;
Y(8,5)=tidal_fluid*l*(l+1)*rhofL/rL^2;
Y(8,6)=0;
Y(8,7)=l*(l+1)/rL^2+1i*omega/k_perm2L*(alphaL^2/(2*muL+lambdaL)+poroL/KfL+(alphaL-poroL)/(KsL));
Y(8,8)=-2/rL;

end 
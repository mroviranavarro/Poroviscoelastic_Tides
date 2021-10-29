%% FUNCTION: propagation_matrix_ocean
% Description: Compute the propagation matrix for the ocean layer 
% Assumption: Incompressible ocean (Eqs. 3.17-3.18 Jara Orue thesis)
% Author: M. Rovira-Navarro 
%% INPUT 
%lin: spherical harmonic degree
%rL: radial position
%rhoL: density of the mantle 
%gL: gravity 


%% OUTPUT
%Y(8,8): propagator matrix 
%% Start of the function
function Y=propagation_matrix_ocean(l,rL,rhoL,gL)
global Gg
Y=zeros(8,8);
% Only z5 and z6 are propagated in the ocean (incompressible assumed) 
%%%%%%%%%%%%%%%%%%%%%%%
Y(5,5)=4*pi*Gg*rhoL/gL-(l+1)/rL;
Y(5,6)=1;
%%%%%%%%%%%%%%%%%%%%%%%
Y(6,5)=2*(l-1)*4*pi*Gg*rhoL/(rL*gL);
Y(6,6)=-(4*pi*Gg*rhoL/gL-(l-1)/rL);


end 
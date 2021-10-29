%% FUNCTION: compute_energy
% Description: Given the strain, stress, variation of fluid content and pore pressure compute energy dissipated  
% Author: M. Rovira-Navarro 
%% INPUT 
%strain: strain tensor 
                %strain(icolat,ilon,ir,1)=\epsilon_r_r;
                %strain(icolat,ilon,ir,2)=\epsilon_theta_theta;
                %strain(icolat,ilon,ir,3)=\epsilon_phi_phi;
                %strain(icolat,ilon,ir,4)=\epsilon_r_theta;
                %strain(icolat,ilon,ir,5)=\epsilon_r_phi;
                %strain(icolat,ilon,ir,6)=\epsilon_theta_phi;  
%stress: stress(icolat,ilon,ir,1)=\sigma_r_r;
                %stress(icolat,ilon,ir,2)=\sigma_theta_theta;
                %stress(icolat,ilon,ir,3)=\sigma_phi_phi;
                %stress(icolat,ilon,ir,4)=\sigma_r_theta;
                %stress(icolat,ilon,ir,5)=\sigma_r_phi;
                %stress(icolat,ilon,ir,6)=\sigma_theta_phi;
% flux: flux vector
                % flux(icolat,ilon,ir,1): radial component of the flux
                % flux(icolat,ilon,ir,2): latitudinal component of the flux
                % flux(icolat,ilon,ir,3): longitudinal component of the flux
%p_fluid: fluid pore pressure
                %p_fluid(icolat,ilon,ir,1)
%C_fluid: variation of fluid content: 
                %C_fluid(icolat,ilon,ir,1)
%omega: forcing frequency 
%etaf:  viscosity of the fluid
%k_perm: permeability 

%% OUTPUT 
%energy_solid: volumetric energy dissipated in the solid matrix 
            %energy_solid(icolat,ilon,ir)--- Eq. (19a)
%energy_solid_pore: energy dissipated in the solid due to pore pressure 
            %energy_solid_pore(icolat,ilon,ir)--- Eq. (19a) term with pressure
%energy_fluid: energy dissipated in due to Darcy's flow 
            %energy_fluid(icolat,ilon,ir,1)--- Eq. (19b)
%energy_solid_surface: energy dissipated in the solid matrix integrated radially 
            %energy_solid_surface(icolat,ilon)
%energy_solid_pore_surface: energy dissipated in the solid due to pore pressure integrated radially  
            %energy_solid_pore_surface(icolat,ilon)
%energy_fluid_surface: energy dissipated in due to Darcy's flow integrated radially  
            %energy_fluid_surface(icolat,ilon,1)            
%energy_solid_total: Total energy dissipated in the solid integrated
%energy_fluid_total: Total energy dissipated in the fluid

%% START OF THE FUNCTION 
function [energy_solid, energy_solid_pore,energy_fluid,energy_solid_surface, energy_solid_pore_surface,energy_fluid_surface,energy_solid_total, energy_fluid_total]=compute_energy(strain,stress,flux,p_fluid,C_fluid,omega,etaf,k_perm,r,colat,lon,poro)
%% Volumetric (W/m^3)
%dissipation in the solid due to strain*stress
%rr (x1)
en_rr=real(stress(:,:,:,1)).*imag(strain(:,:,:,1))-imag(stress(:,:,:,1)).*real(strain(:,:,:,1));
%thetatheta (x1)
en_thetatheta=real(stress(:,:,:,2)).*imag(strain(:,:,:,2))-imag(stress(:,:,:,2)).*real(strain(:,:,:,2));
%phiphi (x1)
en_phiphi=real(stress(:,:,:,3)).*imag(strain(:,:,:,3))-imag(stress(:,:,:,3)).*real(strain(:,:,:,3));
%rtheta (x2)
en_rtheta=real(stress(:,:,:,4)).*imag(strain(:,:,:,4))-imag(stress(:,:,:,4)).*real(strain(:,:,:,4));
%rphi (x2)
en_rphi=real(stress(:,:,:,5)).*imag(strain(:,:,:,5))-imag(stress(:,:,:,5)).*real(strain(:,:,:,5));
%thetaphi (x2)
en_thetaphi=real(stress(:,:,:,6)).*imag(strain(:,:,:,6))-imag(stress(:,:,:,6)).*real(strain(:,:,:,6));
% add components together
energy_solid_strain_stress=-omega/2*(en_rr+en_thetatheta+en_phiphi+2*en_rtheta+2*en_rphi+2*en_thetaphi);
%dissipation in the solid due to pore pressure
energy_solid_pore=-omega/2*(real(p_fluid).*imag(C_fluid)-imag(p_fluid).*real(C_fluid));
%dissipation in the liquid
lay_poro=find(poro>0);
if isempty(lay_poro)==0
energy_fluid=1/2*etaf(lay_poro)/k_perm(lay_poro)*(abs(flux(:,:,:,1)).^2+abs(flux(:,:,:,2)).^2+abs(flux(:,:,:,3)).^2);
else
energy_fluid=1/2*etaf(1)/k_perm(1)*(abs(flux(:,:,:,1)).^2+abs(flux(:,:,:,2)).^2+abs(flux(:,:,:,3)).^2);
end
% energy in the solid
energy_solid=energy_solid_strain_stress+energy_solid_pore;
%% Integrated to the surface (W/m^2)
% integrate upwards
e_aux=(energy_solid(:,:,1:end-1)+energy_solid(:,:,2:end))/2;
r_aux=((r(1:end-1)+r(2:end))/2)';
Delta_r=diff(r)';
r_aux2=r_aux.^2.*Delta_r;
r_aux3(1,1,:)=r_aux2;
e_aux2=e_aux.*r_aux3;
energy_solid_surface=sum(e_aux2,3);
%pore
e_aux=(energy_solid_pore(:,:,1:end-1)+energy_solid_pore(:,:,2:end))/2;
e_aux2=e_aux.*r_aux3;
energy_solid_pore_surface=sum(e_aux2,3);
%fluid
e_aux=(energy_fluid(:,:,1:end-1)+energy_fluid(:,:,2:end))/2;
e_aux2=e_aux.*r_aux3;
energy_fluid_surface=sum(e_aux2,3);
%% Integrated total (W)
%compute average energy
% solid
e_aux=(energy_solid_surface(:,1:end-1)+energy_solid_surface(:,2:end))/2;
e_aux=(e_aux(1:end-1,:)+e_aux(2:end,:))/2;
colat_aux=(colat(1:end-1)+colat(2:end))/2;
colat_aux=repmat(colat_aux',1,length(e_aux(1,:)));
colat_dif=repmat(diff(colat)',1,length(e_aux(1,:)));
lon_dif=repmat(diff(lon),length(e_aux(:,1)),1);
e_aux=e_aux.*lon_dif.*colat_dif.*sin(colat_aux);
energy_solid_total=sum(e_aux(:));
%pore
e_aux=(energy_solid_pore_surface(:,1:end-1)+energy_solid_pore_surface(:,2:end))/2;
e_aux=(e_aux(1:end-1,:)+e_aux(2:end,:))/2;
e_aux=e_aux.*lon_dif.*colat_dif.*sin(colat_aux);
energy_solid_pore_total=sum(e_aux(:));
%fluid
e_aux=(energy_fluid_surface(:,1:end-1)+energy_fluid_surface(:,2:end))/2;
e_aux=(e_aux(1:end-1,:)+e_aux(2:end,:))/2;
e_aux=e_aux.*lon_dif.*colat_dif.*sin(colat_aux);
energy_fluid_total=sum(e_aux(:));
end
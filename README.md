# Poroviscoelastic_Tides

**AUTHOR**: M. Rovira-Navarro.  
**DATE**: October 2021.   
**CODE**: Code accompanying the paper: "The Tides of Enceladus' Core". It can be used to btain the viscoelastic or poroviscoelastic tidal response of a moon due to eccentricity tides.   
The following information is given here:     
1. INCLUDED FUNCTIONS: Functions included, their use, input and outputs   
2. HOW TO USE? Description of how the code is normally used.  
3.  EXAMPLES: Description of the examples accompanying the code.    


## INCLUDED FUNCTIONS

#### tidal.m  
**Description**: Propagate the solution from the center to the surface    
**Inputs**:   
Required  
- l: spherical harmonic degree.
- R: vector containing the upper boundary of each layer. 
- rho: vector containing the average density of each layer. 
- rhof: vector containing the fluid average density of each layer. 
- mu: vector containing the shear modulus of each layer. 
- Ks: vector containing the solid bulk modulus of each layer. 
- etas: vector containing the viscosity of the solid.  
- alpha: vector conatining alpha. 
- poro: vector containing porosity for each layer.  
- k_perm: vector containing permeability for each layer. 
- etaf: vector containing fluid viscosity for each layer. 
- Kf: vector containing the fluid bulk modulus of each layer. 
- liquid: vector stating if layer is fluid (1) or not (0). Here used to indicate the model has a fluid core. 
- omega: forcing frequency. 

Optional
- radial_points: number of radial points
- resample: To solve the PDE it is good to have  a lot of points in the  radial direction, but then to build the solution u need to deal with big matrices that slows things down. In such case the resample option is useful. radial_points X lat_points X lon_points
- pressure_BC: use pressure boundary condition instead of constant flux
- strain_BC: use prescribed strain as in Liao et al. 2020
- print_results: print some results on the screen solution at the surface, Love numbers and y plots
- self_gravity: (1) if self-gravity is used in the momentum equation, (0) if not
- tidal_fluid: (1) if tidal potential affects the fluid, (0) if not
- gravity_on: turn gravity on (1) and off (0). Default 1


**Outputs**
- y(1:8,1:nrr,1:Nlayers): solution vector for each layer, y functions are given in Appendix A
	- y1: normal displacement 
	- y2: tangential displacements
	- y3: normal stress
	- y4: tangential stress
	- y5: perturbing potential 
	- y6: potential stress
	- y7: pore pressure
	- y8: radial flux
- r: radial points where solution is obtained 


#### tidal_ocean.m 
Description: Propagate the solution from the core to the surface but with the possibility of a subsurface ocean 

**Inputs: **
Required  
- l: spherical harmonic degree
- R: vector containing the upper boundary of each layer
- rho: vector containing the average density of each layer
- rhof: vector containing the fluid average density of each layer
- mu: vector containing the shear modulus of each layer
- Ks: vector containing the solid bulk modulus of each layer
- etas: vector containing the viscosity of the solid 
- alpha: vector conatining alpha
- poro: vector containing porosity for each layer 
- k_perm: vector containing permeability for each layer
- etaf: vector containing fluid viscosity for each layer
- Kf: vector containing the fluid bulk modulus of each layer
- liquid: vector stating if layer is fluid (1) or not (0). Here used to indicate the model has a fluid core, or a subsurface oceam
- omega: forcing frequency

Optional 
- radial_points: number of radial points
- resample: To solve the PDE it is good to have  a lot of points in the  radial direction, but then to build the solution u need to deal with big matrices that slows things down. In such case the resample option is useful.
- pressure_BC: use pressure boundary condition instead of constant flux
- strain_BC: use strain as boundary condition as in Liao et al. 2020
- print_results: print some results on the screen solution at the surface, love numbers and y plots
- self_gravity: (1) if self-gravity is used in the momentum equation, (0) if not
- tidal_fluid: (1) if tidal potential affects the fluid, (0) if not
- gravity_on: turn gravity on (1) and off (0). Default 1

**Outputs**
- y(1:8,1:nrr,1:Nlayers): solution vector for each layer,  y functions are given in Appendix A
	- y1: normal displacement 
	- y2: tangential displacements
	- y3: normal stress
	- y4: tangential stress
	- y5: perturbing potential 
	- y6: potential stress
	- y7: pore pressure
	- y8: radial flux
- r: radial points where solution is obtained 

#### propagation_matrix_porosity.m 
Description: Compute the propagation matrix for the poroviscoelastic normal mode problem  (Eq. (B1))

**Inputs: **
- lin: spherical harmonic degree
- rL: radial position
- rhoL: average density
- rhofL: fluid density
- muL: shear modulus
- lambdaL: Lame parameter
- KsL: solid bulk modulus
- KuL: undrained bulk modulus
- KdL: drained bulk modulus
- KfL: fluid bulk modulus
- alphaL: alpha
- poroL: porosity
- k_perm2L: permeability/viscosity
- gL: gravity 
- omega: forcing frequency
- self_gravity: (1) if self-gravity is used in the momentum equation, (0) if not
- tidal_fluid: (1) if tidal potential affects the fluid, (0) if not

**Outputs:**
- Y: Propagation matrix 

#### propagation_matrix_ocean.m 
**Description:** Compute the propagation matrix for the ocean layer (Eq. (B2))  
**Inputs:** 
- lin: spherical harmonic degree
- rL: radial position
- rhoL: density of the mantle 
- gL: gravity 
**Outputs:**
- Y: Propagation matrix 
	
#### build_solution.m  
**Description:** Given the solution vector y(r), compute the solution in a colat-lon-r grid for a given order and degree (see Appendix C) 

**Inputs:** 
Required
- y: solution vector y(1:8,1:nrr,1:Nlayers)
- r: radial points corresponding to y(1:8,1:nrr,1:Nlayers)
- rhof: fluid density
- rhos: solid density
- Ks: bulk modulus of the solid 
- Kf: bulk modulus of the fluid
- mu: shear modulus of the solid
- etas: viscosity of the solid	
- etaf: viscosity of the fluid
- liquid: 1 if layer is liquid
- k_perm: permeability 
- w_moon: forcing frequency of the moon
- alpha: Biot constant 
- poro: porosity
- l: degree of the tidal forcing
- m: order of the tidal forcing: %0,-2 or 2.

Optional 
- tidal_fluid: (1) if tides affect the fluid, (0) if not. 
- lat_points: number of points used in the latitide,longitude grid. Default 70

**Outputs:**
- colat: colatitude where solution is given 
- lon: longitudes where solution is given 
- rr: radial points where solution is returned
- displacements: displacement vector  
	- displacements(icolat,ilon,ir,1): radial component of the displacement
	- displacements(icolat,ilon,ir,2): latitudinal component of the
	- displacement
	- displacements(icolat,ilon,ir,3): longitudinal component of the diplacement
- flux: flux vector
	- flux(icolat,ilon,ir,1): radial component of the flux
	- flux(icolat,ilon,ir,2): latitudinal component of the flux
	- flux(icolat,ilon,ir,3): longitudinal component of the flux
- stress:
	- stress(icolat,ilon,ir,1)=\sigma_r_r;
	- stress(icolat,ilon,ir,2)=\sigma_theta_theta;
	- stress(icolat,ilon,ir,3)=\sigma_phi_phi;
	- stress(icolat,ilon,ir,4)=\sigma_r_theta;
	- stress(icolat,ilon,ir,5)=\sigma_r_phi;
	- stress(icolat,ilon,ir,6)=\sigma_theta_phi; 
- strain: strain tensor 
	- strain(icolat,ilon,ir,1)=\epsilon_r_r;
	- strain(icolat,ilon,ir,2)=\epsilon_theta_theta;
	- strain(icolat,ilon,ir,3)=\epsilon_phi_phi;
	- strain(icolat,ilon,ir,4)=\epsilon_r_theta;
	- strain(icolat,ilon,ir,5)=\epsilon_r_phi;
	- strain(icolat,ilon,ir,6)=\epsilon_theta_phi;                 
- gravpot: perturbing gravitational potential, gravpot(icolat,ilon,ir,1)
- p_fluid: fluid pore pressure, p_fluid(icolat,ilon,ir,1)
- C_fluid: variation of fluid content, C_fluid(icolat,ilon,ir,1)
- varargout can be used to get some extra outputs
	- varargout{1}: porosity change 
	- varargout{2}: divergence of the displacements
	- varargout{3}: divergence of the flux
                

#### compute_energy.m 

**Description:** Given the strain, stress, variation of fluid content and pore pressure compute energy dissipated 
**Inputs:** 
- strain: strain tensor 
	- strain(icolat,ilon,ir,1)=\epsilon_r_r;
	- strain(icolat,ilon,ir,2)=\epsilon_theta_theta;
	- strain(icolat,ilon,ir,3)=\epsilon_phi_phi;
	- strain(icolat,ilon,ir,4)=\epsilon_r_theta;
	- strain(icolat,ilon,ir,5)=\epsilon_r_phi;
	- strain(icolat,ilon,ir,6)=\epsilon_theta_phi;  
- stress: stress(icolat,ilon,ir,1)=\sigma_r_r;
	- stress(icolat,ilon,ir,2)=\sigma_theta_theta;
	- stress(icolat,ilon,ir,3)=\sigma_phi_phi;
	- stress(icolat,ilon,ir,4)=\sigma_r_theta;
	- stress(icolat,ilon,ir,5)=\sigma_r_phi;
	- stress(icolat,ilon,ir,6)=\sigma_theta_phi;
- flux: flux vector
	- flux(icolat,ilon,ir,1): radial component of the flux
	- flux(icolat,ilon,ir,2): latitudinal component of the flux
	- flux(icolat,ilon,ir,3): longitudinal component of the flux
- p_fluid: fluid pore pressure, p_fluid(icolat,ilon,ir,1)
- C_fluid: variation of fluid content, C_fluid(icolat,ilon,ir,1)
- omega: forcing frequency 
- etaf:  viscosity of the fluid
- k_perm: permeability

**Outputs:**
- energy_solid(icolat,ilon,ir): volumetric energy dissipated in the solid matrix, computed using Eq. (19a)
- energy_solid_pore(icolat,ilon,ir): energy dissipated in the solid due to pore pressure, Eq. (19a) term with pressure
- energy_fluid(icolat,ilon,ir,1): energy dissipated in due to Darcy's flow, computed using  Eq. (19b)
- energy_solid_surface energy_solid_surface(icolat,ilon): energy dissipated in the solid matrix integrated radially 
-  energy_solid_pore_surface(icolat,ilon): energy dissipated in the solid due to pore pressure integrated radially
- energy_fluid_surface(icolat,ilon): energy dissipated in due to Darcy's flow integrated radially
- energy_solid_total: Total energy dissipated in the solid integrated
- energy_fluid_total: Total energy dissipated in the fluid


## HOW TO USE? 

1. Define interior model 
2. Obtain radial functions using tidal.m or tidal_ocean
3. Obtain the complex fields a(colat,lon,r)
4. Compute tidal dissipation

## EXAMPLES 

The following example scripts are included

#### Europa.m:
Example to compute the tidal response of a viscoelastic Europa. The interior model of Europa is given by that of Beuthe, 2013 "Spatial patterns of tidal dissipation" https://doi.org/10.1016/j.icarus.2012.11.020
% We use the interior structure given in Table 5

#### Io.m
Example to compute the tidal response of a viscoelastic Io. 
The interior model of Io corresponds to models A and B of Steinke et al. 2020. https://doi.org/10.1016/j.icarus.2019.05.001.

#### Enceladus_Only_Core.m
Example to compute the proviscoelastic response of Enceladus. 
We consider the model of Liao et al. 2020: https://doi.org/10.1029/2019JE006209. In such case the boundary conditions are given by Eq. (A16)

#### Enceladus.m
Example to compute the proviscoelastic response of Enceladus. Interior model parameters given in Table 1 of the paper

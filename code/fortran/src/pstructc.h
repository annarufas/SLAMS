! Attributes of a particle, or the Lagrangian state variables. 
! NOTE: a cluster is made of identical (clone) particles.

integer(C_INT) :: id,          & ! identification number
				  initType,    & ! initial type: 1=faecal aggregate, 2=non-faecal aggregate, 3=living phyto, 4=dead phyto, 5=TEP, 6=clay, 7=dead zoo, 8=inorganic, single, non-faecal, 9=organic, single, non-faecal that is not dead phyto (fragment)
				  initPft,     & ! 1=diat, 2=flagel, 3=cocco, 4=pico
				  phase,       & ! 1=in water column, 2=at seafloor, 3=emptied/beyond detection limit/< 0.2 um diameter
				  living,      & ! 1=living phytoplankton(has not aggregated), 0=dead phytoplankton/TEP/terrigenous/zoo/faecal/agg
				  faecal,      & ! 1=faecal pellet, 0=non-faecal particle
				  tstepCreat     ! time step of creation

real(C_DOUBLE) :: nPxC,  & 		 ! no. particles per cluster / SLAMS-1.0: n (code), n (paper)
				  nPpxP, & 		 ! np. primary particles per aggregate particle / SLAMS-1.0: Nn (code), p (paper)
				  nPpxC    		 ! no. primary particles per cluster / SLAMS-1.0: p (code) (=n x Nn) 	

real(C_DOUBLE) 		         :: molesOrgC,     & ! mol organic C per particle
		  				        molesTepC,     & ! mol TEP-C per particle
          				        massOrgMatter, & ! mass organic matter (=prot+carb+lip+nucleics), g per particle
		  				        massTep          ! mass TEP-C, g per particle					 						 
real(C_DOUBLE), dimension(3) :: molesMineral,  & ! 1=opal(diat+radiol), 2=calcite(cocco), 3=clay
						        massMineral   

real(C_DOUBLE) :: depth,         & ! cluster depth, m
				  mass,          & ! particle mass, g
				  stickiness,    & ! particle stickness, unitless (0-1)
				  fracDim,       & ! particle fractal dimension, unitless (1-3)			  
				  radiusPp,      & ! radius of the primary particles of the particle, um
				  radius,        & ! radius of the particle, um 				  
				  solidVolume, 	 & ! material volume of the particle, um3
				  porosity,      & ! particle porosity, unitless (0-1)
				  density,       & ! particle density, g cm-3 
				  excessDensity, & ! particle excess density (density particle - density water), g cm-3 
				  velocity         ! particle settling velocity, m d-1
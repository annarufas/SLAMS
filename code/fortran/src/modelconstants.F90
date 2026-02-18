module modelconstants

! ----------------------------------------------------------------------------------------
! This section defines the universal constants used in SLAMS. Model constants are
! written in block capitals and grouped by their type.
! ----------------------------------------------------------------------------------------

implicit none

public

! Material physical constants
real*8, parameter :: MOLAR_MASS_CARBON   = 12.01d0, & ! g mol-1
				     MOLAR_MASS_CACO3    = 100.1d0, & ! g mol-1
				   	 MOLAR_MASS_OPAL     = 67.3d0,  & ! g mol-1, hydrated amorphous silica (SiO2·0.4(H2O)), after Mortlock & Froelich (1989), use a factor of 2.4 for biogenic silica/silica
				     MOLAR_MASS_CLAY     = 389.34d0,& ! g mol-1, illite ((K,H3O)(Al,Mg,Fe)2(Si,Al)4O10[(OH)2,(H2O)]) (Journet et al. 2008 indicate in Table 1 it's the most abundant type of clay)         					 
					 RHO_ORGMATTER       = 1.06d0,  & ! g mol-1, marine organic matter (proteins + carbohydrates + lipids + nucleic acids) 
					 RHO_TEP             = 0.800d0, & ! g cm-3, after Mari et al. 2017, Azetsu-Scott & Passow 2004
					 RHO_CALCITE         = 2.72d0,  & ! g cm-3
					 RHO_OPAL            = 2.0d0,   &  ! g cm-3, hydrated amorphous silica (SiO2·nH2O), after DeMaster, the density of biogenic silica equals 2.0 g cm-3 for diatoms and sponges, whereas radiolarian densities range from 1.7-2.0 g cm-3. The water content of biogenic silica varies from 8 wt.% to 17 wt.%, depending on the type of siliceous biota as well as their age (Hurd and Theyer, 1977).
					 RHO_CLAY            = 2.70d0,  & ! g cm-3 
					 MOLAR_VOLUME_OXYGEN = 22.4d0,  & ! L mol-1	
					 MESOZOO_DENSITY     = 1.0d0,  &  ! g cm-3, after Kiorboe 2013	
					 
! Physical constants
 
				     GRAVITY_CNT         = 9.81d0,   & ! m s-2
				     BOLTZMANN_CNT       = 1.38d-23, & ! J K-1 (= m2 kg s-2 K-1)
				     AVOGADRO_CNT 	     = 6.02d23,  & ! photons, 1 mol light = 6.02d23 mol photons
				     ENERGY_PHOTON       = 3.90d-19, & ! J per photon	
				     REYNOLDS_LAMINAR_LIMIT = 1d-1,  & ! Re below which Stokes approx used	
			     	
! Conversion factors

				     PI              = 3.14159265358979323846d0, &
				     DEG_TO_RAD      = PI/180d0, & ! from degrees to radians
				     RAD_TO_DEG      = 180d0/PI, & ! from radians to degrees
				     WATT_TO_PHOTON  = 1d6/(ENERGY_PHOTON*AVOGADRO_CNT), & ! from W m-2 to umol photons m-2 s-1
				     SECONDS_PER_DAY = 24d0*3600d0
				     				     				   
end module modelconstants

! Unused:
!
!MOLAR_MASS_SILICON  = 28.09d0 
!MOLAR_MASS_OXYGEN   = 15.99d0
!MOLAR_MASS_HYDROGEN = 1.008d0	
!RHO_CARBON          = 2.27d0
!RHO_ARAGONITE       = 2.90d0	
!			     
!ENERGY_DISSIPATION_RATE          = 1.00d-8  ! m2 s-3
!K_EPPLEY                         = 0.0633d0 ! (ºC)-1 
!PHOTOSYNTHESIS_ACTIVATION_ENERGY = 0.32d0   ! eV (methabolic theory of ecology)
!SCATTERING_COEFF_SEAWATER        = 0.05d0,   & ! m-1 (Kirk 1994),0.1 m-1 Acevedo-Trejos 2016 (PSFD model) (https://oceancolor.gsfc.nasa.gov/atbd/kd_490/)
!ABSORPTION_COEFF_CHLA            = 1.5d-2,   & ! chl a-specific light-absorption coefficient, m2 (mg Chla)-1 (Anderson et al. 2015) 
!MU0                              = 5.165d-6, & ! g cm-1 s-1, water viscosity parameter
!VISCO_EXPONENT                   = 2240d0,   & ! exponent of the curve that relates water temperature with viscosity, mu = mu0*exp(b/tempKelvin)
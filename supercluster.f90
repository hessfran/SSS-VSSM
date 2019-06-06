program toy_cluster_contraction
	implicit none
	!-------------------------------------------------------------------------------------------------------------------------------
	!
	! This program is used for testing the performance of different Supercluster Contraction models.
	! It assumes a rectangular lattice (C2v symmetry) such as the RuO2(110) surface with four species + vacant site (five in total).
	! Add your own energy functions to find the optimum Supercluster Contraction.
	! If you use this program in testing or optimizing Supercluster Contractions, please cite the parent paper.
	! Author: F. Hess (dr.franziska.hess@gmail.com)
	! Feel free to contact me if there is any problem and please notify me of any bugs you find.
	! This program was tested for ifort 11.1, ifort 18.0.1 and gfortran 7.1.0.
	! Compilation with gfortran requires the -ffree-line-length-200 flag.
	! No other flags should be necessary.
	!
	!-------------------------------------------------------------------------------------------------------------------------------
	! Instructions
	!-------------------------------------------------------------------------------------------------------------------------------
	! Adding a Supercluster Contraction requires (example)
	! - a subroutine for pre-calculation of the supercluster energies (precalculate_SCA_clusters) as well as
	! - a function for energy evaluation during production (energy_SCA), 
	! - arrays for storing the precalculated supercluster energies (SCA_I, SCA_II_1, SCA_II_2), and
	! - definitions of the Supercluster shapes (SCA_I_sites, SCA_II_1_sites, SCA_II_2_sites). For storing the results, 
	! - reals are required for storing timing information (time_pre_SCA, time_SCA), 
	! - a real array for storing the resulting energy (res_energy_SCA).
	! 
	! The computation of energies for random configurations is invoked via the subroutine
	! calculate_energies (real function energy_function, character(16) function_name, integer(:,:,:) allconfs, real elapsed_time, real(:) res_energy),
	! which takes the SC function as its first argument.
	! In the end the energies calculated with the Supercluster Contraction should be compared to the literal result 
	! to verify that the Supercluster Contraction gives the same result as the originial energy function.
	!
	!-------------------------------------------------------------------------------------------------------------------------------
	!
	integer::x,y,i     ! counters
	integer::n_int, R  ! other parameters
	real::   r_rand    ! a random number
	
	! Timing (add your own)
	real:: time_pre_SCA, time_pre_SCB, time_pre_SCC, time_pre_SCD                          ! pre-calculate times
	real:: time_literal_long, time_literal_short, time_SCA, time_SCB, time_SCC, time_SCD   ! actual calculation time
	real:: time_blind, time_ref                                                            ! reference times for determining relative speedup
	
	! "Cluster Expansion"
	real,    allocatable, dimension(:,:)		 :: interaction_base_values     ! just some number to generate a "Cluster Expansion"
	real,    allocatable, dimension(:)			 :: e0                          ! adsorption energy at zero coverage
	real,    allocatable, dimension(:,:,:,:)	 :: pairwise_interactions       ! Pairwise interaction parameters
	real,    allocatable, dimension(:,:,:)		 :: Lshape_1,Lshape_2, Lshape_3 ! Three-body interaction parameters -- not actually assigned in the present program, used only for simulating the time for literal evaluation
	integer, allocatable, dimension(:,:)		 :: interacting_sites           ! A tweak to speed up literal evaluation on some architectures, cf energy_literal_long vs. energy_literal_short
	integer, allocatable, dimension(:,:)		 :: i_rank_2_temporary          ! dummy
	integer, allocatable, dimension(:,:,:)       :: test_configuration          ! random configurations used for testing
	! energy results for cross-checking if the Supercluster Contraction gives correct results (add your own)
	real,    allocatable, dimension(:)			 :: res_energy_literal_long, res_energy_literal_short, res_energy_SCA, res_energy_SCB, res_energy_SCC, res_energy_SCD,res_energy_blind 
	
	!-------------------------------------------------------------------------------------------------------------------------------
	!
	! Definition of cluster contractions
	! (add your own)
	!
	!-------------------------------------------------------------------------------------------------------------------------------
	! Instructions
	!-------------------------------------------------------------------------------------------------------------------------------
	! There are different ways to implement a Supercluster Contraction. In this implementation we store
	! - the cluster energies in a real array with rank equal to the number of sites in the supercluster. The first index gives the
	!   occupation of the central site, the following indices the occupations of the other sites with sequence the same as in the 
	!   Cluster shape definition.
	! - The cluster shape by defining the coordinates of the sites in an integer array with shape (1:2, 1: number of sites-1).
	!   The first dimension corresponds to the (x,y) positions relative to the center (0,0), the second dimension to the site.
	!   The center is defined as (0,0) and does not need to be given explicitly (hence number of sites-1). The sequence of sites
	!   in the cluster shape corresponds to the the dimensions in the cluster energy definition. 
	!   The reshape() function is used for compatibility with gfortran. 
	!   For defining cluster shapes, remember that in c2v symmetry top/down is the same as left/right, so - and + can be reversed.
	! 
	!-------------------------------------------------------------------------------------------------------------------------------
	
	! SCA
	real,    allocatable, dimension(:,:,:,:,:,:) :: SCA_I
	real,    allocatable, dimension(:,:,:,:)     :: SCA_II_1,SCA_II_2
	integer, dimension(1:2,1:5)::SCA_I_sites    = reshape((/(/-1,0/),(/-1,1/),(/0,1/),(/1,1/),(/1,0/)/),shape(sca_I_sites))
	integer, dimension(1:2,1:3)::SCA_II_1_sites = reshape((/(/-1, 1/),(/0,1/),(/0,2/)/),shape(sca_II_1_sites))
	integer, dimension(1:2,1:3)::SCA_II_2_sites = reshape((/(/ 1, 1/),(/1,0/),(/2,0/)/),shape(sca_II_2_sites))
	
	! SCB
	real,    allocatable, dimension(:,:,:,:,:,:) :: SCB_I
	real,    allocatable, dimension(:,:,:,:)     :: SCB_II
	integer, dimension(1:2,1:5)::SCB_I_sites  = reshape((/(/-1,0/),(/-1,1/),(/0,1/),(/1,1/),(/1,0/)/),shape(SCB_I_sites)) ! same as in SCA
	integer, dimension(1:2,1:3)::SCB_II_sites = reshape((/(/ 2, 0/),(/ 1, 0/),(/ 0, 2/)/),shape(SCB_II_sites))
	
	! SCC
	real,    allocatable, dimension(:,:,:,:,:,:) :: SCC_I
	real,    allocatable, dimension(:,:,:)       :: SCC_II
	integer, dimension(1:2,1:5)::SCC_I_sites  = reshape((/(/-1,0/),(/-1,1/),(/0,1/),(/1,1/),(/1,0/)/),shape(SCC_I_sites)) ! same as in SCA
	integer, dimension(1:2,1:2)::SCC_II_sites = reshape((/(/ 2, 0/),(/ 0, 2/)/),shape(SCC_II_sites))
	
	! SCD
	real,    allocatable, dimension(:,:,:,:) :: SCD_I
	real,    allocatable, dimension(:,:,:)       :: SCD_II
	integer, dimension(1:2,1:3)::SCD_I_sites  = reshape((/(/-1,0/),(/-1,1/),(/0,1/)/),shape(SCD_I_sites))
	integer, dimension(1:2,1:2)::SCD_II_sites = reshape((/(/ 2, 0/),(/ 0, 2/)/),shape(SCD_II_sites))
	
	! (add your own)
	! real,    allocatable, dimension(:,:,:,:) :: cluster_I
	! real,    allocatable, dimension(:,:,:)       :: cluster_II
	! integer, dimension(1:2,1:3)::cluster_I_sites  = reshape((/(/-1,0/),(/-1,1/),(/0,1/)/),shape(cluster_I_sites)) 
	! integer, dimension(1:2,1:2)::cluster_II_sites = reshape((/(/ 2, 0/),(/ 0, 2/)/),shape(cluster_II_sites))
	
	!-------------------------------------------------------------------------------------------------------------------------------
	
	!-------------------------------------------------------------------------------------------------------------------------------
	! Other things you can play with
	
	integer:: n_species  = 4   ! number of species
	integer:: verbosity  = 0   ! how much output to write ( 0 = only main results, 1 = additional output, 2 = all output)
	integer:: n_rnd_conf = 10000000 ! number of random configurations for timing
	
	!-------------------------------------------------------------------------------------------------------------------------------
	! defining some lateral interactions 
	! base values used for generating an interaction parameter set (this is not a real cluster expansion, just some numbers)
	! define only upper triangular matrix due to symmetry (interaction of 1 with 2 is equal to interaction of 2 with 1)
	! !!! IF YOU CHANGE n_species YOU NEED TO CHANGE THIS AS WELL !!!
	!-------------------------------------------------------------------------------------------------------------------------------
	allocate (e0(0:n_species))
	e0 = (/0.0, -1.0, -1.1, -1.2, -1.3/)
	allocate (interaction_base_values	(0:n_species, 0:n_species))
	interaction_base_values=0
	interaction_base_values(1,1:n_species)=(/0.10, 0.05, 0.08, 0.09/)
	interaction_base_values(2,2:n_species)=(/0.02, 0.03, 0.04/)
	interaction_base_values(3,3:n_species)=(/0.11, 0.12/)
	interaction_base_values(4,4:n_species)=(/0.07/)
	! 
	!--------------------------------------------------------------------------------------------------------------------------------
	
	! definition of other constants
	R = 2 ! Range of lateral interactions, can be changed in principle, but supercluster definitions are hardcoded
	
	! Initialize a few things 
	! no changes should be necessary here, unless you want to test a different Cluster Expansion
	call random_seed() ! initialize RNG
	write(*,*)'Initializing data structures'
	call initialize_data_structures
	call generate_interaction_parameter_set(n_species, R, interaction_base_values, pairwise_interactions)
	call identify_interacting_sites(R,n_int)
	write(*,*)'Generating random configurations'
	call generate_random_configurations (R,n_rnd_conf)
	
	! Precalculation timing
	! (add your own)
	write(*,*)'Precalculating superclusters'
	call precalculate_SCA_clusters(time_pre_SCA)
	call precalculate_SCB_clusters(time_pre_SCB)
	call precalculate_SCC_clusters(time_pre_SCC)
	call precalculate_SCD_clusters(time_pre_SCD)
	
	! actual timing starts here
	! (add your own)
	write(*,*)'Timing Hamiltonians'
	call calculate_energies(energy_literal_long, 'literal long    ', test_configuration, time_literal_long,  res_energy_literal_long )
	call calculate_energies(energy_literal_short,'literal short   ', test_configuration, time_literal_short, res_energy_literal_short)
	call calculate_energies(energy_sca,          'SCA             ', test_configuration, time_SCA,           res_energy_SCA)
	call calculate_energies(energy_scb,          'SCB             ', test_configuration, time_SCB,           res_energy_SCB)
	call calculate_energies(energy_scc,          'SCC             ', test_configuration, time_SCC,           res_energy_SCC)
	call calculate_energies(energy_scd,          'SCC             ', test_configuration, time_SCD,           res_energy_SCD)
	call calculate_energies(energy_blind,        'blind           ', test_configuration, time_blind,         res_energy_blind)
	
	! Compute reference time
	time_ref = minval((/time_literal_short, time_literal_long/))-time_blind
	
	! Main computation results
	! (add your own)
	write(*,*)
	write(*,*)'Pre-compute times'
	write(*,*)'SCA           ', time_pre_SCA
	write(*,*)'SCB           ', time_pre_SCB
	write(*,*)'SCC           ', time_pre_SCC
	write(*,*)'SCD           ', time_pre_SCD
	write(*,*)
	
	! (add your own)
	write(*,*)'absolute time and relative speed-up for energy evaluation'
	write(*,*)'literal long  ', time_literal_long,  time_ref/(time_literal_long -time_blind) 
	write(*,*)'literal short ', time_literal_short, time_ref/(time_literal_short-time_blind)
	write(*,*)'SCA           ', time_SCA,           time_ref/(time_SCA          -time_blind) 
	write(*,*)'SCB           ', time_SCB,           time_ref/(time_SCB          -time_blind) 
	write(*,*)'SCC           ', time_SCC,           time_ref/(time_SCC          -time_blind) 
	write(*,*)'SCD           ', time_SCD,           time_ref/(time_SCD          -time_blind) 
	write(*,*)'blind         ', time_blind
	write(*,*)
	
	! (add your own)
	write(*,*)'Average energy errors'
	write(*,*)'(should be only floating-point error, < 1E-6)'
	write(*,*)'short',sum(abs(res_energy_literal_short)-abs(res_energy_literal_long))/n_rnd_conf
	write(*,*)'SCA  ',sum(abs(res_energy_SCA)-abs(res_energy_literal_long))/n_rnd_conf
	write(*,*)'SCB  ',sum(abs(res_energy_SCB)-abs(res_energy_literal_long))/n_rnd_conf
	write(*,*)'SCC  ',sum(abs(res_energy_SCC)-abs(res_energy_literal_long))/n_rnd_conf
	write(*,*)'SCD  ',sum(abs(res_energy_SCD)-abs(res_energy_literal_long))/n_rnd_conf
	write(*,*)
	
	
	stop
	contains
	
	! ----------------------------------------------------------------------------------------------------------------------------------
	!
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	! Energy computation
	! No changes should be necessary here
	
	subroutine calculate_energies (energy_function,function_name,allconfs,elapsed_time, res_energy)
		implicit none
		integer,dimension(1:n_rnd_conf,-R:R,-R:R),intent(in)::allconfs ! avoid accessing global variables in the loop below
		real                                      :: energy_function ! argument function
		character(len=16),            intent(in)  ::function_name
		real,                         intent(out) :: elapsed_time
		real,dimension(1:n_rnd_conf), intent(out) ::res_energy

		integer :: c ! counter
		real    :: time_start, time_end
		
		! actual energy computation
		call cpu_time(time_start)
		do c=1, n_rnd_conf
			res_energy(c)=energy_function(allconfs(c,-R:R,-R:R))
		end do
		call cpu_time(time_end)
		elapsed_time = time_end - time_start
		
		! output energies of configurations
		if(verbosity >= 2) then
			write(*,*)function_name, ' energy calculation results'
			do c = 1, n_rnd_conf
				write(*,'(i4,f14.8)') c,res_energy(c)
			end do
			write(*,*)'time: ', elapsed_time
			write(*,*)
		end if
		
	end subroutine calculate_energies
	
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	!
	! Actual energy functions
	! (add your own)
	!
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_literal_long(configuration)
		implicit none
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer :: sc, x, y
		real    :: e
		
		! This is a possible, but sometimes not optimal algorithm for energy computation (again this depends on architecture)
		sc = configuration(0,0)
		e  = e0(sc) ! adsorption energy at zero coverage
		do x = -R,R
			do y = -R,R
				e = e + pairwise_interactions( sc, configuration(x,y), abs(x), abs(y) )
			end do
		end do
		
		! Three-body interactions not actually considered, only simulated for timing
		! 8 terms for L-shaped interactions, 
		! untested feature, WILL give wrong results
		e = e + Lshape_1(sc,configuration(-1,0),configuration(-1, 1)) 
		e = e + Lshape_1(sc,configuration( 1,0),configuration( 1, 1)) 
		e = e + Lshape_1(sc,configuration(-1,0),configuration(-1,-1)) 
		e = e + Lshape_1(sc,configuration( 1,0),configuration(-1,-1)) 
		e = e + Lshape_2(sc,configuration(0,-1),configuration(-1, 1)) 
		e = e + Lshape_2(sc,configuration(0, 1),configuration( 1, 1)) 
		e = e + Lshape_2(sc,configuration(0,-1),configuration(-1,-1)) 
		e = e + Lshape_2(sc,configuration(0, 1),configuration(-1,-1)) 
		e = e + Lshape_3(sc,configuration(0, 1),configuration( 1, 0))
		e = e + Lshape_3(sc,configuration(0,-1),configuration( 1, 0))
		e = e + Lshape_3(sc,configuration(0, 1),configuration(-1, 0))
		e = e + Lshape_3(sc,configuration(0,-1),configuration(-1, 0))
		
		energy_literal_long=e
	end function energy_literal_long
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_literal_short(configuration)
		implicit none
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer :: sc, x, y, i
		real    :: e
		
		! This algorithm is a little better than energy_literal_long, but does the same thing in principle
		sc = configuration(0,0)
		e  = e0(sc) ! adsorption energy at zero coverage
		do i=1,n_int
			x=interacting_sites(i,1)
			y=interacting_sites(i,2)
			e = e + pairwise_interactions( sc, configuration(x,y), abs(x), abs(y) )
		end do
		
		! Three-body interactions not actually considered, only simulated for timing
		! 8 terms for L-shaped interactions, 
		! untested feature, WILL give wrong results
		e = e + Lshape_1(sc,configuration(-1,0),configuration(-1, 1)) 
		e = e + Lshape_1(sc,configuration( 1,0),configuration( 1, 1)) 
		e = e + Lshape_1(sc,configuration(-1,0),configuration(-1,-1)) 
		e = e + Lshape_1(sc,configuration( 1,0),configuration(-1,-1)) 
		e = e + Lshape_2(sc,configuration(0,-1),configuration(-1, 1)) 
		e = e + Lshape_2(sc,configuration(0, 1),configuration( 1, 1)) 
		e = e + Lshape_2(sc,configuration(0,-1),configuration(-1,-1)) 
		e = e + Lshape_2(sc,configuration(0, 1),configuration(-1,-1)) 
		e = e + Lshape_3(sc,configuration(0, 1),configuration( 1, 0))
		e = e + Lshape_3(sc,configuration(0,-1),configuration( 1, 0))
		e = e + Lshape_3(sc,configuration(0, 1),configuration(-1, 0))
		e = e + Lshape_3(sc,configuration(0,-1),configuration(-1, 0))
		
		energy_literal_short=e
	end function energy_literal_short
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_SCA(configuration)
		implicit none
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer :: sc
		real    :: e
		
		! uses three different clusters
		
		! cluster shape definition 
		! integer, dimension(1:2,1:5)::SCA_I_sites=(/(/-1,0/),(/-1,1/),(/0,1/),(/1,1/),(/1,0/)/)
		! integer, dimension(1:2,1:3)::SCA_II_1_sites=(/(/-1, 1/),(/0,1/),(/0,2/)/)
		! integer, dimension(1:2,1:3)::SCA_II_2_sites=(/(/ 1, 1/),(/1,0/),(/2,0/)/)
		
		sc = configuration(0,0) ! occupancy of center site
		e  = e0(sc)             ! adsorption energy at zero coverage
		
		! Cluster shape definitions are hardcoded here for improving performance in production, but cluster shape definitions 
		! can also be used. 
		
		! SCA_I 
		e = e + sca_I(sc, configuration(-1,0),configuration(-1, 1),configuration(0, 1),configuration(1, 1),configuration(1,0))
		! y reversed (c2v symmetry)
		e = e + sca_I(sc, configuration(-1,0),configuration(-1,-1),configuration(0,-1),configuration(1,-1),configuration(1,0))
		
		! SCA_II_1
		e = e + sca_II_1(sc, configuration(-1, 1),configuration(0, 1),configuration(0, 2))
		! x and y reversed (c2v symmetry)
		e = e + sca_II_1(sc, configuration( 1,-1),configuration(0,-1),configuration(0,-2))
		
		! SCA_II_2
		e          = e + sca_II_2(sc, configuration( 1, 1),configuration( 1, 0),configuration( 2, 0))
		! x and y reversed (c2v symmetry)
		energy_SCA = e + sca_II_2(sc, configuration(-1,-1),configuration(-1, 0),configuration(-2, 0))
		
	end function energy_SCA
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_SCB(configuration)
		implicit none
		
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer :: sc
		real    :: e
		
		! uses two different clusters
		
		! cluster shape definition 
		! integer, dimension(1:2,1:5)::SCB_I_sites  = (/(/-1,0/),(/-1,1/),(/0,1/),(/1,1/),(/1,0/)/) ! same as in SCA
		! integer, dimension(1:2,1:4)::SCB_II_sites = (/(/ 2, 0/),(/ 1, 0/),(/ 0, 2/)/)
		
		sc = configuration(0,0)
		
		! SCB_I 
		e =     scb_I(sc, configuration(-1,0),configuration(-1, 1),configuration(0, 1),configuration(1, 1),configuration(1,0))
		e = e + scb_I(sc, configuration(-1,0),configuration(-1,-1),configuration(0,-1),configuration(1,-1),configuration(1,0))
		
		! SCB_II
		e = e + scb_II(sc, configuration( 2, 0),configuration( 1, 0),configuration( 0, 2))
		e = e + scb_II(sc, configuration(-2, 0),configuration(-1, 0),configuration( 0,-2))
		
		energy_SCB=e
		
	end function energy_SCB
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_SCC(configuration)
		implicit none
		
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer :: sc
		real    :: e
		
		! uses two different clusters of smaller size
		
		! cluster shape definition 
		! integer, dimension(1:2,1:5)::SCC_I_sites  = (/(/-1,0/),(/-1,1/),(/0,1/),(/1,1/),(/1,0/)/) ! same as in SCA
		! integer, dimension(1:2,1:2)::SCC_II_sites = (/(/ 2, 0/),(/ 0, 2/)/)
		
		sc = configuration(0,0)
		
		! SCC_I 
		e =     SCC_I(sc, configuration(-1,0),configuration(-1, 1),configuration(0, 1),configuration(1, 1),configuration(1,0))
		e = e + SCC_I(sc, configuration(-1,0),configuration(-1,-1),configuration(0,-1),configuration(1,-1),configuration(1,0))
		
		! SCC_II
		e = e + SCC_II(sc, configuration( 2, 0),configuration( 0, 2))
		e = e + SCC_II(sc, configuration(-2, 0),configuration( 0,-2))
		
		energy_SCC=e
		
	end function energy_SCC
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_SCD(configuration)
		implicit none
		
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer :: sc
		real    :: e
		
		! makes more efficient use of symmetry by using two different clusters of even smaller size, but increases number of terms
		
		! cluster shape definition 
		! integer, dimension(1:2,1:3)::SCD_I_sites  = (/(/-1,0/),(/-1,1/),(/0,1/)/) ! same as in SCA
		! integer, dimension(1:2,1:2)::SCD_II_sites = (/(/ 2, 0/),(/ 0, 2/)/)
		
		sc = configuration(0,0)
		
		! SCD_I 
		e =     SCD_I(sc, configuration(-1,0),configuration(-1, 1),configuration(0, 1))
		e = e + SCD_I(sc, configuration(-1,0),configuration(-1,-1),configuration(0,-1))
		e = e + SCD_I(sc, configuration( 1,0),configuration( 1, 1),configuration(0, 1))
		e = e + SCD_I(sc, configuration( 1,0),configuration( 1,-1),configuration(0,-1))
		
		! SCD_II
		e = e + SCD_II(sc, configuration( 2, 0),configuration( 0, 2))
		e = e + SCD_II(sc, configuration(-2, 0),configuration( 0,-2))
		
		energy_SCD=e
		
	end function energy_SCD
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_blind(configuration)
		implicit none
		
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer:: sc
		real:: e
		
		! used only for determining time reference
		
		energy_blind=0
		
	end function energy_blind
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function energy_only_pairwise(configuration)
		implicit none
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		integer:: sc, x, y, i
		real:: e
		
		! Takes only pairwise interactions, used in the precalculation of Supercluster energies, but otherwise similar to energy_literal_short
		
		sc = configuration(0,0)
		e = 0
		do i=1,n_int
			x=interacting_sites(i,1)
			y=interacting_sites(i,2)
			e = e + pairwise_interactions( sc, configuration(x,y), abs(x), abs(y) )
		end do
		energy_only_pairwise=e
	end function energy_only_pairwise
	
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	!
	! Supercluster precalculation functions
	! (add your own)
	! 
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine precalculate_SCA_clusters(elapsed_time)
		implicit none
		real:: time_start, time_end
		real, intent(out):: elapsed_time
		integer,dimension(-R:R,-R:R):: sca_conf
		integer:: s1, s2, s3, s4
		
		call cpu_time(time_start)
		allocate(res_energy_SCA (1:n_rnd_conf))
		allocate(SCA_I    (0:n_species,0:n_species,0:n_species,0:n_species,0:n_species,0:n_species))
		allocate(SCA_II_1 (0:n_species,0:n_species,0:n_species,0:n_species))
		allocate(SCA_II_2 (0:n_species,0:n_species,0:n_species,0:n_species))
		
		
		! SCA_I is actually 0 because three-body interactions are not considered...
		SCA_I = 0
		
		! SCA_II_1
		sca_conf=0
		do s1 = 1,n_species
			sca_conf(0,0) = s1
			do s2 = 0,n_species
				sca_conf(sca_II_1_sites(1,1),sca_II_1_sites(2,1)) = s2
				do s3 = 0,n_species
					sca_conf(sca_II_1_sites(1,2),sca_II_1_sites(2,2)) = s3
					do s4 = 0,n_species
						sca_conf(sca_II_1_sites(1,3),sca_II_1_sites(2,3)) = s4
						sca_II_1(s1,s2,s3,s4) = energy_only_pairwise(sca_conf)
						!write(*,'(5(i2,1x))')(sca_conf(:,y),y=-2,2)
						!write(*,*)
						!read(*,*)
					end do
				end do
			end do
		end do
		
		! SCA_II_2
		sca_conf=0
		do s1 = 1,n_species
			sca_conf(0,0) = s1
			do s2 = 0,n_species
				sca_conf(sca_II_2_sites(1,1),sca_II_2_sites(2,1)) = s2
				do s3 = 0,n_species
					sca_conf(sca_II_2_sites(1,2),sca_II_2_sites(2,2)) = s3
					do s4 = 0,n_species
						sca_conf(sca_II_2_sites(1,3),sca_II_2_sites(2,3)) = s4
						sca_II_2(s1,s2,s3,s4) = energy_only_pairwise(sca_conf)
					end do
				end do
			end do
		end do
		call cpu_time(time_end)
		elapsed_time=time_end-time_start
		
	end subroutine precalculate_SCA_clusters
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine precalculate_SCB_clusters(elapsed_time)
		implicit none
		real:: time_start, time_end
		real, intent(out):: elapsed_time
		integer,dimension(-R:R,-R:R):: scb_conf
		integer:: s1, s2, s3, s4, s5, s6
		
		call cpu_time(time_start)
		allocate(res_energy_SCB (1:n_rnd_conf))
		allocate(SCB_I  (0:n_species,0:n_species,0:n_species,0:n_species,0:n_species,0:n_species))
		allocate(SCB_II (0:n_species,0:n_species,0:n_species,0:n_species))
				
		! again three-body interactions are not considered...
		scb_conf=0
		scb_I=0
		do s1 = 1, n_species
			scb_conf(0,0) = s1
			do s2 = 0, n_species
				scb_conf(scb_I_sites(1,1),scb_I_sites(2,1)) = s2
				do s3 = 0, n_species
					scb_conf(scb_I_sites(1,2),scb_I_sites(2,2)) = s3
					do s4 = 0, n_species
						scb_conf(scb_I_sites(1,3),scb_I_sites(2,3)) = s4
						do s5 = 0, n_species
							scb_conf(scb_I_sites(1,4),scb_I_sites(2,4)) = s5
							do s6 = 0, n_species
								scb_conf(scb_I_sites(1,5),scb_I_sites(2,5)) = s6
								scb_I(s1,s2,s3,s4,s5,s6) = energy_only_pairwise(scb_conf)-(pairwise_interactions(s1,scb_conf(-1,0),1,0)+pairwise_interactions(s1,scb_conf(1,0),1,0)) ! remove these 1NN interactions to prevent double counting with SCB_II		
							end do
						end do
					end do
				end do
			end do
		end do
		
		! SCB_II
		scb_II=0
		scb_conf=0
		do s1 = 1,n_species
			scb_conf(0,0) = s1
			do s2 = 0,n_species
				scb_conf(scb_II_sites(1,1),scb_II_sites(2,1)) = s2
				do s3 = 0,n_species
					scb_conf(scb_II_sites(1,2),scb_II_sites(2,2)) = s3
					do s4 = 0,n_species
						scb_conf(scb_II_sites(1,3),scb_II_sites(2,3)) = s4
						scb_II(s1,s2,s3,s4) = energy_only_pairwise(scb_conf)+e0(s1)/2.
					end do
				end do
			end do
		end do
		
		call cpu_time(time_end)
		elapsed_time=time_end-time_start
		
	end subroutine precalculate_SCB_clusters
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine precalculate_SCC_clusters(elapsed_time)
		implicit none
		real:: time_start, time_end
		real, intent(out):: elapsed_time
		integer,dimension(-R:R,-R:R):: scc_conf 
		integer:: s1, s2, s3, s4, s5, s6
		
		call cpu_time(time_start)
		allocate(res_energy_SCC (1:n_rnd_conf))
		allocate(scc_I  (0:n_species,0:n_species,0:n_species,0:n_species,0:n_species,0:n_species))
		allocate(scc_II (0:n_species,0:n_species,0:n_species))
				
		! again  three-body interactions are not considered...
		scc_conf=0
		scc_I=0
		do s1 = 1, n_species
			scc_conf(0,0) = s1
			do s2 = 0, n_species
				SCC_conf(SCC_I_sites(1,1),SCC_I_sites(2,1)) = s2
				do s3 = 0, n_species
					SCC_conf(SCC_I_sites(1,2),SCC_I_sites(2,2)) = s3
					do s4 = 0, n_species
						SCC_conf(SCC_I_sites(1,3),SCC_I_sites(2,3)) = s4
						do s5 = 0, n_species
							SCC_conf(SCC_I_sites(1,4),SCC_I_sites(2,4)) = s5
							do s6 = 0, n_species
								SCC_conf(SCC_I_sites(1,5),SCC_I_sites(2,5)) = s6
								SCC_I(s1,s2,s3,s4,s5,s6) = energy_only_pairwise(SCC_conf)-(pairwise_interactions(s1,SCC_conf(-1,0),1,0)+pairwise_interactions(s1,SCC_conf(1,0),1,0))/2 ! split this 1NN interaction between the two SCD_I clusters
							end do
						end do
					end do
				end do
			end do
		end do
		
		! SCC_II
		SCC_II=0
		SCC_conf=0
		do s1 = 1,n_species
			SCC_conf(0,0) = s1
			do s2 = 0,n_species
				SCC_conf(SCC_II_sites(1,1),SCC_II_sites(2,1)) = s2
				do s3 = 0,n_species
					SCC_conf(SCC_II_sites(1,2),SCC_II_sites(2,2)) = s3
					SCC_II(s1,s2,s3) = energy_only_pairwise(SCC_conf)+e0(s1)/2.					
				end do
			end do
		end do
		
		call cpu_time(time_end)
		elapsed_time=time_end-time_start
		
	end subroutine precalculate_SCC_clusters
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine precalculate_SCD_clusters(elapsed_time)
		implicit none
		real:: time_start, time_end
		real, intent(out):: elapsed_time
		integer,dimension(-R:R,-R:R):: SCD_conf 
		integer:: s1, s2, s3, s4, s5, s6
		
		call cpu_time(time_start)
		allocate(res_energy_SCD (1:n_rnd_conf))
		allocate(SCD_I (0:n_species,0:n_species,0:n_species,0:n_species))
		allocate(SCD_II(0:n_species,0:n_species,0:n_species))
				
		! again three-body interactions are not considered...
		SCD_conf=0
		SCD_I=0
		do s1 = 1, n_species
			SCD_conf(0,0) = s1
			do s2 = 0, n_species
				SCD_conf(SCD_I_sites(1,1),SCD_I_sites(2,1)) = s2
				do s3 = 0, n_species
					SCD_conf(SCD_I_sites(1,2),SCD_I_sites(2,2)) = s3
					do s4 = 0, n_species
						SCD_conf(SCD_I_sites(1,3),SCD_I_sites(2,3)) = s4
						SCD_I(s1,s2,s3,s4) = energy_only_pairwise(SCD_conf)- (pairwise_interactions(s1,SCD_conf(-1,0),1,0)+pairwise_interactions(s1,SCD_conf(0,1),0,1))/2 ! split this 1NN interaction between the two SCD_I clusters
					end do
				end do
			end do
		end do
		
		! SCD_II
		SCD_II=0
		SCD_conf=0
		do s1 = 1,n_species
			SCD_conf(0,0) = s1
			do s2 = 0,n_species
				SCD_conf(SCD_II_sites(1,1),SCD_II_sites(2,1)) = s2
				do s3 = 0,n_species
					SCD_conf(SCD_II_sites(1,2),SCD_II_sites(2,2)) = s3
					SCD_II(s1,s2,s3) = energy_only_pairwise(SCD_conf)+e0(s1)/2.	! split e0 between the two SCD_II clusters				
				end do
			end do
		end do
		
		call cpu_time(time_end)
		elapsed_time=time_end-time_start
		
	end subroutine precalculate_SCD_clusters
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	!
	! Auxiliary functions
	! No changes should be necessary here
	!
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine generate_interaction_parameter_set(n_species, R, interaction_base_values, pairwise_interactions)
		implicit none
		integer,intent(in)												:: n_species, R
		real, dimension(0:n_species,0:n_species), intent(inout)			:: interaction_base_values
		real, dimension(0:n_species,0:n_species,0:R,0:R), intent(inout)	:: pairwise_interactions
		integer:: s1, s2
		integer :: x, y
		
		! initialize
		pairwise_interactions = 0
		! L-shaped 3-body interactions are not assigned, just added for timing purposes. 
		Lshape_1 = 0 ! not assigned
		Lshape_2 = 0 ! not assigned
		Lshape_3 = 0 ! not assigned
		
		! symmetrize interactions
		do s1=1,4
			do s2=1,s1-1
				interaction_base_values(s1,s2) = interaction_base_values(s2,s1)
			end do
		end do
		
		! interaction base values are weighted by distance to center squared (1/R^2)
		do s1=1,n_species
			do s2=1,n_species
				do x=0,R
					do y=0,R
						if( (x==0 .and. y==0) .or. x*x + y*y > R*R) then
							! center and sites outside interaction range are skipped
							cycle
						else
							! Interactions decay as 1/R^2
							! Substrate symmetry is c2v => x and y units can be different -------vvv ------vvv 
							pairwise_interactions(s1,s2,x,y) = interaction_base_values(s1, s2)/ (1.0*x*x + 1.2*y*y)
						end if
					end do
				end do
			end do
		end do
		
		! output of interaction parameters
		! pay attention to symmetry of lateral interactions
		if(verbosity >= 1)then
			do s1 = 1, n_species
				do s2 = 1, n_species
					write(*,*) 'interaction of species ', s1, ' with ', s2
					write(*,'(3(f6.4,1x))')(pairwise_interactions(s1,s2,:,y),y=0,R)
					write(*,*)
				end do
			end do
		end if
	end subroutine generate_interaction_parameter_set
	
	! ----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine identify_interacting_sites(R,n_int)
		implicit none
		integer,intent(inout):: R, n_int
		integer:: i, x, y
		
		interacting_sites=0 !initialize
		
		! Determine sites within radius of lateral interactions
		i=0
		do x = -R,+R
			do y = -R,+R
				if(x*x+y*y <= R*R)then
					i=i+1
					interacting_sites(i,1:2)=(/x,y/)
				end if
			end do
		end do
		n_int=i ! total nuber of interacting sites
		
		! reallocate interacting sites array using array temporary
		allocate   (i_rank_2_temporary(1:n_int,1:2))
		i_rank_2_temporary = interacting_sites(1:n_int, 1:2)
		deallocate (interacting_sites)
		allocate   (interacting_sites(1:n_int,1:2))
		interacting_sites = i_rank_2_temporary
		deallocate (i_rank_2_temporary)
		
		!output result
		if(verbosity >=1) then
			write(*,"(2(i5,1x))")(interacting_sites(i,1:2),i=1,n_int)
		end if
		
	end subroutine identify_interacting_sites
	
	!----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine initialize_data_structures
		implicit none
		! allocate arrays
		! Interacting sites within radius R, only 0:R due to c2v substrate symmetry
		allocate (pairwise_interactions		(0:n_species, 0:n_species, 0:R, 0:R)) ! define only one quarter of interaction matrix for c2v substrate symmetry
		allocate (Lshape_1(0:n_species,0:n_species,0:n_species))
		allocate (Lshape_2(0:n_species,0:n_species,0:n_species))
		allocate (Lshape_3(0:n_species,0:n_species,0:n_species))
		allocate (interacting_sites(1:4*R*R, 1:2)) ! temporary allocation, will be changed later
		allocate (test_configuration(1:n_rnd_conf,-R:R,-R:R))
		allocate (res_energy_literal_long  (1:n_rnd_conf))
		allocate (res_energy_literal_short (1:n_rnd_conf))
		allocate (res_energy_blind         (1:n_rnd_conf))
		
	end subroutine initialize_data_structures
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine generate_random_configurations (R, n_rnd_conf)
		implicit none
		integer,intent(in)::R, n_rnd_conf
		integer::c
		real,dimension(-R:R,-R:R)::tempconf
		
		do c = 1,n_rnd_conf
			call random_number (tempconf)
			test_configuration(c,-R:R,-R:R)=int(tempconf*5) ! equal probability of occupancies for species 0 to 4
		end do
		
		if(verbosity>=2)then
			do c = 1, n_rnd_conf
			write(*,*)'configuration ',c
			write(*,'(5(i2,1x))') (test_configuration(c,:,y),y=-R,R)
			write(*,*)
			end do
		end if
		
	end subroutine generate_random_configurations	
	
	
	
	
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	
	
	
end program toy_cluster_contraction
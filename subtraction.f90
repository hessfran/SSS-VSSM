program subtraction_scheme
	implicit none
	!-------------------------------------------------------------------------------------------------------------------------------
	!
	! This program is used for comparing the performance of approaches to recalculate Gamma. 
	! It assumes a rectangular lattice (C2v symmetry) such as the RuO2(110) surface with four species + vacant site (five in total).
	! You can exchange/add different Supercluster Contractions to determine the combined speedup.
	! If you use this code for your work, please cite the parent paper.
	! Author: F. Hess (dr.franziska.hess@gmail.com)
	! Feel free to contact me if there is any problem and please notify me of any bugs you find.
	! This program was tested for ifort 11.1 and ifort 18.0.1  (no flags required).
	! Sorry, it doesn't work with gfortran.
	!
	!-------------------------------------------------------------------------------------------------------------------------------
	! Instructions
	!-------------------------------------------------------------------------------------------------------------------------------
	! 
	! The timing is called by a line
	! call time_subtractive_scheme (update_function, gamma_function, partial_gamma_function, search_function, energy_function, size_list, time_blind, amma_blind)
	! It takes functions as arguments and runs the timing test calling these functions, but
	! under otherwise identical conditions. 
	!
	!-------------------------------------------------------------------------------------------------------------------------------
	!
	real,allocatable, dimension(:)::random_array
	!-------------------------------------------------------------------------------------------------------------------------------
	! "Cluster Expansion"
	! For details on this implementation of the cluster expansion and the supercluster contraction
	! please refer to the the program cluster_contraction.f90.
	!-------------------------------------------------------------------------------------------------------------------------------
	real,    allocatable, dimension(:,:)		 :: interaction_base_values     ! just some number to generate a "Cluster Expansion"
	real,    allocatable, dimension(:)			 :: e0                          ! adsorption energy at zero coverage
	real,    allocatable, dimension(:,:,:,:)	 :: pairwise_interactions       ! Pairwise interaction parameters
	real,    allocatable, dimension(:,:,:)		 :: Lshape_1,Lshape_2, Lshape_3 ! Three-body interaction parameters -- not actually assigned in the present program, used only for simulating the time for literal evaluation
	integer, allocatable, dimension(:,:)		 :: interacting_sites           ! A tweak to speed up literal evaluation on some architectures, cf energy_literal_long vs. energy_literal_short
	integer, allocatable, dimension(:,:)		 :: i_rank_2_temporary          ! dummy
	integer, allocatable, dimension(:,:,:)       :: startconfs
	
	integer, allocatable, dimension(:,:) :: site_occupation
	integer, allocatable, dimension(:)   :: wrap_site
	real,    allocatable, dimension(:,:,:) :: site_rate
	real,    allocatable, dimension(:):: time_blind, time_total_gamma, time_total_random, time_local, time_global, time_partial, time_SCD, time_zero
	real,    allocatable, dimension(:):: gamma_blind, gamma_total_gamma, gamma_total_random, gamma_local, gamma_global, gamma_partial, gamma_SCD, gamma_zero
	
	integer, parameter:: R = 2
	integer:: i_rand =  1
	integer:: nx, ny, n_sizes
	integer:: n_int
	real:: old_gamma, new_gamma
	!-------------------------------------------------------------------------------------------------------------------------------
	! Physical constants
	
	real:: kB = 1.38E-23
	real:: h  = 6.626E-34
	real:: T  = 300
	real:: eV = 1.602E-19
	
	!-------------------------------------------------------------------------------------------------------------------------------
	! Other things you can play with
	
	integer:: nr = 5           ! rates per site
	integer:: n_species  = 4   ! number of species
	integer:: verbosity  = 0   ! how much output to write ( 0 = only main results, 1 = additional output, 2 = all output)
	integer:: n_steps    = 100000
	
	! Fibonacci-like numbers. This algorithm works with any size, though
	integer, dimension(12) :: size_list =(/5,8,13,21,34,55,89,144,233,377,610,987/)
	
	! SCD (supercluster-related)
	real,    allocatable, dimension(:,:,:,:) :: SCD_I
	real,    allocatable, dimension(:,:,:)       :: SCD_II
	integer, dimension(1:2,1:3)::SCD_I_sites  = reshape((/(/-1,0/),(/-1,1/),(/0,1/)/),shape(SCD_I_sites))
	integer, dimension(1:2,1:2)::SCD_II_sites = reshape((/(/ 2, 0/),(/ 0, 2/)/),shape(SCD_II_sites))
	
	!-------------------------------------------------------------------------------------------------------------------------------
	! defining some lateral interactions 
	! base values used for generating an interaction parameter set (this is not a real cluster expansion, just some numbers)
	! define only upper triangular matrix due to symmetry (interaction of 1 with 2 is equal to interaction of 2 with 1)
	! !!! IF YOU CHANGE n_species YOU NEED TO CHANGE THIS AS WELL !!!
	! Also please refer to cluster_contraction.f90 for more details.
	!-------------------------------------------------------------------------------------------------------------------------------
	integer, dimension(-R:R)::rseq
	allocate (e0(0:n_species))
	e0 = (/0.0, -1.0, -1.1, -1.2, -1.3/)
	allocate (interaction_base_values	(0:n_species, 0:n_species))
	interaction_base_values=0
	interaction_base_values(1,1:n_species)=(/0.10, 0.05, 0.08, 0.09/)
	interaction_base_values(2,2:n_species)=(/0.02, 0.03, 0.04/)
	interaction_base_values(3,3:n_species)=(/0.11, 0.12/)
	interaction_base_values(4,4:n_species)=(/0.07/)
	
	
	
	
	call random_seed() ! initialize RNG
	call initialize_data_structures
	call generate_interaction_parameter_set(n_species, R, interaction_base_values, pairwise_interactions)
	call identify_interacting_sites(R, n_int)
	call precalculate_SCD_clusters(time_SCD(1))
	call generate_starting_configurations
	
	! This part calls the timing for different combinations of algorithms to determine the combined speedup
	
	write(*,*) 'No calculation'
	call time_subtractive_scheme (no_update,     no_gamma,              no_partial_gamma, no_search, energy_literal_short,  size_list, time_blind,   gamma_blind)
	
	write(*,*)
	write(*,*) 'Standard algorithm'
	write(*,*) 'Total Gamma, random search, local update'
	call time_subtractive_scheme (local_update_small,  total_gamma_summation, no_partial_gamma, random_site, energy_literal_short, size_list, time_local,   gamma_local)
	
	write(*,*)
	write(*,*) 'Subtraction Scheme'
	write(*,*) 'Partial Gamma, random search, local update, literal Cluster Expansion'
	call time_subtractive_scheme (local_update_small,  subtractive_gamma,     partial_gamma,    random_site, energy_literal_short, size_list, time_partial, gamma_partial)
	
	write(*,*) 'Partial Gamma, random search, local update, SC D '
	call time_subtractive_scheme (local_update_small,  subtractive_gamma,     partial_gamma,    random_site, energy_SCD, size_list, time_SCD, gamma_SCD)
	
	write(*,*) 'Partial Gamma, random search, local update, zero_energy'
	call time_subtractive_scheme (local_update_small,  subtractive_gamma,     partial_gamma,    random_site, no_energy, size_list, time_SCD, gamma_SCD)
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	stop
	contains
	!-----------------------------------------------------------------------------------------------------------------------------------
	subroutine time_subtractive_scheme(update_algorithm, accounting_algorithm, partial_gamma_algorithm, search_algorithm, energy_function, size_list, elapsed_time, final_gamma)
		implicit none
		!-----------------------------------------------------------------------------------------------------------------------------------
		! This subroutine times different gamma update and energy computation schemes. If you want to test
		! search as well, please use the program supersite.f90
		! It essentially runs and times a short (n_steps) KMC simulation (without proper search) for different 
		! lattice sizes.
		!-----------------------------------------------------------------------------------------------------------------------------------
		integer :: is, step
		integer, dimension(1:n_sizes), intent(in) :: size_list
		real, dimension(1:n_sizes), intent(out) :: elapsed_time
		real, dimension(1:n_sizes), intent(out) :: final_gamma
		real:: accounting_algorithm, partial_gamma_algorithm
		real:: energy_function
		integer:: search_algorithm, update_algorithm
		integer:: dummy, x_sel, y_sel
		integer, dimension(1:2)::site_dummy = (/1,1/), upsite
		real, dimension(1:n_sizes):: relative_speedup, gamma_error
		real:: time_start, time_end, old_partial_gamma, new_partial_gamma, real_dummy=0
		
		! Interfaces for function arguments
		!interface
		!	integer function update_algorithm(energy_function, x0,y0, R, interacting_sites) 
		!		real:: energy_function
		!		integer:: x0,y0,R
		!		integer, dimension(:,:):: interacting_sites
		!		end function update_algorithm
		!end interface
		
		do is = 1, n_sizes
			write(*,*)size_list(is)
			call cpu_time(time_start)
			i_rand = 1
			nx = size_list(is)
			ny = nx
			allocate (site_occupation(1:nx, 1:ny))
			allocate (wrap_site(1-R:nx+R))
			allocate (site_rate(1:nx, 1:ny, 1:nr))
			call assign_wrap_sites
			site_occupation = startconfs(is,1:nx,1:ny)
			
			! first rate calculation must be global
			dummy = global_update(energy_function,site_dummy,interacting_sites)
			old_partial_gamma = 0
			new_partial_gamma = 0
			new_gamma = total_gamma_summation(0.,0.,0., gamma_error(is))
			
			do step = 1, n_steps
				
				!Update GAMMA
				old_gamma = new_gamma 
				new_gamma = accounting_algorithm(new_gamma, old_partial_gamma, new_partial_gamma, gamma_error(is))
				
				! Search
				dummy = search_algorithm(new_gamma,random_array(i_rand),upsite)
				i_rand = i_rand + 1 
				
				! Calculate gamma_old
				old_partial_gamma = partial_gamma_algorithm(upsite, interacting_sites)
				
				! Execute reaction
				call change_site(site_occupation, upsite)
				
				! Update
				dummy = update_algorithm(energy_function, upsite, interacting_sites)
				
				! Calculate gamma_new
				new_partial_gamma = partial_gamma_algorithm(upsite, interacting_sites)
				i_rand = i_rand+1
			end do
			
			deallocate (site_occupation)
			deallocate (site_rate)
			deallocate (wrap_site)
			
			call cpu_time(time_end)
			elapsed_time(is) = time_end - time_start
			final_gamma(is)  = new_gamma
		end do
		write(*,*)
		relative_speedup=(elapsed_time-time_blind) / (elapsed_time(1) - time_blind(1))
		write(*,'(i5, 1x, f14.6,1x,f14.6,1x , e14.7,1x,e14.7,1x)') (size_list(is), elapsed_time(is), relative_speedup(is), final_gamma(is), gamma_error(is), is = 1, n_sizes)
		write(*,*)
	end subroutine time_subtractive_scheme
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	! GAMMA summation schemes
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function total_gamma_summation(dummy1, dummy2, dummy3, gamma_error)
		implicit none
		! Calculates GAMMA from scratch in every step
		real, intent(in) :: dummy1, dummy2, dummy3
		real, intent(inout):: gamma_error
		integer:: x, y,ri
		
		! Just add all the rates together
		total_gamma_summation =  sum(site_rate(:,:,:))
		
		! Error for single precision
		gamma_error = total_gamma_summation*1E-7
		
	end function total_gamma_summation
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function subtractive_gamma(old_gamma, old_partial_gamma, new_partial_gamma, gamma_error)
		implicit none
		! This is the subtractive GAMMA scheme.
		real, intent(in) :: old_gamma, old_partial_gamma, new_partial_gamma
		real, intent(inout):: gamma_error
		integer:: x, y,ri
		
		! Subtract the old partial gamma, then add the new partial gamma 
		subtractive_gamma = old_gamma - old_partial_gamma + new_partial_gamma
		
		! ACCUMULATE Error (single precision)
		gamma_error = gamma_error + maxval((/subtractive_gamma, old_partial_gamma, new_partial_gamma/))*1E-7
		
	end function subtractive_gamma
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	real function no_gamma(dummy1, dummy2, dummy3, gamma_error)
		implicit none
		! This function does nothing and returns 0.
		real, intent(in)::dummy1, dummy2, dummy3
		real, intent(inout)::gamma_error
		
		gamma_error = 0
		no_gamma = 0
		
	end function no_gamma
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	! Partial gamma summation schemes
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	pure real function partial_gamma(site, interacting_sites)
		implicit none
		
		! Calculates gamma_old and gamma_new
		
		integer, dimension(1:n_int,1:2), intent(in)::interacting_sites
		integer, dimension(1:2), intent(in) :: site
		integer::x,y, i_int
		integer::x0,y0
		partial_gamma = 0
		x0 = site(1)
		y0 = site(2)
		
		! Sum reaction rates within the range of lateral interactions
		do i_int = 1, n_int
			x = wrap_site (x0+interacting_sites(i_int,1))
			y = wrap_site (y0+interacting_sites(i_int,2))
			partial_gamma = partial_gamma + sum(site_rate(x,y,:))
		end do
	end	function partial_gamma
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	pure real function no_partial_gamma(site, interacting_sites)
		implicit none
		! This function does nothing and returns 0.
		integer, dimension(1:n_int,1:2), intent(in)::interacting_sites
		integer, dimension(1:2), intent(in) :: site
		integer:: x0,y0
		no_partial_gamma = 0
	end function no_partial_gamma
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	! update_schemes
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	pure integer function no_update(real_dummy,site,int_dummy)
		implicit none
		! This function does nothing and returns 0.
		integer, intent(in) :: int_dummy
		real,    intent(in) :: real_dummy
		integer, dimension(1:2), intent(in):: site
		
		no_update = 0
		
	end function no_update
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	integer function local_update_small(energy_function, site, interacting_sites)
		implicit none
		! Updates reaction rates around reaction event
		integer, dimension(1:2), intent(in):: site ! site of reaction event
		integer, dimension(1:n_int,1:2), intent(in)::interacting_sites
		integer, dimension(1:2):: csite ! site to be updated
		real:: energy_function
		real:: energy
		integer::x0,y0
		integer::x,y, wy, i_int
		integer, dimension(-R:R,-R:R):: local_surrounding_sites
		
		x0 = site(1)
		y0 = site(2)
		
		! loop over interacting sites
		do i_int = 1, n_int
			x = wrap_site(x0+interacting_sites(i_int,1))
			y = wrap_site(y0+interacting_sites(i_int,2))
			csite=(/x,y/)
			
			! extract piece of configuration
			call get_surrounding_sites(site_occupation,csite, local_surrounding_sites(-R:R,-R:R))
			
			! calculate energy
			energy = energy_function(local_surrounding_sites)
			
			! reset rates, just because.
			site_rate(x,y,1:5)=0
			! calculate new rate (only desorption)
			site_rate(x,y,1)=desorption_rate(energy)
		end do
		local_update_small = 0
		
	end function local_update_small
	

	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine get_surrounding_sites(whole_conf,site, surrounding_sites)
		implicit none
		! Extracts a piece of the configuration
		integer, dimension(1:nx,1:ny), intent(in)::whole_conf
		integer, dimension(-R:R, -R:R), intent(out)::surrounding_sites
		integer, dimension(1:2),intent(in)::site
		integer::x0,y0
		integer::x,y
		x0 = site(1)
		y0 = site(2)
		do x = -R, R
			do y = -R, R
				surrounding_sites(x,y) = whole_conf(wrap_site(x + x0),wrap_site(y + y0))
			end do
		end do
	end subroutine get_surrounding_sites
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	integer function global_update(energy_function, site, interacting_sites)
		implicit none
		! "We don't do global update here." Please use local_update.
		integer, dimension(1:2), intent(in)::site
		integer, dimension(1:n_int,1:2), intent(in)::interacting_sites
		integer::x0,y0
		integer::x,y
		real:: energy
		real:: energy_function
		integer, dimension(-R:R,-R:R):: surrounding_sites
		
		do x = 1, nx
			do y = 1, ny
				call get_surrounding_sites(site_occupation, (/x,y/), surrounding_sites)
				energy = energy_function(surrounding_sites)
				site_rate(x,y,:)=0
				site_rate(x,y,1)=desorption_rate(energy)
			end do
		end do
		global_update = 0
	end function global_update
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	! search_schemes
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	integer function no_search(gamma,r_rand, site)
		implicit none
		! This function does nothing and returns 0.
		! selected site is always (/1,1/)
		integer, dimension(1:2),intent(out)::site
		real, intent(in)::gamma,r_rand
		
		site = (/1,1/)
		no_search = 0
	end function no_search
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	integer function random_site(gamma, r_rand, site) 
		implicit none
		! Select random site
		integer, dimension(1:2),intent(out)::site
		real, intent(in)::gamma,r_rand
		integer:: r_int, nsites
		integer:: x,y
		
		nsites = nx*ny
		r_int  = int(nsites*r_rand+1)
		
		y = mod(r_int-1,nx)+1
		x = (r_int-y)/nx+1
		random_site = 0
		site=(/x,y/)
	end function random_site
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine initialize_data_structures
		implicit none
		integer::i
		n_sizes = size(size_list)
		allocate (pairwise_interactions		(0:n_species, 0:n_species, 0:R, 0:R)) ! define only one quarter of interaction matrix for c2v substrate symmetry
		allocate (Lshape_1(0:n_species,0:n_species,0:n_species))
		allocate (Lshape_2(0:n_species,0:n_species,0:n_species))
		allocate (Lshape_3(0:n_species,0:n_species,0:n_species))
		allocate (interacting_sites(1:4*R*R, 1:2)) ! temporary allocation, will be changed later
		
		allocate(time_blind(1:n_sizes))
		time_blind = 0
		allocate(gamma_blind(1:n_sizes))
		allocate(time_total_gamma(1:n_sizes))
		allocate(gamma_total_gamma(1:n_sizes))
		allocate(time_total_random(1:n_sizes))
		allocate(gamma_total_random(1:n_sizes))
		allocate(time_local(1:n_sizes))
		allocate(gamma_local(1:n_sizes))
		allocate(time_global(1:n_sizes))
		allocate(gamma_global(1:n_sizes))
		allocate(time_partial(1:n_sizes))
		allocate(gamma_partial(1:n_sizes))
		allocate(time_SCD(1:n_sizes))
		allocate(gamma_SCD(1:n_sizes))
		allocate(time_zero(1:n_sizes))
		allocate(gamma_zero(1:n_sizes))
		allocate(random_array(1:n_steps*3))
		allocate(startconfs(1:n_sizes,1:maxval(size_list),1:maxval(size_list)))
		call random_number(random_array)
		if(verbosity >=2)write(*,*) random_array
		rseq=(/(i , i = -R, +R)/)
	end subroutine initialize_data_structures
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine generate_starting_configurations
		implicit none
		integer::is,x,y
		real:: r_rand
		startconfs = 0
		do is = 1, n_sizes
			nx = size_list(is)
			ny = size_list(is)
			do x = 1, nx
				do y = 1, ny
					call random_number(r_rand)
					startconfs(is, x, y) = int(r_rand*5)
				end do
			end do
			if(verbosity >=1)then
				write(*,*)'starting configuration',is, ', size ', nx
				do y = 1, ny
					write(*,'(1000(1x,i1))')(startconfs(is,x,y),x = 1,nx)
				end do
			end if
		end do
		
	end subroutine generate_starting_configurations
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	! Things related to energy computation. Please refer to supercluster.f90 for details on these.
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine generate_interaction_parameter_set(n_species, R, interaction_base_values, pairwise_interactions)
		implicit none
		integer,intent(in)												:: n_species, R
		real, dimension(0:n_species,0:n_species), intent(inout)			:: interaction_base_values
		real, dimension(0:n_species,0:n_species,0:R,0:R), intent(inout)	:: pairwise_interactions
		integer :: s1, s2
		integer :: x, y
				
		! initialize
		pairwise_interactions = 0
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
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine identify_interacting_sites(R,n_int)
		implicit none
		integer, intent(in)::R
		integer,intent(out):: n_int
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
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	pure real function desorption_rate(energy)
		implicit none
		real, intent(in)::energy
		desorption_rate = kB*T/h* exp( + energy*eV/(kB*T)) ! + energy for "desorption"
	end function desorption_rate
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	pure real function energy_literal_short(configuration)
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
		
		e = e + Lshape_1(sc,configuration(-1,0),configuration(-1, 1)) !b
		e = e + Lshape_1(sc,configuration( 1,0),configuration( 1, 1)) !b
		e = e + Lshape_1(sc,configuration(-1,0),configuration(-1,-1)) !b
		e = e + Lshape_1(sc,configuration( 1,0),configuration(-1,-1)) !b
		e = e + Lshape_2(sc,configuration(0,-1),configuration(-1, 1)) !c
		e = e + Lshape_2(sc,configuration(0, 1),configuration( 1, 1)) !c
		e = e + Lshape_2(sc,configuration(0,-1),configuration(-1,-1)) !c
		e = e + Lshape_2(sc,configuration(0, 1),configuration(-1,-1)) !c
		e = e + Lshape_3(sc,configuration(0, 1),configuration( 1, 0))
		e = e + Lshape_3(sc,configuration(0,-1),configuration( 1, 0))
		e = e + Lshape_3(sc,configuration(0, 1),configuration(-1, 0))
		e = e + Lshape_3(sc,configuration(0,-1),configuration(-1, 0))
		
		energy_literal_short=e
	end function energy_literal_short
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine change_site(conf, site)
		
		! Simple adsorption/desorption mechanism
		! If a site is occupied, the molecule is removed
		! Otherwise a random molecule is put on the site
		! Again -- this is not a real KMC simulation
		
		implicit none
		integer, dimension(1:nx,1:ny), intent(inout) :: conf  ! the current configuration
		integer, dimension(1:2), intent(in)::site
		integer::x,y
		x = site(1)
		y = site(2)
		if(conf(x,y)==0)then ! site free
			conf(x,y)=int(random_array(i_rand)*4)+1  ! put random molecule on site
		else
			conf(x,y)=0  ! remove molecule
		end if
		i_rand=i_rand+1 ! 1 random number used
		
	end subroutine change_site
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine precalculate_SCD_clusters(elapsed_time)
		implicit none
		real:: time_start, time_end
		real, intent(out):: elapsed_time
		integer,dimension(-R:R,-R:R):: SCD_conf 
		integer:: s1, s2, s3, s4, s5, s6
		
		call cpu_time(time_start)
		!allocate(res_energy_SCD (1:n_rnd_conf))
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
	
	pure real function energy_only_pairwise(configuration)
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
	
	pure real function energy_SCD(configuration)
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
	
	pure real function no_energy(configuration)
		implicit none
		integer, dimension(-R:R,-R:R), intent(in):: configuration
		no_energy = 0
	end function no_energy
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	recursive integer function wrap_around (x,nx)
		! Implementation of periodic boundary conditions.
		! This wraps around an arbitrary number of times in the recursive implementation.
		! Different wrap-around algorithms were tested, and all were found too slow for KMC purposes as the function 
		! has to be called every time a site is referenced in the code. Therefore, instead of relying on an actual
		! function, we decided to precalculate a wrap-around table (array wrap_sites), which is assigned at the start
		! of the program (see below subroutine assign_wrap_sites).
		implicit none
		integer,intent(in):: x  ! coordinate
		integer,intent(in):: nx ! range (1:nx)
		if(x<1)then         ! exceeding lower bound
			wrap_around = wrap_around(x+nx, nx) ! add range and wrap around again
		else if(x>nx)then   ! exceeding upper bound
			wrap_around = wrap_around(x-nx, nx)	! subtract range and wrap around again
		else 
			wrap_around = x ! not exceeding bounds
		end if
		
	end function wrap_around
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
	subroutine assign_wrap_sites
		! Assigns the wrap-around table wrap_sites using the above wrap_around function.
		! Precalculating the wrap-around table significantly speeds up evaluation of other parts of the code.
		implicit none
		integer :: x ! coordinate
		
		do x = 1 - R, nx + R ! precalculate for the range of lateral interactions
			wrap_site(x) = wrap_around(x,nx) ! use the above wrap_around function
		end do
	end subroutine assign_wrap_sites
	
	!-----------------------------------------------------------------------------------------------------------------------------------
	
end program subtraction_scheme
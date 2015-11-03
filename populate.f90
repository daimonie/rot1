!f90
!f2py --fcompiler=gfortran --f90flags="-fopenmp" -lgomp -m -c populate populate.f90
!Earlier, I wrote a `dirty' function in python that is significantly slowing things down. Because I see no other implementation
! that works, and stackexchange hasn't given me an answer either, I've chosen to write that function in openMP fortran instead.

module populate
	contains
	subroutine fermi_integrand( first, second, third, fourth, fermi_surface, radius, angle, angle_array, fermi)
	!fermi has to be in the argument list because it is going out. Fortran is weird that way.
		implicit none
		
		!flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);
		integer, intent (in):: first, second, third, fourth
		double precision, intent (in), dimension(fourth) :: fermi_surface, angle_array
		double precision, intent (in), dimension(first,second,third,fourth):: radius, angle
		
		double precision, intent(out), dimension(first, second, third, fourth) :: fermi
		
		double precision :: present_surface
		integer :: i, j, n, m, a, omp_get_max_threads 
		!while I would prefer using openMP, it causes a segfault I can't solve.
		call omp_set_num_threads(omp_get_max_threads())
		!$omp parallel do  &
		!$omp default(shared) &
		!$omp private(present_surface) 
		do m = 1, fourth 
			do n = 1, third
				do j = 1, second
					do i = 1, first
						!current angle is phi[i,j,ii,jj]
						!we have to find the current fermi surface
						present_surface = 0.
						do a = 1, fourth    
							if (angle(i,j,n,m) == angle_array(a)) then
								present_surface = fermi_surface(a)
							end if
						end do 
						if (radius(i,j,n,m) >= present_surface) then 
							fermi(i,j,n,m) = 0.
						else 
							fermi(i,j,n,m) = 1.
						end if			
					end do		
				end do		
			end do
		end do  
		!$omp end parallel do
	end subroutine fermi_integrand
	
	subroutine fermi_contour( first, second, fermi_surface, radius, angle, angle_array, fermi)
	!fermi has to be in the argument list because it is going out. Fortran is weird that way.
		implicit none
		
		!flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);
		integer, intent (in):: first, second
		double precision, intent (in), dimension(first) :: fermi_surface, angle_array
		double precision, intent (in), dimension(first,second):: radius, angle
		
		double precision, intent(out), dimension(first, second) :: fermi
		
		double precision :: present_surface
		integer :: i, j, n, m, a, omp_get_max_threads 
		!while I would prefer using openMP, it causes a segfault I can't solve.
		call omp_set_num_threads(omp_get_max_threads())
		!$omp parallel do  &
		!$omp default(shared) &
		!$omp private(present_surface) 
		do j = 1, second
			do i = 1, first
				!current angle is phi[i,j,ii,jj]
				!we have to find the current fermi surface
				present_surface = 0.
				do a = 1, second    
					if (angle(i,j) == angle_array(a)) then
						present_surface = fermi_surface(a)
					end if
				end do 
				if (radius(i,j) >= present_surface) then 
					fermi(i,j) = 0.
				else 
					fermi(i,j) = 1.
				end if			
			end do		
		end do		
		!$omp end parallel do
	end subroutine fermi_contour
end module populate
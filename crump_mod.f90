!
! Computes a potential field in a Cartesian box from the 
! normal component of Bn on the six boundaries. 
!
! BSD 3-Clause License
!
! Copyright (c) 2024, S.A. Gilchrist

! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its
!  contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.!

MODULE CRUMP_MOD

  USE ISO_C_BINDING,ONLY: C_INT64_T,C_DOUBLE,C_INT

  IMPLICIT NONE

  ! Set kind parameters. To use with Python, these must be C_DOUBLE
  ! and C_INT
  INTEGER,PARAMETER :: FP = C_DOUBLE
  INTEGER,PARAMETER :: IT = C_INT64_T
   
  ! MATH CONSTANTS
  REAL(FP),PARAMETER :: DPI = ACOS(REAL(-1,FP))
  
  ! INTEGER(IT) OPTIONS LENGTH AND FIELDS
  !
  INTEGER(C_INT),PARAMETER :: OPT_LENGTH = 8         !< Length of IOPT vector
  INTEGER(C_INT),PARAMETER :: IOPT_NMAX  = 0         !< Max number of iterations
  INTEGER(C_INT),PARAMETER :: IOPT_LITER = 1         !< Last iteration 

  ! REAL OPTION LENGTH A FIELDS
  !
  INTEGER(C_INT),PARAMETER :: ROPT_TOL   = 0   !< Solution tolerance
  INTEGER(C_INT),PARAMETER :: ROPT_DELTA = 1   !< Change in solution at last iteration

  ! No error flag
  INTEGER(C_INT),PARAMETER :: NOERROR = 0
    
  ! Interface 
  INTERFACE zero
    MODULE PROCEDURE zero_3D
  END INTERFACE
  
 CONTAINS

! ---------------------------------------------------------------------
!
!>@name crump_solve
!
!>@brief Solves Laplace equation to get scalar potential and magnetic field
!
!
FUNCTION crump_solve(nshape,iopt,ropt,b) BIND(C) RESULT(ierr)
  
  IMPLICIT NONE

  ! INPUT 
  INTEGER(IT) ,DIMENSION(4),INTENT(IN) :: nshape

  ! INPUT/OUTPUT
  REAL(FP)   ,DIMENSION(nshape(1),nshape(2),nshape(3),nshape(4)),INTENT(INOUT) :: b
  REAL(FP)   ,DIMENSION(0:OPT_LENGTH-1)                         ,INTENT(INOUT) :: ropt
  INTEGER(IT),DIMENSION(0:OPT_LENGTH-1)                         ,INTENT(INOUT) :: iopt 

  ! RETURN
  INTEGER(IT) :: ierr

  ! LOCAL
  REAL(FP)    :: dx,dy,dz
  INTEGER(IT) :: i0,j0,k0,i
  INTEGER(IT) :: return_stat 
  REAL(FP)    :: tstart,tend,tol

  ! LOCAL: ARRAYS
  REAL(FP),DIMENSION(:,:,:),ALLOCATABLE :: bcx,bcy,bcz,U_pad
  REAL(FP),DIMENSION(:)    ,ALLOCATABLE :: xarr,yarr,zarr
  
  ! ================
  ! EXTRACT ROPTS
  ! ================
  !
  ! Real-type solution parameters
  !
  tol = ropt(ROPT_TOL)

  ! ================
  ! CONSTRUCT MESH
  ! ================
  !
  ! Construct a uniformly spaced mesh
  !
  
  ! Extract mesh size
  i0 = nshape(1)
  j0 = nshape(2)
  k0 = nshape(3)

  ! Compute mesh spacing
  dx = REAL(1,FP)/(i0-REAL(1,FP))
  dy = dx
  dz = dx
      
  ! Allocate memory for mesh
  ALLOCATE(xarr(i0),yarr(j0),zarr(k0),STAT=return_stat)
  IF(return_stat .ne. 0) THEN 
    STOP "ERROR(main_code): Memory allocation error (xarr,yarr,zarr)"
  ENDIF
  
  ! Construct mesh
  xarr = [(i,i=0,i0-1)]*dx
  yarr = [(i,i=0,j0-1)]*dy
  zarr = [(i,i=0,k0-1)]*dz
    
  ! ===========================
  ! ALLOCATE LOCAL MEMORY
  ! ===========================
  !
  ! Allocate memory for local arrays
  !
  
  ! Allocate memory for Uc. Variable Uc has a buffer around the edges 
  ! which enforces the boundary conditions. Hence, it is teo points 
  ! larger in each direction 
  !
  ALLOCATE(U_pad(i0+2,j0+2,k0+2),STAT=return_stat)
  IF(return_stat .ne. 0) STOP "ERROR(main_code): Memory allocation error (U_pad)"
  
  ! Allocate memory for boundary conditions on x boundaries 
  !
  ALLOCATE(bcx(j0,k0,2),STAT=return_stat)
  IF(return_stat .ne. 0) STOP "ERROR(main_code): Memory allocation error (bcx)"
  
  ! Allocate memory for boundary conditions on y boundaries 
  !
  ALLOCATE(bcy(i0,k0,2),STAT=return_stat)
  IF(return_stat .ne. 0) STOP "ERROR(main_code): Memory allocation error (bcy)"

  ! Allocate memory for boundary conditions on z boundaries 
  !
  ALLOCATE(bcz(i0,j0,2),STAT=return_stat)
  IF(return_stat .ne. 0) STOP "ERROR(main_code): Memory allocation error (bcz)"

  ! Extract boundary conditions from B 
  bcx(:,:,1) = B(1,:,:,1); bcx(:,:,2) = B(i0,:,:,1)
  bcy(:,:,1) = B(:,1,:,2); bcy(:,:,2) = B(:,j0,:,2)
  bcz(:,:,1) = B(:,:,1,3); bcz(:,:,2) = B(:,:,k0,3)

  ! =============
  ! SOLVE SYSTEM
  ! =============
  !
  ! Call SOR subroutine to find solution. The calculation time 
  ! is computed and printed to screen
  !
    
  ! Solve for the vector potential 
  ierr = sor_laplace(i0,j0,k0,xarr,yarr,zarr,bcx,bcy,bcz,iopt,ropt,U_pad)

  ! Compute new magentic field from grad of potential. Set buffer 
  CALL grad(U_pad,dx,B(:,:,:,1:3))

  ! =============
  ! POST PROCESS
  ! =============
  !
  ! 1: Replace normal component at boundaries with boundary
  !    conditions
  !

  ! Set exact boundary conditions 
  B(1,:,:,1) = bcx(:,:,1); B(i0,:,:,1) = bcx(:,:,2)
  B(:,1,:,2) = bcy(:,:,1); B(:,j0,:,2) = bcy(:,:,2)
  B(:,:,1,3) = bcz(:,:,1); B(:,:,k0,3) = bcz(:,:,2)

  ! Copy out potential
  B(:,:,:,4) = U_pad(2:i0-1,2:j0-1,2:k0-1)
 
END FUNCTION
        
! ----------------------------------------------------------------------
!
!>@name sor_laplace
!
!>@brief Solve Laplace's equation using successive over-relaxation (SOR)
!
!>@details
!! Use Check-board SOR to solve linear system. 
!!
!! Output U, must be padded, i.e. it must be of size 
!! [nx+2,ny+2,nz+2].
!!
!
INTEGER(IT) FUNCTION sor_laplace(nx,ny,nz,x,y,z,bcx,bcy,bcz,iopt,ropt,u) RESULT(ierr)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT),INTENT(IN) :: nx,ny,nz
  REAL(FP)   ,INTENT(IN) :: x(nx),y(ny),z(nz)
  REAL(FP)   ,INTENT(IN) :: bcx(ny,nz,2),bcy(nx,nz,2),bcz(nx,ny,2)
  
  ! INPUT/OUTPUT
  REAL(FP)   ,DIMENSION(0:OPT_LENGTH-1),INTENT(INOUT) :: ropt 
  INTEGER(IT),DIMENSION(0:OPT_LENGTH-1),INTENT(INOUT) :: iopt 
  
  ! OUTPUT
  REAL(FP),DIMENSION(nx+2,ny+2,nz+2),INTENT(OUT) :: u ! Note: padded.
  
  ! LOCAL
  REAL(FP),DIMENSION(:,:,:),ALLOCATABLE :: uc ! Local copy
  REAL(FP) :: w,rho,dx,du
  INTEGER(IT) :: i,j,k,nmax
  LOGICAL :: converged 
  
  ! Compute grid spacing in x, y, and z directions respectively. 
  dx = x(2)-x(1)
 
  ! Allocate array to hold old copy and initialize to zero
  ALLOCATE(uc(nx+2,ny+2,nz+2))
  CALL zero(uc)

  ! Compute rho for Gauss-Sidel relaxation
  rho = 1 - DPI**2/REAL(nx,FP)**2
  w   = 1

  ! Set u to zero 
  CALL zero(u)

  ! Set NMAX. Use HUGE, if IOPT_NMAX is negative
  IF(iopt(IOPT_NMAX) > 0) THEN
    nmax = iopt(IOPT_NMAX)
  ELSE
    nmax = HUGE(nmax)
  ENDIF

  ! Iterate SOR until du_tol is reached or the maximum number of 
  ! iterations is reached
  ! 
  converged = .FALSE.
  SOR_LOOP : DO i=1,nmax

    ! Perform relaxation
    uc = u
    CALL orelax_cb_laplace(nx+2,ny+2,nz+2,bcx,bcy,bcz,dx,w,rho,u)
       
    ! Compute max change in solution
    du = maxabs_diff(uc,u)

    ! Break out of loop if change is below the tolerance
    IF(du < ropt(ROPT_TOL) ) THEN
      converged = .TRUE.
      EXIT SOR_LOOP
    ENDIF

  ENDDO SOR_LOOP

  ! Record final change
  ropt(ROPT_DELTA) = du

  ! Print last iteration
  iopt(IOPT_LITER) = i

  ! Set IERR flag
  IF(converged) THEN
    ierr = NOERROR
  ELSE
    ierr = 1
  ENDIF

END FUNCTION

! ---------------------------------------------------------------------
!
!>@name orelax_cb_laplace
!
!>@brief Perform a single sweep of red-black over-relaxation
!
!>@details
!!
!! Performs successive over relaxation on a red-black (checkerboard)
!! grid with Neumann boundary conditions on all six planar 
!! boundaries. 
!

SUBROUTINE orelax_cb_laplace(nx,ny,nz,bcx,bcy,bcz,dx,w,rho,u)

  IMPLICIT NONE
 
  ! INPUT
  INTEGER(IT),INTENT(IN) :: nx,ny,nz
  REAL(FP)   ,INTENT(IN) :: bcx(ny-2,nz-2,2),bcy(nx-2,nz-2,2),bcz(nx-2,ny-2,2)
  REAL(FP)   ,INTENT(IN) :: dx
  REAL(FP)   ,INTENT(IN) :: rho 
 
  ! OUTPUT
  REAL(FP),DIMENSION(nx,ny,nz),INTENT(OUT) :: u
  
  ! LOCAL
  REAL(FP) :: w,meanu
  INTEGER(IT)  :: i,j,k
  REAL(FP),PARAMETER :: inv6 = REAL(1,FP)/REAL(6,FP)
  REAL(FP),PARAMETER :: inv4 = REAL(1,FP)/REAL(4,FP)
  REAL(FP),PARAMETER :: inv2 = REAL(1,FP)/REAL(2,FP)
  REAL(FP),PARAMETER :: o1   = REAL(1,FP)
 
  ! Integer constants
  INTEGER(IT),PARAMETER :: i1 = 1
  INTEGER(IT),PARAMETER :: i2 = 2

  
  ! Update buffer 
  !
  u(1,2:ny-1,2:nz-1) = u(3,2:ny-1,2:nz-1)-2*bcx(:,:,1)*dx; u(nx,2:ny-1,2:nz-1) = u(nx-2,2:ny-1,2:nz-1)+2*bcx(:,:,2)*dx
  u(2:nx-1,1,2:nz-1) = u(2:nx-1,3,2:nz-1)-2*bcy(:,:,1)*dx; u(2:nx-1,ny,2:nz-1) = u(2:nx-1,ny-2,2:nz-1)+2*bcy(:,:,2)*dx
  u(2:nx-1,2:ny-1,1) = u(2:nx-1,2:ny-1,3)-2*bcz(:,:,1)*dx; u(2:nx-1,2:ny-1,nz) = u(2:nx-1,2:ny-1,nz-2)+2*bcz(:,:,2)*dx

  ! =======================
  ! UPDATE SQUARES: BLACK
  ! =======================
  
  ! Update interior points
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  DO k=2,nz-1 
    DO j=2,ny-1
      DO i=2+MOD(j+MOD(k,i2),i2),nx-1,2  
	       

        ! Update U at point (i,j,k)
        u(i,j,k) =u(i,j,k)*(1-w)+w*(u(i+1,j,k)+u(i-1,j,k)+u(i,j-1,k)+u(i,j+1,k)+&
                                      u(i,j,k-1)+u(i,j,k+1))*inv6
                              
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO 
  
  ! Update buffer 
  !
  u(1,2:ny-1,2:nz-1) = u(3,2:ny-1,2:nz-1)-2*bcx(:,:,1)*dx; u(nx,2:ny-1,2:nz-1) = u(nx-2,2:ny-1,2:nz-1)+2*bcx(:,:,2)*dx
  u(2:nx-1,1,2:nz-1) = u(2:nx-1,3,2:nz-1)-2*bcy(:,:,1)*dx; u(2:nx-1,ny,2:nz-1) = u(2:nx-1,ny-2,2:nz-1)+2*bcy(:,:,2)*dx
  u(2:nx-1,2:ny-1,1) = u(2:nx-1,2:ny-1,3)-2*bcz(:,:,1)*dx; u(2:nx-1,2:ny-1,nz) = u(2:nx-1,2:ny-1,nz-2)+2*bcz(:,:,2)*dx

  ! Update over-relaxation parameter  
  !
  IF(w .eq. 1) THEN 
    w = o1/(1 -inv2*rho) 
  ELSE
    w = o1/(1 -w*rho*inv4)
  ENDIF

  ! =======================
  ! UPDATE SQUARES: RED
  ! =======================

  ! Update interior points
  !$omp parallel do private(i,j,k)
  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2+MOD(j+i1+MOD(k,i2),i2),nx-1,2 
	
        ! Update U at point (i,j,k)
        u(i,j,k) =u(i,j,k)*(1 -w)+w*(u(i+1,j,k)+u(i-1,j,k)+u(i,j-1,k)+u(i,j+1,k)+&
                                       u(i,j,k-1)+u(i,j,k+1))*inv6

                             
      ENDDO
    ENDDO
  ENDDO
  !$omp end parallel do 
   
  ! Update over-relaxation parameter 
  !
  w = o1/(1 -w*rho*inv4)

  ! Remove mean. Solution is unique up to the addition of a 
  ! constant, i.e.  u->u+c. Enforcing mean(u) = 0 at each iteration
  ! picks the unique solution with c = 0.
  !
  meanu = 0
  
  ! Compute sum in parallel 
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:meanu) ORDERED
  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-1      
        meanu = meanu + u(i,j,k)      
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
  ! Divide by total
  meanu = meanu/REAL((nx-2)*(ny-2)*(nz-2),FP)

  ! Subtract off
  !$OMP PARALLEL DO PRIVATE(i,j,k) 
  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-1      
        u(i,j,k) = u(i,j,k) - meanu      
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name grad
!
!>@brief Compute the gradient of the magnetic scalar potential
!
SUBROUTINE grad(u,dx,B)

  IMPLICIT NONE 

  ! INPUT
  REAL(FP),DIMENSION(:,:,:),INTENT(IN) :: u
  REAL(FP)                 ,INTENT(IN) :: dx
  
  ! OUTPUT
  REAL(FP),DIMENSION(:,:,:,:),INTENT(OUT) :: B
  
  ! LOCAL
  REAL(FP) :: fac
  INTEGER(IT)  :: i,j,k

  ! Compute spacing factor 
  fac = REAL(1,FP)/(REAL(2,FP)*dx)
  
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  DO k=1,SIZE(B,3)
    DO j=1,SIZE(B,2)
      DO i=1,SIZE(B,1)
        B(i,j,k,1) = fac*(u(i+2,j+1,k+1)-u(i  ,j+1,k+1))
        B(i,j,k,2) = fac*(u(i+1,j+2,k+1)-u(i+1,j  ,k+1))
        B(i,j,k,3) = fac*(u(i+1,j+1,k+2)-u(i+1,j+1,k  ))
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO 

END SUBROUTINE

! ----------------------------------------------------------------------
!
! OPEN MP ASSIGNMENT SUBROUTINES
!
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!
!>@name copy_3D
!
!>@brief Copy 3D real array using OpenMP loop.
!
!>@details
!! Copies B to A. Same as A = B assignment.
!! Both arrays must be the same size.
!
!>@param[in]  B Array that is copied
!>@param[out] A Array that is replaced with copy of B
!
SUBROUTINE copy_3D(A,B)

  IMPLICIT NONE
 
  ! INPUT
  REAL(FP),DIMENSION(:,:,:),INTENT(IN) :: B
  
  ! OUTPUT
  REAL(FP),DIMENSION(:,:,:),INTENT(OUT) :: A
 
  ! LOCAL
  INTEGER(IT) :: i,j,k

  !$OMP PARALLEL DO PRIVATE(i,j,k)
  DO k=1,SIZE(A,3)
    DO j=1,SIZE(A,2)
      DO i=1,SIZE(A,1)
        A(i,j,k) = B(i,j,k)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name maxabs_diff
!>@brief Compute max. abs. difference between A and B
!>@details
!! Equivalent of MAXVAL(ABS(A-B)), but parallel 
!!
!!

FUNCTION maxabs_diff(A,B)

  IMPLICIT NONE

  ! INPUT
  REAL(FP),DIMENSION(:,:,:),INTENT(IN) :: A,B

  ! RETURN
  REAL(FP) :: maxabs_diff

  ! LOCAL
  INTEGER :: i,j,k
  REAL(FP) :: du_abs

  ! Initialize
  maxabs_diff = 0

  !$OMP PARALLEL DO PRIVATE(i,j,k)  REDUCTION(max:maxabs_diff) ORDERED
  DO k=1,SIZE(A,3)
    DO j=1,SIZE(A,2)
      DO i=1,SIZE(A,1)
        du_abs      = ABS(A(i,j,k)-B(i,j,k))
        maxabs_diff = MAX(maxabs_diff,du_abs)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

END FUNCTION

! ----------------------------------------------------------------------
!
!>@name zero_3D
!
!>@brief Set a real array to zero using OpenMP loop.
!
!>@details
!! Should be used to initialize arrays in order for good NUMA.
!!
!!
!>@param[inout] A Array that is zeroed
!
SUBROUTINE zero_3D(A)

  IMPLICIT NONE
 
  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(:,:,:),INTENT(INOUT) :: A
 
  ! LOCAL
  INTEGER(IT) :: i,j,k

  !$OMP PARALLEL DO PRIVATE(i,j,k)
  DO k=1,SIZE(A,3)
    DO j=1,SIZE(A,2)
      DO i=1,SIZE(A,1)
        A(i,j,k) = 0
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
END SUBROUTINE

! ----------------------------------------------------------------------

!>@name get_opt_length
INTEGER(C_INT) FUNCTION get_opt_length() BIND(C) RESULT(res)
  IMPLICIT NONE 
  res = OPT_LENGTH
END FUNCTION

!>@name get_ropt_tol
INTEGER(C_INT) FUNCTION get_ropt_tol() BIND(C) RESULT(res)
  IMPLICIT NONE 
  res = ROPT_TOL
END FUNCTION

!>@name get_ropt_delta
INTEGER(C_INT) FUNCTION get_ropt_delta() BIND(C) RESULT(res)
  IMPLICIT NONE 
  res = ROPT_DELTA
END FUNCTION

!>@name get_iopt_nmax
INTEGER(C_INT) FUNCTION get_iopt_nmax() BIND(C) RESULT(res)
  IMPLICIT NONE 
  res = IOPT_NMAX
END FUNCTION

!>@name get_iopt_liter
INTEGER(C_INT) FUNCTION get_iopt_liter() BIND(C) RESULT(res)
  IMPLICIT NONE 
  res = IOPT_LITER
END FUNCTION

!>@name get_noerror
INTEGER(C_INT) FUNCTION get_noerror() BIND(C) RESULT(res)
  IMPLICIT NONE 
  res = NOERROR
END FUNCTION

END MODULE

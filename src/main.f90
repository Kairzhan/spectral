! This program is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
! more details.
!
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <https://www.gnu.org/licenses/>. 
program main
  use, intrinsic :: iso_c_binding 
  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_io
  use mpi_f08
  implicit none

  ! Mesh resolution, cubic mesh is assumed
  integer, parameter :: n=64, nx=n, ny=n, nz=n
  integer, parameter :: nxnynz=nx*ny*nz
  integer, parameter :: nhp1=n/2+1
  integer, parameter :: ncube=n**3
  
  ! Processor grid configuration
  integer, parameter :: p_row=2, p_col=2
  
  integer :: ierr
  character(len=16) :: filename
  
  ! Arrays in physical space
  real(8), allocatable, dimension(:, :, :) :: ux, uy, uz
  real(8), allocatable, dimension(:, :, :) :: vortx, vorty, vortz
  real(8), allocatable, dimension(:, :, :) :: velvortx, velvorty, velvortz
  real(8), allocatable, dimension(:, :, :) :: uxg, uyg, uzg

  ! Arrays in spectral space
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat, uy_hat, uz_hat
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: ux_hat2, uy_hat2, uz_hat2
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: vortx_hat, vorty_hat, vortz_hat
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: velvortx_hat, velvorty_hat, velvortz_hat
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :, :) :: rhsx, rhsy, rhsz, rhsx0, rhsy0, rhsz0
  real(8), allocatable, dimension(:, :, :) :: kx, ky, kz, k2, k2inv, kabs, ksqr
  integer, allocatable, dimension(:, :, :) :: ik2
  
  integer, dimension(3) :: fft_start, fft_end, fft_size

  real(8), parameter :: pi=3.1415927410125732, pi2=2*pi

  ! Kinematic viscosity
  real(8), parameter :: visc = 1/300.
  
  real(8) :: dx, dy, dz
  integer :: istep

  ! Number of timesteps
  integer, parameter :: nsteps=2000
  ! Dump velocity each ... steps
  integer, parameter :: ndump=50
  ! Output to ASCII files every ... steps
  integer, parameter :: noutput=50
  ! \Delta t
  real(8), parameter :: dt=0.005

  real(8) :: klimit=n/3.
  integer :: i, j, k
  real :: start, finish
  
  call MPI_Init(ierr)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init
  
  dx=pi2/nx
  dy=pi2/ny
  dz=pi2/nz

  ! Allocate arrays
  allocate(uxg(nx, ny, nz))
  allocate(uyg(nx, ny, nz))
  allocate(uzg(nx, ny, nz))

  call alloc_x(ux, opt_global=.true.)
  call alloc_x(uy, opt_global=.true.)
  call alloc_x(uz, opt_global=.true.)
  call alloc_x(vortx, opt_global=.true.)
  call alloc_x(vorty, opt_global=.true.)
  call alloc_x(vortz, opt_global=.true.)
  call alloc_x(velvortx, opt_global=.true.)
  call alloc_x(velvorty, opt_global=.true.)
  call alloc_x(velvortz, opt_global=.true.)
  
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

  if (nrank==1) then
     print *, xstart(1), xend(1)
     print *, xstart(2), xend(2)
     print *, xstart(3), xend(3)
     print *, fft_start(1),fft_end(1)
     print *, fft_start(2),fft_end(2)
     print *, fft_start(3),fft_end(3)
  endif
 
  ! ux_hat
  allocate (ux_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (uy_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (uz_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  ! ux_hat2
  allocate (ux_hat2( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (uy_hat2( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (uz_hat2( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  ! vortx_hat
  allocate (vortx_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (vorty_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (vortz_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  ! velvortx_hat
  allocate (velvortx_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (velvorty_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (velvortz_hat( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  ! rhs
  allocate (rhsx( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (rhsy( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (rhsz( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  allocate (rhsx0( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (rhsy0( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (rhsz0( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  ! kx
  allocate (kx( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (ky( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (kz( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  allocate (k2( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (k2inv( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (kabs( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  allocate (ik2( &
       fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))

  ! wavenumbers
  do i=fft_start(1), fft_end(1)
     if (i<=nx/2) then
        kx(i, :, :)=i-1
     else
        kx(i, :, :)=-n+i-1
     endif
  enddo

  do j=fft_start(2), fft_end(2)
     if (j<=ny/2) then
        ky(:, j, :)=j-1
     else
        ky(:, j, :)=-n+j-1
     endif
  enddo

  do k=fft_start(3), fft_end(3)
     if (k<=nz/2) then
        kz(:, :, k)=k-1
     else
        kz(:, :, k)=-n+k-1
     endif
  enddo

  do k=fft_start(3), fft_end(3)
     do j=fft_start(2), fft_end(2)
        do i=fft_start(1), fft_end(1)
           k2(i,j,k)=kx(i,j,k)**2+ky(i,j,k)**2+kz(i,j,k)**2
           kabs(i,j,k)=sqrt(k2(i,j,k))
           if (k2(i,j,k)/=0) then
              k2inv(i,j,k)=1./k2(i,j,k)
           else
              k2inv(i,j,k)=1.
           endif
        enddo
     enddo
  enddo
  
  ik2 = int ( sqrt(k2) + 0.5 )
    
  ! Initialize velocity field
  do k=xstart(3), xend(3)
     do j=xstart(2), xend(2)
        do i=xstart(1), xend(1)
           ux(i,j,k) = sin((i-0.5)*dx)*cos((j-0.5)*dy)*cos((k-0.5)*dz)
           uy(i,j,k) = -cos((i-0.5)*dx)*sin((j-0.5)*dy)*cos((k-0.5)*dz)
           uz(i,j,k) = 0
        enddo
     enddo
  enddo

  call decomp_2d_fft_3d(ux, ux_hat)
  call decomp_2d_fft_3d(uy, uy_hat)
  call decomp_2d_fft_3d(uz, uz_hat)

  ux_hat=ux_hat/nxnynz
  uy_hat=uy_hat/nxnynz
  uz_hat=uz_hat/nxnynz

  !------------------------------------------------------------
  ! Main time loop
  !------------------------------------------------------------
  
  call cpu_time(start)
  start=MPI_Wtime()
  
  do istep=1, nsteps
     if (mod(istep, 10)==0 .and. nrank==0) then
        write(*,"(I7,1X,1pe14.6)") istep, istep*dt
     endif
     
     ! Dealiasing
     where (ik2>klimit)
        ux_hat=0
        uy_hat=0
        uz_hat=0
     endwhere
     
     ! Get the latest velocity
     ux_hat2=ux_hat
     uy_hat2=uy_hat
     uz_hat2=uz_hat
     call decomp_2d_fft_3d(ux_hat2, ux)
     call decomp_2d_fft_3d(uy_hat2, uy)
     call decomp_2d_fft_3d(uz_hat2, uz)

     ! Compute vorticity in spectral space and transfer to physical space
     vortx_hat=cmplx(0, 1)*(ky*uz_hat-kz*uy_hat)
     vorty_hat=cmplx(0, 1)*(kz*ux_hat-kx*uz_hat)
     vortz_hat=cmplx(0, 1)*(kx*uy_hat-ky*ux_hat)
     call decomp_2d_fft_3d(vortx_hat, vortx)
     call decomp_2d_fft_3d(vorty_hat, vorty)
     call decomp_2d_fft_3d(vortz_hat, vortz)
     
     ! Compute u x \omega in physical space
     do k=xstart(3), xend(3)
        do j=xstart(2), xend(2)
           do i=xstart(1), xend(1)
              velvortx(i, j, k) = uy(i,j,k)*vortz(i,j,k)-uz(i,j,k)*vorty(i,j,k)
              velvorty(i, j, k) = uz(i,j,k)*vortx(i,j,k)-ux(i,j,k)*vortz(i,j,k)
              velvortz(i, j, k) = ux(i,j,k)*vorty(i,j,k)-uy(i,j,k)*vortx(i,j,k)
           enddo
        enddo
     enddo

     call decomp_2d_fft_3d(velvortx, velvortx_hat)
     call decomp_2d_fft_3d(velvorty, velvorty_hat)
     call decomp_2d_fft_3d(velvortz, velvortz_hat)
     velvortx_hat=velvortx_hat/nxnynz
     velvorty_hat=velvorty_hat/nxnynz
     velvortz_hat=velvortz_hat/nxnynz

     ! Dealias modes
     where (ik2>klimit)
        velvortx_hat=0
        velvorty_hat=0
        velvortz_hat=0
     endwhere
     
     if (istep==1) then
        rhsx0=velvortx_hat
        rhsy0=velvorty_hat
        rhsz0=velvortz_hat
     endif

     ux_hat=ux_hat*(1-0.5*visc*dt*k2)+3./2.*dt*velvortx_hat-1./2.*dt*rhsx0
     uy_hat=uy_hat*(1-0.5*visc*dt*k2)+3./2.*dt*velvorty_hat-1./2.*dt*rhsy0
     uz_hat=uz_hat*(1-0.5*visc*dt*k2)+3./2.*dt*velvortz_hat-1./2.*dt*rhsz0

     rhsx=kx*ux_hat+ky*uy_hat+kz*uz_hat
     ux_hat=ux_hat-kx*(rhsx)*k2inv
     uy_hat=uy_hat-ky*(rhsx)*k2inv
     uz_hat=uz_hat-kz*(rhsx)*k2inv

     ux_hat=ux_hat/(1+0.5*visc*dt*k2)
     uy_hat=uy_hat/(1+0.5*visc*dt*k2)
     uz_hat=uz_hat/(1+0.5*visc*dt*k2)

     ! Dealias modes
     where (ik2>klimit)
        ux_hat=0
        uy_hat=0
        uz_hat=0
     endwhere

     rhsx0=velvortx_hat
     rhsy0=velvorty_hat
     rhsz0=velvortz_hat
        
     if (mod(istep, ndump)==0) then
        write(filename, "('uuu',i9.9,'.dat')") istep
        call decomp_2d_write_one(1, ux, filename)
        write(filename, "('vvv',i9.9,'.dat')") istep
        call decomp_2d_write_one(1, uy, filename)
        write(filename, "('www',i9.9,'.dat')") istep
        call decomp_2d_write_one(1, uz, filename)
     endif
     if (mod(istep, noutput)==0) then
        call output_fstep(nx, ny, nz, ux, uy, uz, istep)
     endif
  enddo ! timeloop

  call cpu_time(finish)
  finish=MPI_Wtime()
  if (nrank==0) print '("Time = ",f16.5," seconds.")', finish-start
  
  deallocate(ux, uy, uz)
  deallocate(ux_hat, uy_hat, uz_hat)
  
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_Finalize(ierr)
endprogram main

! This subroutine is used to output velocity in TECPLOT ASCII format
subroutine output_fstep(nx, ny, nz, fx, fy, fz, istep)
  implicit none
  integer :: nx, ny, nz
  integer :: istep
  real(8), dimension(nx, ny, nz) :: fx, fy, fz
  integer :: i, j, k
  real :: dx, dy, dz
  real :: pi
  character(len=16) :: filename

  pi=4.0*atan(1.0)
  dx=pi*2/nx
  dy=dx
  dz=dx

  ! Z plane
  write(filename, "('s2Z',i9.9,'.tec')") istep
  open(unit=100, file=filename, status="replace")

  write(100, *) 'TITLE = "2D Z-slice"'
  write(100, *) 'VARIABLES = "X", "Y", "Z", "U", "V", "W"'
  write(100, "('ZONE I=',i4,1x,', J=',i4,1x,', K=',i4,1x,'F=POINT')") nx, ny, 1

  k=nz/2
  do j=1,ny
     do i=1,nx
        write(100, "(3(f10.5,2x),3(1pe15.7,2x))") &
             (i-0.5)*dx, (j-0.5)*dy, (k-0.5)*dz, &
             fx(i,j,k), fy(i,j,k), fz(i,j,k)
     enddo
  enddo
  close(100)

  ! Y plane
  write(filename, "('s2Y',i9.9,'.tec')") istep
  open(unit=100, file=filename, status="replace")

  write(100, *) 'TITLE = "2D Y-slice"'
  write(100, *) 'VARIABLES = "X", "Y", "Z", "U", "V", "W"'
  write(100, "('ZONE I=',i4,1x,', J=',i4,1x,', K=',i4,1x,'F=POINT')") nx, 1, nz

  j=ny/2
  do k=1, nz
     do i=1,nx
        write(100, "(3(f10.5,2x),3(1pe15.7,2x))") &
             (i-0.5)*dx, (j-0.5)*dy, (k-0.5)*dz, &
             fx(i,j,k), fy(i,j,k), fz(i,j,k)
     enddo
  enddo
  close(100)

  ! X plane
  write(filename, "('s2X',i9.9,'.tec')") istep
  open(unit=100, file=filename, status="replace")

  write(100, *) 'TITLE = "2D X-slice"'
  write(100, *) 'VARIABLES = "X", "Y", "Z", "U", "V", "W"'
  write(100, "('ZONE I=',i4,1x,', J=',i4,1x,', K=',i4,1x,'F=POINT')") 1, ny, nz

  i=nx/2
  do k=1,nz
     do j=1,ny
        write(100, "(3(f10.5,2x),3(1pe15.7,2x))") &
             (i-0.5)*dx, (j-0.5)*dy, (k-0.5)*dz, &
             fx(i,j,k), fy(i,j,k), fz(i,j,k)
     enddo
  enddo
  close(100)

  ! Uncomment to save whole 3D arrays in TECPLOT ASCII format
  ! write(filename, "('sim',i9.9,'.tec')") istep
  ! open(unit=100, file=filename, status="replace")

  ! write(100, *) 'TITLE = "3D data"'
  ! write(100, *) 'VARIABLES = "X", "Y", "Z", "U", "V", "W"'
  ! write(100, "('ZONE I=',i4,1x,', J=',i4,1x,', K=',i4,1x,'F=POINT')") nx, ny, nz

  ! do k=1,nx
  !    do j=1,ny
  !       do i=1,nz
  !          write(100, "(3(f10.5,2x),3(1pe15.7,2x))") &
  !               (i-0.5)*dx, (j-0.5)*dy, (k-0.5)*dz, &
  !               fx(i,j,k), fy(i,j,k), fz(i,j,k)
  !       enddo
  !    enddo
  ! enddo
  ! close(100)
endsubroutine output_fstep

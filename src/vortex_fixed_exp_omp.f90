!    this program calculates the evolution of a vortex sheet
!    in a periodic domain of length l using the approach of 
!    krasny(1986).

!     z=x+iy is the complex conjugate of position vector
!       in this program
!    before running the program, set the environment variable run_dir
!    to be the program running directory. e.g
!    export run_dir=/home/documents/vortex/run/

!   vortex_fixed_exp_mpi.f90 => Simulates a fixed number of vortices
!                               using explicit RK4 time advancement
!                               uses OpenMP for parallel communications
!----------------------------------------------------------------------
    module digit
    integer,parameter :: qp=selected_real_kind(32)
    end module digit
!----------------------------------------------------------------------
    module vortex_data
    use digit
    implicit none
    save
    complex(qp),allocatable,dimension(:) :: f,temp,z,zn,zbar,zbarn
    real(qp),allocatable,dimension(:) :: x,y
    real(qp),allocatable,dimension(:) :: gama
    real(qp) :: delx,delta,ggama,l,eps,delta2,kw
    real(qp),parameter :: pi=4.0_qp*atan(1.0_qp)
    integer :: n
    end module vortex_data    
!----------------------------------------------------------------------
    program vortex
    use digit
    use vortex_data
    implicit none
    integer :: i,j,k,nt,kk,jj,restart,nt_start
    real(qp) :: delt,finalt,c1,c2
    real(qp),allocatable,dimension(:) :: si
    complex(qp),allocatable,dimension(:) :: k1,k2,k3,k4
    character(len=255) :: filename,tempfile,dir_path



    l=1.0_qp     ! length of the domain
    n=32000        ! no. of vortices
    delx=l/n     ! initial spacing between vortices
    n=nint(l/delx)     ! no. of vortices
    ggama=1.0_qp     ! total circulation of the vortex sheet
    eps=0.01_qp     ! radius of vortex blobs (epsilon)
    delta=eps/l     ! the scaled radius of vortex blobs(=epsilon/l)
    delta2=delta*delta
    delt=0.001_qp     ! time step    
    finalt=0.98_qp     ! final time for the simulation
    nt=nint(finalt/delt) ! no. of time steps 
    restart=1 ! if set equal to 0, its a fresh run. if set to 1, then its a restart from nt_start
    nt_start=552    ! starting point for the present calculation - the input file must be 'datant_start' 
!              with the position of the vortices (x_i,y_i) as two columns 
!    set nt_start to zero if it is a fresh start
    kw=2.0_qp*pi/l  ! wave number
    c1=1.0_qp/6.0_qp
    c2=1.0_qp/3.0_qp


    allocate(z(n),zbar(n),zn(n),zbarn(n),gama(n),si(n),x(n),y(n))
    allocate(k1(n),k2(n),k3(n),k4(n),f(n),temp(n))

!    initial conditions section - for both fresh run and restart in serc

    if (restart .eq. 0) then

!    initial conditions
    do i=1,n
    si(i)=(i-1)*delx
    z(i)=cmplx(si(i)+0.01_qp*sin(kw*si(i)),-0.01_qp*sin(kw*si(i))+0.01_qp*sin(4.0_qp*kw*si(i)),qp)
!    z(i)=cmplx(si(i)+0.01_qp*sin(kw*si(i)),-0.01_qp*sin(kw*si(i)),qp)
    enddo
    
    else
    
!    read in initial vortex positions

    write(tempfile,'(a,i0,a)') "data",nt_start
    filename=trim(adjustl(tempfile))

    open(unit=22,defaultfile='/localscratch/aserajsh/fixed_vortex/',file=filename,status='old',action='read')
    do i=1,n
       read(22,300) x(i),y(i)
    enddo
    close(22)
300    format (e42.33e3,e42.33e3)
!300    format (e,e)  ! works only with intel compiler

    do i=1,n
    z(i)=cmplx(x(i),y(i),qp)
    enddo

    endif

    do i=1,n
    gama(i)=ggama*l/n
    enddo


!    call system('gnuplot -p vortex_init_plot.plt')

!       call get_environment_variable("run_dir",tempfile)
!    dir_path=trim(adjustl(tempfile))

!    write and display initial conditions
    open(unit=38,file='/localscratch/aserajsh/fixed_vortex/vortex_init_data',position='append',status='old')
!    open(unit=38,defaultfile=dir_path,file='vortex_init_data',position='append',status='replace')
    do j=1,n
    write(38,300) real(z(j)), aimag(z(j))
    enddo
    close(38)

    zbar=conjg(z)

!    time advancing loop
    do k=nt_start+1,nt 
    
!    rk4 method 

!    deallocate(temp,k1,k2,k3,k4,f,zbarn)
!    allocate(temp(n),k1(n),k2(n),k3(n),k4(n),f(n),zbarn(n))
!     uncomment the above lines for the general case

    temp=zbar   ! temp is the variable that holds y in y'=f(y,t) 
    call rkfun
    k1=delt*f
    temp=zbar+0.5_qp*k1
    call rkfun
    k2=delt*f
    temp=zbar+0.5_qp*k2
    call rkfun
    k3=delt*f
    temp=zbar+k3        
    call rkfun
    k4=delt*f

    zbarn=zbar+c1*k1+c2*(k2+k3)+c1*k4 
    
    open(unit=48,file='/localscratch/aserajsh/fixed_vortex/run_status',position='append',status='old')
!    open(unit=48,defaultfile=dir_path,file='run_status',position='append',status='old')
    write(48,*) 'time step completed=',k*delt 
    close(48)

    zn=conjg(zbarn)
!    call vortex_intro ! enable for the introduction of new vortices
    z=zn 
    zbar=zbarn  !modify to conjg(z) for the general case


!    write to file
    write(tempfile,'(a,i0,a)') "data",k
    filename=trim(adjustl(tempfile))

    open(unit=58,defaultfile='/localscratch/aserajsh/fixed_vortex',file=filename,status='new')
!    open(unit=58,defaultfile=dir_path,file=filename,status='new')
    do j=1,n
    write(58,300) real(z(j)),aimag(z(j))
    enddo
    close(58)

!    call system('gnuplot -p vortex_plot.plt')
    enddo    ! end of time advancement loop

    open(unit=78,file='/localscratch/aserajsh/fixed_vortex/run_status',position='append',status='old')
!    open(unit=78,defaultfile=dir_path,file='run_status',position='append',status='old')
    write(78,*) 'end of computation'
    close(78)

    end program vortex 

!----------------------------------------------------------------------
    subroutine rkfun
!    this subroutine calculates the rhs function for the rk4 method
    use digit
    use vortex_data
    use omp_lib
    implicit none
    integer :: i,j
    real(qp) :: t1,t2,t3    
    complex(qp) :: dtemp,ess

    do i=1,n
    f(i)=cmplx(0.0_qp,0.0_qp,qp)
    enddo


    call omp_set_nested(.true.)

    if(omp_get_nested().eq.1) then
    write(*,*) 'openmp nested - enabled'
    endif

!$omp parallel shared(f) num_threads(16)
!$omp do private(ess) 
    do i=1,n
    ess=cmplx(0.0_qp,0.0_qp,qp)


!$omp parallel num_threads(2)
!$omp do private(dtemp,t1,t2,t3) firstprivate(kw,gama,temp,delta2,n) reduction(+:ess)
    do j=1,n
    if(j.ne.i) then
    dtemp=temp(i)-temp(j)
    t1=sinh(-kw*aimag(dtemp))
    t2=sin(kw*real(dtemp))
    t3=cosh(-kw*aimag(dtemp))-cos(kw*real(dtemp))+(delta2)
    ess=ess+(-0.5_qp*gama(j))*cmplx(t1/t3,t2/t3,qp)
    endif
    enddo
!$omp end do
!$omp end parallel


    f(i)=ess
    enddo
!$omp end do
!$omp end parallel

    end subroutine rkfun
!-----------------------------------------------------------------------    





    

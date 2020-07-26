!    this program calculates the evolution of a vortex sheet
!    in a periodic domain of length l using the approach of 
!    krasny(1986).

!     z=x+iy is the complex conjugate of position vector
!       in this program

!   vortex_fixed_exp_mpi.f90 => Simulates a fixed number of vortices
!                               using explicit RK4 time advancement
!                               uses MPI for parallel communications
!----------------------------------------------------------------------
    module vortex_data
    implicit none
    save
    complex(kind=16),allocatable,dimension(:) :: f,temp,z,zn,zbar,zbarn,fn
    real(kind=16),allocatable,dimension(:) :: x,y
    real(kind=16),allocatable,dimension(:) :: gama
    real(kind=16) :: delx,delta,ggama,l,eps,delta2,kw
    real(kind=16),parameter :: pi=4.0q0*atan(1.0q0)
    integer :: n,nn
    end module vortex_data    
!-----------------------------------------------------------------------
    program vortex
    use vortex_data
    include 'mpif.h'
    implicit none
    integer :: i,j,k,nt,kk,jj,restart,nt_start
    real(kind=16) :: delt,finalt,c1,c2
    real(kind=16),allocatable,dimension(:) :: si
    complex(kind=16),allocatable,dimension(:) :: k1,k2,k3,k4
    character(len=255) :: filename,tempfile

    integer :: ierr,num_procs,proc_num,quad_real,quad_complex
    integer, dimension(mpi_status_size) :: status

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, num_procs, ierr)
    call mpi_comm_rank(mpi_comm_world, proc_num, ierr)

    call mpi_type_contiguous(2,mpi_real8,quad_real,ierr)
    call mpi_type_commit(quad_real,ierr)

    call mpi_type_contiguous(2,mpi_complex16,quad_complex,ierr)
    call mpi_type_commit(quad_complex,ierr)

    n=32000
    nn=nint(real(n)/real(num_procs))
    n=nn*num_procs
    
    l=1.0q0     ! length of the domain
!    n=32000        ! no. of vortices
    delx=l/n     ! initial spacing between vortices
!    n=nint(l/delx)     ! no. of vortices
    ggama=1.0q0     ! total circulation of the vortex sheet
    eps=0.01q0     ! radius of vortex blobs (epsilon)
    delta=eps/l     ! the scaled radius of vortex blobs(=epsilon/l)
    delta2=delta*delta
    delt=0.0001q0     ! time step    
    finalt=0.98q0     ! final time for the simulation
    nt=nint(finalt/delt) ! no. of time steps 
    restart=0 ! if set equal to 0, its a fresh run. if set to 1, then its a restart from nt_start
    nt_start=0    ! starting point for the present calculation - the input file must be 'datant_start' 
!              with the position of the vortices (x_i,y_i) as two columns 
!    set nt_start to zero if it is a fresh start
    kw=2.0q0*pi/l  ! wave number
    c1=1.0q0/6.0q0
    c2=1.0q0/3.0q0

!-----------------------------------------------------------------------------------------------------------
! only the master initializes the data
    if (proc_num == 0) then

    allocate(z(n),zbar(n),zn(n),zbarn(n),gama(n),si(n),x(n),y(n))
    allocate(k1(n),k2(n),k3(n),k4(n),f(n),temp(n),fn(nn))

!    initial conditions section - for both fresh run and restart in serc

    if (restart .eq. 0) then

!    initial conditions
    do i=1,n
    si(i)=(i-1)*delx
    z(i)=cmplx(si(i)+0.01q0*sin(kw*si(i)),-0.01q0*sin(kw*si(i))+0.01q0*sin(4.0q0*kw*si(i)),16)
!    z(i)=cmplx(si(i)+0.01q0*sin(kw*si(i)),-0.01q0*sin(kw*si(i)),16)
    enddo
    
    else
    
!    read in initial vortex positions

    write(tempfile,'(a,i0,a)') "data",nt_start
    filename=trim(adjustl(tempfile))

    open(unit=22,file=filename,status='old',action='read')
    do i=1,n
       read(22,400) x(i),y(i)
    enddo
    close(22)
400    format (e41.32e3,e41.32e3)

    do i=1,n
    z(i)=cmplx(x(i),y(i),16)
    enddo

    endif

    do i=1,n
    gama(i)=ggama*l/n
    enddo


!    write and display initial conditions
    open(unit=38,file='/gpfs/bglscratch/aserajsh/fixed_vortex/vortex_init_data',position='append',status='old')
    do j=1,n
    write(38,300) real(z(j)), aimag(z(j))
    enddo
    close(38)
300    format (e41.32e3,e41.32e3)
    zbar=conjg(z)

    endif  ! the master finishes initializing the data
!-----------------------------------------------------------------------------------------------------------
! workers need only temp(n) and gama(n), fn(nn)

    if (proc_num /= 0) then
    allocate(gama(n),temp(n),fn(nn))    
    endif

!    since this is a fixed_vortex case, we can broadcast the strength of vortices once and for all
    call mpi_bcast(gama,n,quad_real,0,mpi_comm_world,ierr)

!-----------------------------------------------------------------------------------------------------------
!    mpi syntax
!    call mpi_bcast(start,count,datatype,root,comm,ierr)
!    call mpi_send(start,count,datatype,dest,tag,comm,ierr) 
!    call mpi_recv(start,count,datatype,source,tag,comm,status,ierr)    
!    call mpi_gather(sendbuf,sendcnt,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!    call mpi_scatter(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!-----------------------------------------------------------------------------------------------------------
!    time advancing loop
    do k=nt_start+1,nt

!    rk4 method 
    if (proc_num == 0) then
!    deallocate(temp,k1,k2,k3,k4,f,zbarn)
!    allocate(temp(n),k1(n),k2(n),k3(n),k4(n),f(n),zbarn(n))
!     uncomment the above lines for the general case
    temp=zbar   ! temp is the variable that holds y in y'=f(y,t)
    endif

!    send temp vector to all other processors - step 1 of rk4
    call mpi_bcast(temp,n,quad_complex,0,mpi_comm_world,ierr)
    call rkfun(proc_num)
    call mpi_gather(fn,nn,quad_complex,f,nn,quad_complex,0,mpi_comm_world,ierr)
    
    if (proc_num == 0) then
    k1=delt*f
    temp=zbar+0.5q0*k1
    endif

!    send temp vector to all other processors - step 2 of rk4
    call mpi_bcast(temp,n,quad_complex,0,mpi_comm_world,ierr)
    call rkfun(proc_num)
    call mpi_gather(fn,nn,quad_complex,f,nn,quad_complex,0,mpi_comm_world,ierr)
    
    if (proc_num == 0) then
    k2=delt*f
    temp=zbar+0.5q0*k2
    endif

!    send temp vector to all other processors - step 3 of rk4
    call mpi_bcast(temp,n,quad_complex,0,mpi_comm_world,ierr)
    call rkfun(proc_num)
    call mpi_gather(fn,nn,quad_complex,f,nn,quad_complex,0,mpi_comm_world,ierr)
    
    if (proc_num == 0) then
    k3=delt*f
    temp=zbar+k3
    endif


!    send temp vector to all other processors - step 4 of rk4
    call mpi_bcast(temp,n,quad_complex,0,mpi_comm_world,ierr)
    call rkfun(proc_num)
    call mpi_gather(fn,nn,quad_complex,f,nn,quad_complex,0,mpi_comm_world,ierr)
    
    if (proc_num == 0) then
    k4=delt*f
    zbarn=zbar+c1*k1+c2*(k2+k3)+c1*k4 
    
!    open(unit=48,file='/gpfs/bglscratch/aserajsh/fixed_vortex/run_status',position='append',status='old')
!    write(48,*) 'time step completed=',k*delt 
!    close(48)

!    zn=conjg(zbarn)
!!    call vortex_intro ! enable for the introduction of new vortices
!    z=zn 
!    zbar=zbarn  !modify to conjg(z) for the general case

    zbar=zbarn

!    write to file
    write(tempfile,'(a,i0,a)') "data",k
    filename=trim(adjustl(tempfile))

! file is created in the current working directory
    open(unit=58,file=filename,status='new')
    do j=1,n
!    write(58,300) real(z(j)),aimag(z(j))
    write(58,300) real(zbar(j)),-aimag(zbar(j))
    enddo
    close(58)

    endif

    enddo    ! end of time advancement loop

!    if (proc_num == 0) then
!    open(unit=78,file='/gpfs/bglscratch/aserajsh/fixed_vortex/run_status',position='append',status='old')
!    write(78,*) 'end of computation'
!    close(78)
!    endif

    call mpi_type_free(quad_real, ierr)
    call mpi_type_free(quad_complex, ierr)

    call mpi_finalize(ierr)
    end program vortex 

!----------------------------------------------------------------------
    subroutine rkfun(proc_num)
!    this subroutine calculates the rhs function for the rk4 method
    use vortex_data
    implicit none
    integer,intent(in) :: proc_num
    integer :: i,j,istart,iend,itemp
    real(kind=16) :: t1,t2,t3    
    complex(kind=16) :: dtemp,ess

!    do i=1,n
!    f(i)=cmplx(0.0q0,0.0q0,16)
!    enddo

!    fn=cmplx(0.0q0,0.0q0,16)
    istart=(proc_num*nn)+1
    iend=(proc_num*nn)+nn
    itemp=1
    do i=istart,iend
        ess=cmplx(0.0q0,0.0q0,16)
            do j=1,n
                if(j.ne.i) then
                    dtemp=temp(i)-temp(j)
                    t1=sinh(-kw*aimag(dtemp))
                    t2=sin(kw*real(dtemp))
                    t3=cosh(-kw*aimag(dtemp))-cos(kw*real(dtemp))+(delta2)
                    ess=ess+(-0.5q0*gama(j))*cmplx(t1/t3,t2/t3,16)
                endif
            enddo
        fn(itemp)=ess    
        itemp=itemp+1
    enddo

    end subroutine rkfun
!-----------------------------------------------------------------------    





    

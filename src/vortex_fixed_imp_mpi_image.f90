!    this program calculates the evolution of a vortex sheet
!    in a periodic domain of length l using the approach of 
!    krasny(1986).

!     z=x+iy is the complex conjugate of position vector
!   vortex_fixed_imp_mpi_image.f90 => Simulates a fixed number of vortices with an image system
!                                      using implicit time advancement
!                                      uses MPI for parallel communications
!                                      uses the Krasny filter to smoother the precision errors
!----------------------------------------------------------------------
! this module contains data that is shared by the subroutines calling it 
    module vortex_data
    implicit none
    save
    complex(kind=16),allocatable,dimension(:) :: f,temp,z,fn
    real(kind=16),allocatable,dimension(:) :: gama
    real(kind=16) :: delx,delta,ggama,l,eps,delta2,kw
    real(kind=16),parameter :: pi=4.0q0*atan(1.0q0)
    integer :: n,nn,nby2,nnn
    end module vortex_data    
!-----------------------------------------------------------------------
!    main program

    program vortex
    use vortex_data
    include 'mpif.h'
    implicit none
    integer :: i,j,k,nt,it,kj,istart,iend,itemp,tag,m
    complex(kind=16) :: tsum
    real(kind=16) :: delt,delt2,finalt,error,a,b
    real(kind=16),allocatable,dimension(:) :: si,pk
    complex(kind=16),allocatable,dimension(:) :: z_1,z_2,z_fft_temp,z_fft,z1,z2
    character(len=255) :: filename,tempfile

    integer :: ierr,num_procs,proc_num,quad_real,quad_complex
    integer, dimension(mpi_status_size) :: status

!    initialize mpi processes

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, num_procs, ierr)
    call mpi_comm_rank(mpi_comm_world, proc_num, ierr)

!    create mpi data types for quadruple precision real and complex numbers
    call mpi_type_contiguous(2,mpi_real8,quad_real,ierr)
    call mpi_type_commit(quad_real,ierr)

    call mpi_type_contiguous(2,mpi_complex16,quad_complex,ierr)
    call mpi_type_commit(quad_complex,ierr)

!    basic initialization done by all processes

    n=32768        ! number of vortices

! rounding off the number of vortices as integer multiple of the number of processes
    nn=nint(real(n)/real(num_procs))
    n=nn*num_procs
    nby2=n/2
    nnn=nby2/num_procs
    
    l=1.0q0     ! length of the domain
    delx=l/nby2     ! initial spacing between vortices
    ggama=1.0q0     ! total circulation of the vortex sheet
    eps=0.01q0     ! radius of vortex blobs (epsilon)
    delta=eps/l     ! the scaled radius of vortex blobs(=epsilon/l)
    delta2=delta*delta
    delt=0.001q0     ! time step
    delt2=0.5q0*delt        
    finalt=0.98q0     ! final time for the simulation
    nt=nint(finalt/delt) ! no. of time steps 
    kw=2.0q0*pi/l  ! wave number
    tag=1
!-----------------------------------------------------------------------------------------------------------
! only the master initializes the vortex position data and other variables of size n
    if (proc_num == 0) then

    allocate(z(n),gama(n),si(n),z_1(n),z_2(n))
    allocate(f(n),temp(n),fn(nn),z_fft_temp(nnn),z_fft(nby2),pk(nby2))
    allocate(z1(nby2),z2(nby2))

!    initial conditions

!    initialize the vortex position
    do i=1,nby2
    si(i)=(i-1)*delx
    z(i)=cmplx(si(i)+0.01q0*sin(kw*si(i)),-0.01q0*sin(kw*si(i)),16)
    enddo

    do i=nby2+1,n
    si(i)=(i-nby2-1)*delx
    z(i)=cmplx(si(i)+0.01q0*sin(kw*si(i)),-0.4q0-0.01q0*sin(kw*si(i)+pi),16)
    enddo

!    initialize the vorticity distribution    
    do i=1,nby2
    gama(i)=ggama*l/nby2
    enddo
    do i=nby2+1,n
    gama(i)=-ggama*l/nby2
    enddo


!    write and display initial conditions
    open(unit=38,file='/gpfs/bglscratch/aserajsh/fixed_vortex_image_32000/vortex_init_data',position='append',status='old')
    do j=1,n
    write(38,300) real(z(j)), aimag(z(j))
    enddo
    close(38)
300    format (e41.32e3,e41.32e3)

    endif  ! the master finishes initializing the data
!-----------------------------------------------------------------------------------------------------------
! workers need only temp(n) and gama(n), fn(nn)

    if (proc_num /= 0) then
    allocate(gama(n),temp(n),fn(nn),z(n),z_fft_temp(nnn),z_fft(nby2))
    allocate(z1(nby2),z2(nby2))    
    endif

!    since this is a fixed_vortex case, we can broadcast the strength of vortices once and for all
    call mpi_bcast(gama,n,quad_real,0,mpi_comm_world,ierr)

!-----------------------------------------------------------------------------------------------------------
!    mpi syntax
!    call mpi_bcast(start,count,datatype,root,comm,ierr)
!    call mpi_send(start,count,datatype,dest,tag,comm,ierr) 
!    call mpi_recv(start,count,datatype,source,tag,comm,status,ierr)    
!    call mpi_gather(sendbuf,sendcnt,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!    call mpi_scatter(sendbuf,sendcnt,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!-----------------------------------------------------------------------------------------------------------
!    time advancing loop
    do k=1,nt

    if (proc_num == 0) then
!    if(k == 601) then
!    do i=1,n
!    z(i)=z(i)+cmplx(0.0q0,0.01q0*sin(4.0q0*kw*si(i)),16)
!    enddo
!    endif
    temp=z   ! temp is the variable that holds y in y'=f(y,t)
    endif

    call mpi_bcast(temp,n,quad_complex,0,mpi_comm_world,ierr)
    call comp_vel(proc_num)        ! compute the point vortices velocity - f vector
    call mpi_gather(fn,nn,quad_complex,f,nn,quad_complex,0,mpi_comm_world,ierr)


    if (proc_num == 0) then
        z_1=z+delt2*f    ! half time step    
        temp=z+delt*f    ! predictor time step
    endif
    
    call mpi_bcast(temp,n,quad_complex,0,mpi_comm_world,ierr)
    call comp_vel(proc_num)    
    call mpi_gather(fn,nn,quad_complex,f,nn,quad_complex,0,mpi_comm_world,ierr)

    if (proc_num == 0) then
        temp=z_1+delt2*f    ! one time step
    endif

    it = 1
    do                                  ! error test and iteration

        if (proc_num == 0) then
            z_2=temp
        endif
        
        call mpi_bcast(temp,n,quad_complex,0,mpi_comm_world,ierr)
        call comp_vel(proc_num)
        call mpi_gather(fn,nn,quad_complex,f,nn,quad_complex,0,mpi_comm_world,ierr)

        if (proc_num == 0) then 
            temp = z_1+delt2*f               ! iterate v(t+dt)
!            error = maxval(abs(z_2-temp))
            error = maxval(abs(real(z_2)-real(temp))+abs(aimag(z_2)-aimag(temp)))
        endif
            it = it + 1
        call mpi_bcast(error,1,quad_real,0,mpi_comm_world,ierr)        

            if ((error .lt. 0.001q0*delta) .or. (it.gt.30)) exit
    enddo    

    if(proc_num == 0) then
    z=temp
    endif

    call mpi_barrier(mpi_comm_world,ierr)    
!-----------------------------------------------------------------
!    krasny filter starts here
    if (tag == 1) then

    if (proc_num == 0) then

    do i=1,nby2
        z1(i)=z(i)-si(i)
    enddo
    
    do i=nby2+1,n
        z2(i-nby2)=z(i)-cmplx(si(i),-0.4q0,16)
    enddo

    endif

!    filtering the original vortex sheet(z1) begins here
    call mpi_bcast(z1,nby2,quad_complex,0,mpi_comm_world,ierr)


    istart=(proc_num*nnn)+1
    iend=(proc_num*nnn)+nnn
    itemp=1

    do m=istart,iend
        tsum=cmplx(0.0q0,0.0q0,16)
        do j=1,nby2
        kj=(m-1)*(j-1)
        tsum=tsum+z1(j)*cmplx(cos(kw*kj/nby2),sin(-kw*kj/nby2),16)
        enddo
    z_fft_temp(itemp)=tsum
    itemp=itemp+1
    enddo

    call mpi_gather(z_fft_temp,nnn,quad_complex,z_fft,nnn,quad_complex,0,mpi_comm_world,ierr)

    if(proc_num == 0) then
    z_fft=z_fft/nby2
    pk=log(abs(z_fft))
    error=-69.0776q0
    j=0
    do i=1,nby2
        if(pk(i) .lt. error) then
            z_fft(i)=cmplx(0.0q0,0.0q0,16)
            j=j+1
        endif
    enddo

    endif

    call mpi_bcast(z_fft,nby2,quad_complex,0,mpi_comm_world,ierr)

    istart=(proc_num*nnn)+1
    iend=(proc_num*nnn)+nnn
    itemp=1
    do m=istart,iend
        tsum=cmplx(0.0q0,0.0q0,16)
        do j=1,nby2
        kj=(m-1)*(j-1)
        tsum=tsum+z_fft(j)*cmplx(cos(kw*kj/nby2),sin(kw*kj/nby2),16)
        enddo
    z_fft_temp(itemp)=tsum
    itemp=itemp+1
    enddo

    call mpi_gather(z_fft_temp,nnn,quad_complex,z1,nnn,quad_complex,0,mpi_comm_world,ierr)

    if(proc_num == 0) then

        do i=1,nby2
            z1(i)=z1(i)+si(i)
        enddo
    endif
    
!    filtering the original vortex sheet (z1) ends here

!    filtering the image vortex sheet (z2) begins here

    call mpi_bcast(z2,nby2,quad_complex,0,mpi_comm_world,ierr)


    istart=(proc_num*nnn)+1
    iend=(proc_num*nnn)+nnn
    itemp=1

    do m=istart,iend
        tsum=cmplx(0.0q0,0.0q0,16)
        do j=1,nby2
        kj=(m-1)*(j-1)
        tsum=tsum+z2(j)*cmplx(cos(kw*kj/nby2),sin(-kw*kj/nby2),16)
        enddo
    z_fft_temp(itemp)=tsum
    itemp=itemp+1
    enddo

    call mpi_gather(z_fft_temp,nnn,quad_complex,z_fft,nnn,quad_complex,0,mpi_comm_world,ierr)

    if(proc_num == 0) then
    z_fft=z_fft/nby2
    pk=log(abs(z_fft))
    error=-69.0776q0
    do i=1,nby2
        if(pk(i) .lt. error) then
            z_fft(i)=cmplx(0.0q0,0.0q0,16)
            j=j+1
        endif
    enddo
    if(j .eq. 0) then
        tag=0
        open(unit=38,file='/gpfs/bglscratch/aserajsh/fixed_vortex_image_32000/filter_status',position='append',status='old')
        write(38,*) 'switching off the filter', k
        close(38)
    endif

    endif

    call mpi_bcast(tag,1,mpi_integer,0,mpi_comm_world,ierr)    
    call mpi_bcast(z_fft,nby2,quad_complex,0,mpi_comm_world,ierr)

    istart=(proc_num*nnn)+1
    iend=(proc_num*nnn)+nnn
    itemp=1
    do m=istart,iend
        tsum=cmplx(0.0q0,0.0q0,16)
        do j=1,nby2
        kj=(m-1)*(j-1)
        tsum=tsum+z_fft(j)*cmplx(cos(kw*kj/nby2),sin(kw*kj/nby2),16)
        enddo
    z_fft_temp(itemp)=tsum
    itemp=itemp+1
    enddo

    call mpi_gather(z_fft_temp,nnn,quad_complex,z2,nnn,quad_complex,0,mpi_comm_world,ierr)
    if(proc_num == 0) then
        do i=nby2+1,n
            z2(i-nby2)=z2(i-nby2)+cmplx(si(i),-0.4q0,16)
        enddo
    endif
!    filtering the image vortex sheet (z2) ends here

!    reconstructing the full variable z including both vortex sheets
    if(proc_num == 0) then
    do i=1,nby2
        z(i)=z1(i)
    enddo
    do i=nby2+1,n
        z(i)=z2(i-nby2)
    enddo
    endif
!    end reconstructing the full variable z including both vortex sheets

    endif
!    krasny filter ends here
!-------------------------------------------------------------------    
    call mpi_barrier(mpi_comm_world,ierr)    

    if(proc_num == 0) then

!    write to file
        write(tempfile,'(a,i0,a)') "data",k
        filename=trim(adjustl(tempfile))

        !     file is created in the current working directory
        open(unit=58,file=filename,status='new')
        do j=1,n
        write(58,300) real(z(j)),aimag(z(j))
        enddo
        close(58)

    endif

    enddo    ! end of time advancement loop

    call mpi_type_free(quad_real, ierr)
    call mpi_type_free(quad_complex, ierr)

    call mpi_finalize(ierr)
    end program vortex 

!----------------------------------------------------------------------
    subroutine comp_vel(proc_num)
!    this subroutine calculates the point vortices velocity from their position
    use vortex_data
    implicit none
    integer,intent(in) :: proc_num
    integer :: i,j,istart,iend,itemp
    real(kind=16) :: t1,t2,t3    
    complex(kind=16) :: dtemp,ess

!    each worker process calculates a part of the f vector, namely the fn vector
    istart=(proc_num*nn)+1
    iend=(proc_num*nn)+nn
    itemp=1
    do i=istart,iend
        ess=cmplx(0.0q0,0.0q0,16)
            do j=1,n
                if(j.ne.i) then
                    dtemp=temp(i)-temp(j)
                    t1=-sinh(kw*aimag(dtemp))
                    t2=sin(kw*real(dtemp))
                    t3=cosh(kw*aimag(dtemp))-cos(kw*real(dtemp))+(delta2)
                    ess=ess+(0.5q0*gama(j))*cmplx(t1/t3,t2/t3,16)
                endif
            enddo
        fn(itemp)=ess    
        itemp=itemp+1
    enddo

    end subroutine comp_vel
!----------------------------------------------------------------------



    

    program energy
    include 'mpif.h'
    implicit none
    integer,parameter :: wp=8
    integer :: i,j,k,m,n,nx,ny,nxn,nxi,grid_size,grid_size_n,nby2
    real(kind=wp) :: dxi,delx,dely,lx,ly,total_e_field,total_e_spect
    complex(kind=wp),allocatable,dimension(:,:) :: vel,vel_temp,grid,grid_n
    complex(kind=wp),allocatable,dimension(:) :: vortex,fft_vel_corr
    complex(kind=wp) :: dtemp,ess
    real(kind=wp) :: t1,t2,t3,kw,delta,delta2,e_0,im,x1,y1
    real(kind=wp),parameter :: pi=4.0_wp*atan(1.0_wp)
    real(kind=wp),allocatable,dimension(:) :: gama,x,y,e1
    character(len=255) :: filename,tempfile

    real(kind=wp),allocatable,dimension(:) :: xi,c,cc,e,temp
    complex(kind=wp),allocatable,dimension(:,:,:) :: vel_corr

    integer :: ierr,num_procs,proc_num
    integer, dimension(mpi_status_size) :: status

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, num_procs, ierr)
    call mpi_comm_rank(mpi_comm_world, proc_num, ierr)


    n=32768        ! total number of vortices
    nby2=n/2
    nx=1024        
    nxn=nx/num_procs    ! no. of x points for each processor
    ny=1024
    grid_size=nx*ny
    grid_size_n=nxn*ny

    ly=3.0_wp
    lx=1.0_wp
    delx=lx/(nx-1)
    dely=ly/(ny-1)
    
    kw=2.0_wp*pi/lx     ! wave number
    delta=0.01_wp/lx     ! the scaled radius of vortex blobs(=epsilon/l)
    delta2=delta*delta

    if (proc_num == 0) then
    allocate(vel(nx,ny),vel_temp(nxn,ny),grid(nx,ny),vortex(n),gama(n),x(n),y(n),e1(ny),grid_n(nxn,ny))
    else
    allocate(gama(n),vortex(n),vel_temp(nxn,ny),grid_n(nxn,ny))
    endif
    
    if (proc_num == 0) then

!     create the xi array which is the correlation distance
!    if the xi array matches the x-coordinate array, then it is easier to calculate correlations
    nxi=nx        ! no. of points in the array xi
    dxi=lx/(nxi-1)        ! jump in the distance
    allocate(xi(nxi),vel_corr(2,nx,ny),c(nxi),cc(ny),fft_vel_corr(nxi),e(nxi),temp(2*nxi))    
    xi(1)=0.449_wp
    do i=2,nxi
        xi(i)=xi(i-1)+dxi
    enddo

!    form the 2-d grid in the region (0:lx,-0.5*ly:0.5*ly)
    do i=1,nx
        do j=1,ny
            grid(i,j) = cmplx(0.4490_wp+(i-1)*delx,(-0.5_wp*ly)+(j-1)*dely,wp)
        enddo
    enddo
!    strength of vortices initialization
    do i=1,nby2
        gama(i)=1.0_wp/nby2
    enddo
    do i=nby2+1,n
        gama(i)=-1.0_wp/nby2
    enddo
    endif
!-----------------------------------------------------------------------------------------------------------
!    mpi syntax
!    call mpi_bcast(start,count,datatype,root,comm,ierr)
!    call mpi_send(start,count,datatype,dest,tag,comm,ierr) 
!    call mpi_recv(start,count,datatype,source,tag,comm,status,ierr)    
!    call mpi_gather(sendbuf,sendcnt,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!    call mpi_scatter(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!-----------------------------------------------------------------------------------------------------------


    call mpi_bcast(gama,n,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_scatter(grid,grid_size_n,mpi_complex16,grid_n,grid_size_n,mpi_complex16,0,mpi_comm_world,ierr)

!    do m=980,1,-1    ! loop going over each data file
    m=980

    if (proc_num == 0) then
!    read in the vortex position data
    write(tempfile,'(a,i0,a)') "data",m
    filename=trim(adjustl(tempfile))

    open(unit=22,file=filename,status='old',action='read',form='formatted')
    do i=1,n
       read(22,400) x(i),y(i)
    enddo
    close(22)

    do i=1,n
    vortex(i)=cmplx(x(i),y(i),wp)
    enddo
    endif    !! end of task 1 for proc 0

!    send vortices position to all procs

    call mpi_bcast(vortex,n,mpi_complex16,0,mpi_comm_world,ierr)

!    everyone calculates a part of the velocity field    
    do i=1,nxn
        do j=1,ny
                   ess=cmplx(0.0_wp,0.0_wp,wp)
            do k=1,n
                dtemp=grid_n(i,j)-vortex(k)
                t1=-sinh(kw*aimag(dtemp))
                t2=sin(kw*real(dtemp))
                t3=cosh(kw*aimag(dtemp))-cos(kw*real(dtemp))+(delta2)
                ess=ess+(0.5_wp*gama(k))*cmplx(t1/t3,t2/t3,wp)
            enddo
                    vel_temp(i,j)=ess
        enddo
    enddo

!    create the image vortex sheet
    
!    im=-0.2_wp

!    do i=1,n
!    x1=real(vortex(i),wp)
!    y1=aimag(vortex(i))
!    y1=2.0_wp*im-y1
!    vortex(i)=cmplx(x1,y1,wp)
!    gama(i)=-gama(i)
!    enddo

!    everyone calculates a part of the image velocity field    
!    do i=1,nxn
!        do j=1,ny
!                   ess=cmplx(0.0_wp,0.0_wp,wp)
!            do k=1,n
!                dtemp=grid_n(i,j)-vortex(k)
!                t1=-sinh(kw*aimag(dtemp))
!                t2=sin(kw*real(dtemp))
!                t3=cosh(kw*aimag(dtemp))-cos(kw*real(dtemp))+(delta2)
!                ess=ess+(0.5_wp*gama(k))*cmplx(t1/t3,t2/t3,wp)
!            enddo
!                       vel_temp(i,j)=vel_temp(i,j)+ess
!        enddo
!    enddo
    

!    reconstruct the full velocity field
    call mpi_gather(vel_temp,grid_size_n,mpi_complex16,vel,grid_size_n,mpi_complex16,0,mpi_comm_world,ierr)



    if (proc_num == 0) then
!    write the velocity field to a file
    if (m == 980) then
    write(tempfile,'(a,i0,a)') "vel_field",m
    filename=trim(adjustl(tempfile))

    open(unit=48,file=filename,status='new')
    do i=1,nx
        do j=1,ny
            write(48,300) real(grid(i,j)),aimag(grid(i,j)),real(vel(i,j)),aimag(vel(i,j))
        enddo
    enddo
    close(48)
    endif


!    estimate the total kinetic energy contained in the flow field
    do j=1,ny
        e1(j)=0.0_wp
        do i=2,nx    
            e1(j)=e1(j)+delx*0.5_wp*(abs(vel(i,j))*abs(vel(i,j))+abs(vel(i-1,j))*abs(vel(i-1,j)))
        enddo
    enddo
    total_e_field=0.0_wp
    do j=2,ny
        total_e_field=total_e_field+dely*0.5_wp*(e1(j)+e1(j-1))
    enddo
    total_e_field=total_e_field/2.0_wp
    
!----------------------------------------------------------------------------------
!     calculate the velocity correlation function
    vel_corr(1,:,:)=vel(:,:)
    vel_corr(2,:,:)=vel(:,:)

    c=0.0_wp

    do k=1,nxi

    if (k .eq. 1) goto 10
        do i=1,nx 
            vel_corr(2,i,:)=vel_corr(1,mod(i,nx)+1,:)
        enddo
        
10    continue

        do j=1,ny
            cc(j)=0.0_wp
            do i=2,nx
                t1=real(vel_corr(2,i,j))*real(vel(i,j))+aimag(vel_corr(2,i,j))*aimag(vel(i,j))
                t2=real(vel_corr(2,i-1,j))*real(vel(i-1,j))+aimag(vel_corr(2,i-1,j))*aimag(vel(i-1,j))
                cc(j)=cc(j)+delx*0.5_wp*(t1+t2)
            enddo
        enddo
        
        do j=2,ny
            c(k)=c(k)+dely*0.5_wp*(cc(j)+cc(j-1))
        enddo
        
        vel_corr(1,:,:)=vel_corr(2,:,:)
    enddo
    c=c/ly

!    find the fast fourier transform of the velocity correlation function c(xi)

    temp=0.0_wp
    j=1
    do i=1,2*nxi-1,2
        temp(i)=c(j)
        j=j+1
    enddo
    call fft(temp,nxi,-1)
    temp=temp/nxi


    j=1
    do k=1,2*nxi-1,2
    fft_vel_corr(j)=cmplx(temp(k),temp(k+1),wp)
    j=j+1
    enddo

    e=abs(fft_vel_corr)/2.0_wp
    
!    estimate the total energy contained in the spectrum
    total_e_spect=0.0_wp
    do j=2,nxi
        total_e_spect=total_e_spect+0.5_wp*(2.0_wp/lx)*(e(j)+e(j-1))
    enddo

    open(unit=58,file='/gpfs/bglscratch/aserajsh/energy/total_energy',status='old',position='append')
    write(58,*) m,', ',total_e_field,', ',total_e_spect
    close(58)

    write(tempfile,'(a,i0,a)') "energy_spec",m
    filename=trim(adjustl(tempfile))

    open(unit=68,file=filename,status='new')
    do k=1,nxi
        write(68,*) e(k)
    enddo
    close(68)
!----------------------------------------------------------------------------------
    endif    !    end of task for proc 0

    call mpi_barrier(mpi_comm_world,ierr)

!    enddo        !    end of reading each data file

300    format (e15.6e3,e15.6e3,e15.6e3,e15.6e3)
400    format (e42.32e3,e42.32e3)


    call mpi_finalize(ierr)
    end program energy

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

    subroutine fft(data,nn,isign)
    integer,parameter :: wp=8
    integer isign,nn
    real(kind=wp) data(2*nn)
!    replaces data(1:2*nn) by its discrete fourier transform, if isign is input as 1; or replaces
!    data(1:2*nn) by nn times its inverse discrete fourier transform, if isign is input as −1.
!    data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
!    must be an integer power of 2 (this is not checked for!).
    integer i,istep,j,m,mmax,n
    real(kind=wp) tempi,tempr
    real(kind=wp) theta,wi,wpi,wpr,wr,wtemp
    real(kind=wp),parameter :: pi=4.0_wp*atan(1.0_wp)
    n=2*nn
    j=1
    do i=1,n,2
!    this is the bit-reversal section of the routine.
        if(j.gt.i)then
            tempr=data(j)
!            exchange the two complex numbers.
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
        endif
        m=n/2
1        if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        goto 1
        endif
        j=j+m
    enddo

    mmax=2    !    here begins the danielson-lanczos section of the routine.
2    if (n.gt.mmax) then    !    outer loop executed log2 nn times.
    istep=2*mmax
    theta=2.0_wp*pi/(isign*mmax)    !    initialize for the trigonometric recurrence.
    wpr=-2.0_wp*sin(0.5_wp*theta)**2
    wpi=sin(theta)
    wr=1.0_wp
    wi=0.0_wp
        do m=1,mmax,2     !    here are the two nested inner loops.

            do i=m,n,istep
                j=i+mmax    !    this is the danielson-lanczos formula:
                tempr=wr*data(j)-wi*data(j+1)
                tempi=wr*data(j+1)+wi*data(j)
                data(j)=data(i)-tempr
                data(j+1)=data(i+1)-tempi
                data(i)=data(i)+tempr
                data(i+1)=data(i+1)+tempi
            enddo 
            wtemp=wr    !    trigonometric recurrence.
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
        enddo 
        mmax=istep
        goto 2    !    not yet done.
    endif        !    all done.
    return
    end subroutine fft
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

















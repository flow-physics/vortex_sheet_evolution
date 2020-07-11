!    This program calculates the evolution of a vortex sheet
!    in a periodic domain of length L using the approach of 
!    Krasny(1986).

!     Z=X+iY is the complex conjugate of position vector

!    Author: Rajesh Venkatesan
!----------------------------------------------------------------------
! This module contains shared data for the simulation
    module vortex_data
    implicit none
    save
    complex(kind=16),allocatable,dimension(:) :: F,temp,Z,FN
    real(kind=16),allocatable,dimension(:) :: X,Y
    real(kind=16),allocatable,dimension(:) :: gama
    real(kind=16) :: delx,delta,ggama,L,eps,delta2,kw
    real(kind=16),PARAMETER :: pi = 4.0Q0*atan(1.0Q0)
    integer :: N,NN
    end module vortex_data    
!-----------------------------------------------------------------------
!    MAIN program

    program vortex
    use vortex_data
    use mpi
    implicit none
    integer :: i,j,k,ntime,it
    real(kind=16) :: delt,delt2,final_time,error,A,B
    real(kind=16),allocatable,dimension(:) :: SI
    complex(kind=16),allocatable,dimension(:) :: Z_1,Z_2
    character(len=255) :: filename,tempfile

    integer :: ierr,num_procs,proc_num,quad_real,quad_complex
    integer, dimension(MPI_STATUS_SIZE) :: status

!    Initialize MPI processes

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

!    Create MPI data types for quadruple precision real and complex numbers
    call MPI_TYPE_CONTIGUOUS(2,MPI_REAL8,quad_real,ierr)
    call MPI_TYPE_COMMIT(quad_real,ierr)

    call MPI_TYPE_CONTIGUOUS(2,MPI_COMPLEX16,quad_complex,ierr)
    call MPI_TYPE_COMMIT(quad_complex,ierr)

!    Basic initialization done by all processes

    N = 32000        ! Number of vortices

! Rounding off the number of vortices as integer multiple of the number of processes
    NN = NINT(real(N)/real(num_procs))
    N  = NN*num_procs
    
    L      = 1.0Q0             ! Length of the domain
    delx   = L/N               ! Initial spacing between vortices
    ggama  = 1.0Q0             ! Total circulation of the vortex sheet
    eps    = 0.01Q0            ! Radius of vortex blobs (epsilon)
    delta  = eps/L             ! The scaled radius of vortex blobs(=epsilon/L)
    delta2 = delta*delta
    delt   = 0.00001Q0         ! Time Step
    delt2  = 0.5Q0*delt        
    final_time = 0.98Q0        ! Final time for the simulation
    ntime  = NINT(final_time/delt) ! No. of time steps 
    kw     = 2.0Q0*pi/L        ! Wave number

!-----------------------------------------------------------------------------------------------------------
! Only the Master initializes the vortex position data and other variables of size N
    if (proc_num == 0) then
        allocate(Z(N),gama(N),SI(N),X(N),Y(N),Z_1(N),Z_2(N))
        allocate(F(N),temp(N),FN(NN))

    !    Initial conditions

    !    Initialize the vortex position
        do i=1,N
            SI(i)=(i-1)*delx
            Z(i)=cmplx(SI(i)+0.01Q0*sin(kw*SI(i)),-0.01Q0*sin(kw*SI(i))+0.01Q0*sin(4.0Q0*kw*SI(i)),16)
        end do

    !    Initialize the vorticity distribution    
        do i=1,N
        gama(i)=ggama*L/N
        end do


    !    Write and display initial conditions
        open(unit=38,file='./vortex_init_data',POSITION='APPEND',status='OLD')
            do j=1,N
                write(38,300) real(Z(j)), aimag(Z(j))
            end do
        close(38)
300     format (E41.32E3,E41.32E3)
    end if  ! The Master finishes initializing the data
!-----------------------------------------------------------------------------------------------------------
! Workers need only temp(N) and gama(N), FN(NN)

    if (proc_num /= 0) then
        allocate(gama(N),temp(N),FN(NN))    
    end if

!    Since this is a case with fixed number of vortices, We can broadcast the strength of vortices once and for all
    call MPI_BCAST(gama,N,quad_real,0,MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------------------------------------------------------
!    MPI Syntax
!    call MPI_BCAST(start,count,datatype,root,comm,ierr)
!    call MPI_SEND(start,count,datatype,dest,tag,comm,ierr) 
!    call MPI_RECV(start,count,datatype,source,tag,comm,status,ierr)    
!    call MPI_GATHER(sendbuf,sendcnt,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!    call MPI_SCATTER(sendbuf,sendcnt,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
!-----------------------------------------------------------------------------------------------------------
!    TIME ADVANCING LOOP
    do k=1,ntime
        if (proc_num == 0) then
            temp=Z   ! temp is the variable that holds y in y'=f(y,t)
        end if

        call MPI_BCAST(temp,N,quad_complex,0,MPI_COMM_WORLD,ierr)
        call COMP_VEL(proc_num)        ! compute the point vortices velocity - F vector
        call MPI_GATHER(FN,NN,quad_complex,F,NN,quad_complex,0,MPI_COMM_WORLD,ierr)


        if (proc_num == 0) then
            Z_1=Z+delt2*F    ! half time step    
            temp=Z+delt*F    ! predictor time step
        end if
        
        call MPI_BCAST(temp,N,quad_complex,0,MPI_COMM_WORLD,ierr)
        call COMP_VEL(proc_num)    
        call MPI_GATHER(FN,NN,quad_complex,F,NN,quad_complex,0,MPI_COMM_WORLD,ierr)

        if (proc_num == 0) then
            temp=Z_1+delt2*F    ! one time step
        end if

        it = 1
        do                                  ! error test and iteration

            if (proc_num == 0) then
                Z_2=temp
            end if
            
            call MPI_BCAST(temp,N,quad_complex,0,MPI_COMM_WORLD,ierr)
            call COMP_VEL(proc_num)
            call MPI_GATHER(FN,NN,quad_complex,F,NN,quad_complex,0,MPI_COMM_WORLD,ierr)

            if (proc_num == 0) then 
                temp = Z_1+delt2*F               ! iterate v(t+dt)
    !            error = maxval(ABS(Z_2-temp))
                error = maxval(abs(real(Z_2)-real(temp))+abs(aimag(Z_2)-aimag(temp)))
            end if
                it = it + 1
            call MPI_BCAST(error,1,quad_real,0,MPI_COMM_WORLD,ierr)        

            if ((error .lt. 0.001Q0*delta) .OR. (it .gt. 30)) exit
        end do    
        
    !    Master writes the data for the finished time step
        if(proc_num == 0) then
            Z=temp
        !    Shift vortices going out of the domain
            do i=1,N
                if(real(Z(i)) .gt. 1.0Q0) THEN
                    A=real(Z(i))
                    B=aimag(Z(i))
                    A=A-1.0Q0
                    Z(i)=cmplx(A,B,16)
                end if
                if(real(Z(i)) .lt. 0.0Q0) THEN
                    A=real(Z(i))
                    B=aimag(Z(i))
                    A=A+1.0Q0
                    Z(i)=cmplx(A,B,16)
                end if            
            end do
        !    End shifting vortices
            if (mod(k,10) .eq. 0) then 
            
                !    Write to file
                write(tempfile,'(A,I0,A)') "data",k
                filename=trim(adjustl(tempfile))

                !     File is created in the current working directory
                open(unit=58,file=filename,status='NEW')
                do j=1,N
                    write(58,300) real(Z(j)),aimag(Z(j))
                end do
                close(58)
            end if
        end if
    end do    ! End of Time advancement loop

    call MPI_TYPE_FREE(quad_real, ierr)
    call MPI_TYPE_FREE(quad_complex, ierr)

    call MPI_FINALIZE(ierr)
    end program vortex 

!----------------------------------------------------------------------
    subroutine COMP_VEL(proc_num)
!    This subroutine calculates the point vortices velocity from their position
    use vortex_data
    implicit none
    integer,intent(in) :: proc_num
    integer :: i,j,istart,iend,itemp
    real(kind=16) :: T1,T2,T3    
    complex(kind=16) :: DTEMP,ESS

!    Each worker process calculates a part of the F vector, namely the FN vector
    istart=(proc_num*NN)+1
    iend=(proc_num*NN)+NN
    itemp=1
    do i = istart,iend
        ESS=cmplx(0.0Q0,0.0Q0,16)
        do j=1,N
            if(j .ne. i) THEN
                DTEMP=temp(i)-temp(j)
                T1=-sinh(kw*aimag(DTEMP))
                T2=sin(kw*real(DTEMP))
                T3=COSH(kw*aimag(DTEMP))-COS(kw*real(DTEMP))+(delta2)
                ESS=ESS+(0.5Q0*gama(j))*cmplx(T1/T3,T2/T3,16)
            end if
        end do
        FN(itemp)=ESS    
        itemp=itemp+1
    end do

    end subroutine COMP_VEL
!-----------------------------------------------------------------------    





    

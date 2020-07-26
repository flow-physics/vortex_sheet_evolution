    program Hamiltonian
    include 'mpif.h'
    integer :: i,j,k,n,nn,nt
    integer :: istart,iend,itemp
    real(kind=16) :: delta,delta2,kw,l,eps,h1,hsum,ess
    real(kind=16),allocatable,dimension(:) :: h,x,y,htemp,htotal
    real(kind=16),parameter :: pi=4.0q0*atan(1.0q0)


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

    l=1.0q0 
    n=32768 
    nn=n/num_procs
    kw=2.0q0*pi/l 
    eps=0.01q0     ! radius of vortex blobs (epsilon)
    delta=eps/l     ! the scaled radius of vortex blobs(=epsilon/l)
    delta2=delta*delta 
    nt=980

    if (proc_num == 0) then
    allocate(h(nt),x(n),y(n),htemp(nn),htotal(num_procs))
    endif

    if (proc_num /= 0) then
    allocate(x(n),y(n),htemp(nn))
    endif
300    format (e41.32e3,e41.32e3)
400    format (e10.3e2,e41.32e3)

    do k=1,nt

    if(proc_num == 0) then
    write(tempfile,'(a,i0,a)') "data",k
    filename=trim(adjustl(tempfile))
    open(unit=22,file=filename,status='old',form='formatted',action='read')
    do i=1,n
       read(22,300) x(i),y(i)
    enddo
    close(22)
    endif

    call mpi_bcast(x,n,quad_real,0,mpi_comm_world,ierr)
    call mpi_bcast(y,n,quad_real,0,mpi_comm_world,ierr)


    istart=(proc_num*nn)+1
    iend=(proc_num*nn)+nn
    itemp=1

    do i=istart,iend
    ess=0.0q0
    do j=1,n
    if (j.ne.i) then
    ess=ess+log(cosh(kw*(y(i)-y(j)))-cos(kw*(x(i)-x(j)))+delta2)
    endif
    enddo
    htemp(itemp)=ess
    itemp=itemp+1
    enddo
    hsum=sum(htemp)

!    call mpi_gather(sendbuf,sendcnt,sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)
    call mpi_gather(hsum,1,quad_real,htotal,1,quad_real,0,mpi_comm_world,ierr)    

!    write to run status file 
    if(proc_num == 0) then
    h1=sum(htotal)
    h(k)=(-h1/(4.0q0*pi*n))

    open(unit=48,file='/gpfs/bglscratch/aserajsh/energy/run_status',position='append',status='old')
    write(48,*) 'time=',k*0.001q0,h(k)
    close(48)
    endif
!    end writing run status

    enddo ! end of time advancement loop


    if(proc_num == 0) then

    open(unit=58,file='/gpfs/bglscratch/aserajsh/energy/hamiltonian',position='append',status='old')
    do j=1,nt
    write(58,400) j*0.001q0,h(j)
    enddo
    close(58)

    endif

    call mpi_type_free(quad_real, ierr)
    call mpi_finalize(ierr)

    end program Hamiltonian

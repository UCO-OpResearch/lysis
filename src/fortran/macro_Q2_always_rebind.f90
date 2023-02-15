program macrolysis

!! BRAD 2023-01-24: This code has been modified in the following ways:
!!                  - Data folder is relative to git repository root
!!                  - Data is stored in subfolders based on expCode
!!                  - Data file codes are now set globally from the (in/out)FileCode variables
!!                  - "passerby molecule" is corrected
!!                  - More console output during runs
!!                  - Parameters can be set from the command line
!!                  - Parameters have been moved to the top of the file
!!

!Runs the macroscale model in a clot with 72.7 nm diameter fibers and pore size. 1.0135 uM. FB conc. = 8.8 uM. This code allows tPA that is forced to unbind by plasmin-mediated degradation of fibrin to immediately rebind. THIS CODE ACCOUNTS FOR MICRO RUNS IN WHICH 50,000 OR 10,000 INDEPENDENT SIMULATIONS WERE DONE. CHANGE LINE 16 (nummicro=) to 500 or 100 depending on if 50,000 or 10,000 micro runs were completed. This code also computes mean first passage time, so it's like a combination of macro_Q2.f90 and macro_Q2_mfpt.f90.
implicit none
character(15)   :: expCode
character(80)   :: inFileCode
character(80)   :: outFileCode

integer  :: N=93!93  !# of lattice nodes in one row in the horizontal direction
integer  :: F=121!121 !71 !81  !# of lattice nodes in one column in the vertical direction
integer  :: Ffree=29!29 !3 !13 !1st node in vertical direction containing fibers. so if Ffree=10, then rows 1-9
                               !have no fibers, there's one more row of fiber-free planar veritcal edges, and then
                               !the row starting with the Ffree-th (e.g. 10th) vertical node is a full row of fibers 
integer  :: stats= 10 !! BRAD 2023-01-04: 10
integer  :: M=43074 !total number of tPA molecules: 21588 is Colin's [tPA]=0.3 nM; 43074 is Colin's [tPA]=0.6 nM; 86148 is Colin's [tPA]=1.2 nM;
integer  :: tf=20*60 !! BRAD 2023-01-06: 20*60!15*60 !final time in sec

integer  :: nummicro=500 !if the number of microscale runs was 50,000, take nummicro=500; if it was 10,000, take nummicro=100
!!!CHANGES MADE FOR FORCED UNBINDING/DIFFUSION/REBINDING:
double precision :: kon = 0.1 !0.1 !tPA binding rate. units of inverse (micromolar*sec). MAKE SURE THIS MATCHES MICROSCALE RUN VALUE!!!!
double precision :: frac_forced =0.0852 !0.5143!0.0054!0.0852 !fraction of times tPA was forced to unbind in microscale model. MAKE SURE THIS MATCHES MICROSCALE RUN VALUE!!!!
double precision :: avgwait = 27.8 !2.78 !27.8 !measured in seconds, this is the average time a tPA molecule stays bound to fibrin. It's 1/koff. For now I'm using 27.8 to be 1/0.036, the value in the absence of PLG

double precision :: q=0.2d+00      !0.2             !q is the probability of moving. Make sure it is small enough that we've converged
double precision :: delx= 1.0135d-04 !10**(-4)             !pore size (distance between nodes), measured in centimeters
double precision :: Diff= 5.0d-07 !5*10**(-7)           !diffusion coefficient, measured in cm**2/s
!! BRAD 2023-01-08: Does this need to be a float, or can it be an integer?
!! BRITT:           Has been an integer for years and will probably stay that way
integer          :: bs = 427              !concentration of binding sites in micromolar
double precision :: dist=1.0862d+00 !microns because distance between nodes is 1.0135 micron and diameter of 1 fiber is 0.0727 micron

integer  :: seed=-2137354075


integer  :: num
integer  :: enoFB !the last edge number without fibrin
integer  :: backrow !defines the first fiber number in the back row of the clot. To calculate mean first passage time, I will record the first time that each tPA molecule diffuses to a fiber with edge number backrow or greater
double precision :: tstep  !4/6/11 CHANGED THIS TO (12*Diff) FROM (8*Diff). SEE WRITTEN NOTES 4/6/11 FOR WHY
double precision :: num_t            !number of timesteps


integer  :: i, istat
integer  :: j, ij, newindex
integer  :: k
integer  :: ii 
integer  :: Nsave, Ninteger, nplt, cNsave
integer  :: count, countij
integer  :: count2, countc, countcolr4 
integer  :: count66 
integer  :: Ntot2
integer  :: colr2, colr4 
integer  :: z
integer  :: numPLi
double precision     :: t

double precision     :: t_bind
double precision     :: percent2, percent4
double precision     :: rmicro, ttPA
double precision     :: r, r1, r3, r2, r4


character(80) :: filetype,formatted
character(80) :: filename1
character(80) :: filename2
character(90) :: filename3
character(95) :: filename4
character(80) :: filename6


integer*1, dimension(:, :), allocatable  :: closeneigh
integer, dimension(:,:), allocatable    :: endpts
integer, dimension(:,:), allocatable    :: neighborc
integer, dimension(8)      :: temp_neighborc
integer, dimension(:,:), allocatable   :: V
double precision, dimension(:), allocatable      :: init_state
double precision, dimension(2)         :: p
double precision, dimension(2)         :: pfit
double precision, dimension(101)       :: CDFtPA, CDFlys
double precision, dimension(101)       :: tsec1, tseclys
double precision, dimension(:), allocatable      :: degrade
double precision, dimension(:), allocatable     :: t_degrade
double precision, dimension(:), allocatable        :: t_leave
double precision, dimension(:), allocatable        :: t_wait
double precision, dimension(:), allocatable        :: bind
integer, dimension(:), allocatable             :: Nsavevect 
integer                              :: name1, name2

logical       :: isBinary = .True.      ! flag for binary output
integer       :: degunit = 20
integer       :: Vunit = 21
integer       :: V2unit = 22
integer       :: Nunit = 23
integer       :: tunit = 24
character(80) :: degfile       ! degradation vector
character(80) :: Vfile
character(80) :: V2file
character(80) :: Nfile
character(80) :: tfile

!stuff for the random # generator. I need to use the file kiss.o when I compile in order for this to work, and kiss.o
!is obtained by compiling the file kiss.c by doing "cc -c kiss.c".
integer :: kiss32, mscw, state(4), old_state(4), ui
double precision :: uf, urcw1
external :: mscw, kiss32, urcw1

double precision, dimension(:), allocatable         :: rvect
double precision, dimension(:,:), allocatable :: degnext
integer, dimension(:,:), allocatable :: Vedgenext
integer, dimension(:,:), allocatable :: Vboundnext
integer, dimension(:), allocatable  :: ind
double precision, dimension(:), allocatable  :: place
double precision, dimension(:), allocatable  :: degold
double precision, dimension(:), allocatable :: tsave
integer  :: zero1
integer, dimension(:,:), allocatable  :: front
integer, dimension(:), allocatable  :: firstdeg
integer, dimension(:), allocatable  :: deglast
integer, dimension(:,:), allocatable  :: lastmove
integer  :: fdeg
integer  :: first0
integer, dimension(:,:), allocatable  :: move
integer  :: temp
integer  :: lasti
integer, dimension(:,:), allocatable  :: plotstuff, totmove,time2plot
double precision, dimension(:,:), allocatable  :: plotstuff2
integer       :: moveunit = 25
integer       :: lastmoveunit = 26
integer       :: plotunit = 27
integer       :: degnextunit = 28
integer       :: Venextunit = 29
integer       :: Vbdnextunit = 30
integer       :: mfptunit = 42
character(80) :: movefile
character(80) :: lastmovefile
character(80) :: plotfile
character(80) :: degnextfile
character(80) :: Venextfile
character(80) :: Vbdnextfile
character(80) :: mfptfile

integer, dimension(:), allocatable  :: intact2
integer  :: countintact2, lenintact2, counth, countv, countpv, countmacrounbd, countmicrounbd
integer  :: jj, iplt, yplace, x2, xplace, y1, y2, yvplace, xvplace, imod
integer  :: jplt, kplt, kjplt, vertplace, Vyvert, Vy1, xVedgeplace, Vedgeplace, Vx 
integer, dimension(:,:), allocatable  :: X1plot, Y1plot
integer, dimension(:,:), allocatable  :: X2plot, Y2plot
integer, dimension(:), allocatable  :: Xvplot, Yvplot
double precision, dimension(:,:), allocatable  :: bdtPA, freetPA
double precision, dimension(:,:), allocatable :: lysismat !(100,100) if only did 10,000 micro runs, (500,100) if did 50,000
integer, dimension(100)  :: lenlysismat
integer  :: r400
integer  :: x1unit = 31
integer  :: x2unit = 32
integer  :: y1unit = 33
integer  :: y2unit = 34
integer  :: xvunit = 35
integer  :: yvunit = 36
integer  :: tPAbdunit = 37
integer  :: tPAfreeunit = 38
integer  :: cbindunit = 39
integer  :: cindunit = 40
integer  :: bind1unit = 41
character(80)  :: x1file
character(80)  :: x2file
character(80)  :: y1file
character(80)  :: y2file
character(80)  :: xvfile
character(80)  :: yvfile
character(80)  :: tPAbdfile
character(80)  :: tPAfreefile
character(80)  :: cbindfile
character(80)  :: cindfile
character(80)  :: bind1file
integer :: countbind, countindep
integer, dimension(:,:), allocatable  :: countbindV, countindepV, bind1V
integer, dimension(:), allocatable :: bind1
integer, dimension(:), allocatable :: forcedunbdbydeg
double precision, dimension(:), allocatable :: mfpt !vector I'll use to save the first passage times of each tPA molecule
integer, dimension(:), allocatable :: yesfpt !vector of 1's and 0's to let me know if the particular tPA molecule has already hit the back edge of the clot or not

!! BRAD 2023-01-06:
integer(8) :: total_regular_moves
integer(8) :: total_restricted_moves
integer :: total_binds
integer :: degraded_fibers
integer :: reached_back_row
real :: time_begin, time_end
real :: rounded_time
real :: degraded_percent
real :: reached_back_row_percent
double precision :: last_degrade_time
real :: temp_len_lysis_mat

integer :: t_degrade_unit = 102
character(80) :: t_degrade_file
integer :: m_location_unit = 103
character(80) :: m_location_file
integer :: m_bound_unit = 104
character(80) :: m_bound_file
!integer :: m_bind_time_unit = 105

logical :: all_fibers_degraded
logical :: most_molecules_passed

integer :: cmd_count, param_i, param_len, param_val_len, cmd_status, io_status
character(80) :: param, param_value

cmd_count = command_argument_count ()
write (*,*) 'number of command arguments = ', cmd_count

param_i = 0
do while (param_i<cmd_count)
    param_i = param_i+1
    call get_command_argument (param_i, param, param_len, cmd_status)
    if (cmd_status .ne. 0) then
        write (*,*) ' get_command_argument failed: status = ', cmd_status, ' arg = ', param_i
        stop
    end if
    write (*,*) 'command arg ', param_i, ' = ', param (1:param_len)
    param_i = param_i+1
    call get_command_argument (param_i, param_value, param_val_len, cmd_status)
    if (cmd_status .ne. 0) then
        write (*,*) ' get_command_argument failed: status = ', cmd_status, ' arg = ', param_i
        stop
    end if
    write (*,*) 'command arg ', param_i, ' = ', param_value (1:param_val_len)
    select case (param(3:param_len))
        case('expCode')
            expCode = param_value (1:param_val_len)
            write (*,*) 'Setting expCode = ', expCode
        case ('inFileCode')
            inFileCode = param_value (1:param_val_len)
            write (*,*) 'Setting inFileCode = ', inFileCode
        case ('outFileCode')
            outFileCode = param_value (1:param_val_len)
            write (*,*) 'Setting outFileCode = ', outFileCode
        case ('N')
            read(param_value,*,iostat=io_status)  N
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting N = ', N
        case ('F')
            read(param_value,*,iostat=io_status)  F
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting F = ', F
        case ('Ffree')
            read(param_value,*,iostat=io_status)  Ffree
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting Ffree = ', Ffree
        case ('stats')
            read(param_value,*,iostat=io_status)  stats
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting stats = ', stats
        case ('M')
            read(param_value,*,iostat=io_status)  M
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting M = ', M
        case ('tf')
            read(param_value,*,iostat=io_status)  tf
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting tf = ', tf
        case ('nummicro')
            read(param_value,*,iostat=io_status)  nummicro
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting nummicro = ', nummicro
        case ('kon')
            read(param_value,*,iostat=io_status)  kon
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting kon = ', kon
        case ('frac_forced')
            read(param_value,*,iostat=io_status)  frac_forced
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting frac_forced = ', frac_forced
        case ('avgwait')
            read(param_value,*,iostat=io_status)  avgwait
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting avgwait = ', avgwait
        case ('q')
            read(param_value,*,iostat=io_status)  q
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting q = ', q
        case ('delx')
            read(param_value,*,iostat=io_status)  delx
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting delx = ', delx
        case ('Diff')
            read(param_value,*,iostat=io_status)  Diff
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting Diff = ', Diff
        case ('bs')
            read(param_value,*,iostat=io_status)  bs
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting bs = ', bs
        case ('dist')
            read(param_value,*,iostat=io_status)  dist
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting dist = ', dist
        case ('seed')
            read(param_value,*,iostat=io_status)  seed
            if (io_status .ne. 0) then
                write (*,*) 'String conversion error'
                stop
            end if
            write (*,*) 'Setting seed = ', seed
        case default
            write (*,*) 'Unrecognized parameter'
            stop
    end select
end do

write (*,*) 'command line processed'


num=(2*N-1)*F+N*(F-1)
enoFB=(3*N-1)*(Ffree-1) !the last edge number without fibrin
backrow=num-(2*N-1)+1 !defines the first fiber number in the back row of the clot. To calculate mean first passage time, I will record the first time that each tPA molecule diffuses to a fiber with edge number backrow or greater
tstep=q*delx**2/(12*Diff)  !4/6/11 CHANGED THIS TO (12*Diff) FROM (8*Diff). SEE WRITTEN NOTES 4/6/11 FOR WHY
num_t=tf/tstep            !number of timesteps

allocate (closeneigh(num,num))
allocate (endpts(2,num))
allocate (neighborc(8,num))
allocate (V(2,M))
allocate (init_state(enoFB))
allocate (degrade(num))
allocate (t_degrade(num))
allocate (t_leave(M))
allocate (t_wait(M))
allocate (bind(M))
allocate (Nsavevect(stats))
allocate (rvect(M))
allocate (degnext(tf+1,num))
allocate (Vedgenext(tf+1,M))
allocate (Vboundnext(tf+1,M))
allocate (ind(F-1))
allocate (place(F-1))
allocate (degold(num))
allocate (tsave(tf+1))
allocate (front(tf,N))
allocate (firstdeg(N))
allocate (deglast(N))
allocate (lastmove(N,stats))
allocate (move(N,N))
allocate (plotstuff(N,N), totmove(N,N),time2plot(N,N))
allocate (plotstuff2(N,N))
allocate (intact2(num))
allocate (X1plot(2,F*(N-1)), Y1plot(2,F*(N-1)))
allocate (X2plot(2,N*(F-1)), Y2plot(2,N*(F-1)))
allocate (Xvplot(N*F), Yvplot(N*F))
allocate (bdtPA(2,M), freetPA(2,M))
allocate (lysismat(nummicro,100))!(100,100) if only did 10,000 micro runs, (500,100) if did 50,000
allocate (countbindV(stats,tf), countindepV(stats,tf), bind1V(stats,tf))
allocate (bind1(num))
allocate (forcedunbdbydeg(M))
allocate (mfpt(M)) !vector I'll use to save the first passage times of each tPA molecule
allocate (yesfpt(M))  !vector of 1's and 0's to let me know if the particular tPA molecule has already hit the back edge of the clot or not


!! BRAD END


if( isBinary ) then
    !filetype = 'unformatted' !if you compile with gfortran or f95
    filetype = 'binary'      !if you compile with ifort
else
    filetype = 'formatted'
end if
write(*,*)' filetype=',filetype


write(*,*)' N=',N
write(*,*)' F=',F
write(*,*)' Ffree=',Ffree
write(*,*)' num=',num
write(*,*)' M=',M
write(*,*)' obtained using code macro_Q2_always_rebind.f90 on data ',expCode
!write(*,*)'fraction of time tPA is forced to unbind',frac_forced



! Initialize the Random Number Generator

ui = kiss32()
uf = urcw1()

!! BRAD 2023-01-26:
if (seed == 0) seed = mscw()
!seed= 1884637428
!seed = -2137354075
!seed = 5784279
write(*,*)' seed=',seed


state(1) = 129281
state(2) = 362436069
state(3) = 123456789
state(4) = seed
call set_kiss32(state)

call get_kiss32(state)

closeneigh=0
neighborc=0

!the corner vertical (3-D) edges
closeneigh(1,2) = 4
closeneigh(1,2*N) = 4

closeneigh(2*N-1,2*N-1-1) = 4
closeneigh(2*N-1,2*N-1+N) = 4

closeneigh((3*N-1)*(F-1)+1,(3*N-1)*(F-1)+1+1) = 4
closeneigh((3*N-1)*(F-1)+1,(3*N-1)*(F-1)+1-N) = 4

closeneigh((3*N-1)*(F-1)+1+2*(N-1),(3*N-1)*(F-1)+1+2*(N-1)-1) = 4
closeneigh((3*N-1)*(F-1)+1+2*(N-1),(3*N-1)*(F-1)+1+2*(N-1)-(2*N-1)) = 4

!the bottom and top left-most and right-most horizontal edges
closeneigh(2,2*N) = 2
closeneigh(2,2*N+1) = 2
closeneigh(2,1) = 2
closeneigh(2,3) = 2

closeneigh(2*N-2,2*N-1) = 2
closeneigh(2*N-2,2*N-3) = 2
closeneigh(2*N-2,2*N-2+N) = 2
closeneigh(2*N-2,2*N-2+N+1) = 2

closeneigh((3*N-1)*(F-1)+2,(3*N-1)*(F-1)+2-1) = 2
closeneigh((3*N-1)*(F-1)+2,(3*N-1)*(F-1)+2+1) = 2
closeneigh((3*N-1)*(F-1)+2,(3*N-1)*(F-1)+2-N) = 2
closeneigh((3*N-1)*(F-1)+2,(3*N-1)*(F-1)+2-(N+1)) = 2

closeneigh((3*N-1)*(F-1)+2+2*(N-2),(3*N-1)*(F-1)+2+2*(N-2)+1) = 2
closeneigh((3*N-1)*(F-1)+2+2*(N-2),(3*N-1)*(F-1)+2+2*(N-2)-1) = 2
closeneigh((3*N-1)*(F-1)+2+2*(N-2),(3*N-1)*(F-1)+2+2*(N-2)-(2*N-1)) = 2
closeneigh((3*N-1)*(F-1)+2+2*(N-2),(3*N-1)*(F-1)+2+2*(N-2)-(2*N-2)) = 2

!right and left bottom-most and top-most vertical planar edges
closeneigh(2*N,1) = 2
closeneigh(2*N,2*N+N) = 2
closeneigh(2*N,2) = 2
closeneigh(2*N,2*N+N+1) = 2

closeneigh(2*N+N-1,2*N+N-1-N) = 2
closeneigh(2*N+N-1,2*N+N-1+2*N-1) = 2
closeneigh(2*N+N-1,2*N+N-1-(N+1)) = 2
closeneigh(2*N+N-1,2*N+N-1+N+(N-2)) = 2

closeneigh((3*N-1)*(F-2)+2*N,(3*N-1)*(F-2)+2*N+N) = 2
closeneigh((3*N-1)*(F-2)+2*N,(3*N-1)*(F-2)+2*N-(2*N-1)) = 2
closeneigh((3*N-1)*(F-2)+2*N,(3*N-1)*(F-2)+2*N+N+1) = 2
closeneigh((3*N-1)*(F-2)+2*N,(3*N-1)*(F-2)+2*N-(2*N-1)+1) = 2

closeneigh((3*N-1)*(F-2)+2*N+(N-1),(3*N-1)*(F-2)+2*N+(N-1)+(2*N-1)) = 2
closeneigh((3*N-1)*(F-2)+2*N+(N-1),(3*N-1)*(F-2)+2*N+(N-1)-N) = 2
closeneigh((3*N-1)*(F-2)+2*N+(N-1),(3*N-1)*(F-2)+2*N+(N-1)+(2*N-1)-1) = 2
closeneigh((3*N-1)*(F-2)+2*N+(N-1),(3*N-1)*(F-2)+2*N+(N-1)-N-1) = 2

!the bottom and top rows of vertical edges on the lattice (not including
!the left-most and right-most edges)
do j=2,N-1
    closeneigh(1+2*(j-1),1+2*(j-1)-1) = 2
    closeneigh(1+2*(j-1),1+2*(j-1)+1) = 2
    closeneigh(1+2*(j-1),2*N-1+j) = 4

    closeneigh(1+2*(j-1)+(3*N-1)*(F-1),1+2*(j-1)+(3*N-1)*(F-1)-1) = 2
    closeneigh(1+2*(j-1)+(3*N-1)*(F-1),1+2*(j-1)+(3*N-1)*(F-1)+1) = 2
    closeneigh(1+2*(j-1)+(3*N-1)*(F-1),(3*N-1)*(F-2)+2*N+(j-1)) = 4

end do

!the left and right columns of vertical edges (not including the top-most
!and bottom-most edges)
do i=2,F-1
    closeneigh((3*N-1)*(i-1)+1,(3*N-1)*(i-1)+1-N) = 2
    closeneigh((3*N-1)*(i-1)+1,(3*N-1)*(i-1)+1+(2*N-1)) = 2
    closeneigh((3*N-1)*(i-1)+1,(3*N-1)*(i-1)+1+1) = 4

    closeneigh((3*N-1)*(i-1)+1+2*(N-1),(3*N-1)*(i-1)+1+2*(N-1)-(2*N-1)) = 2
    closeneigh((3*N-1)*(i-1)+1+2*(N-1),(3*N-1)*(i-1)+1+2*(N-1)+N) = 2
    closeneigh((3*N-1)*(i-1)+1+2*(N-1),(3*N-1)*(i-1)+1+2*(N-1)-1) = 4
end do


!the bottom and top rows of horizontal edges on the lattice (not including
!the left-most and right-most edges)
do j=2,N-2
    closeneigh(2+(j-1)*2,2+(j-1)*2-1) = 2
    closeneigh(2+(j-1)*2,2+(j-1)*2+1) = 2
    closeneigh(2+(j-1)*2,2*N-1+j) = 2
    closeneigh(2+(j-1)*2,2*N-1+j+1) = 2

    closeneigh(2+(j-1)*2+(3*N-1)*(F-1),2+(j-1)*2+(3*N-1)*(F-1)-1) = 2
    closeneigh(2+(j-1)*2+(3*N-1)*(F-1),2+(j-1)*2+(3*N-1)*(F-1)+1) = 2
    closeneigh(2+(j-1)*2+(3*N-1)*(F-1),(3*N-1)*(F-2)+2*N+j-1) = 2
    closeneigh(2+(j-1)*2+(3*N-1)*(F-1),(3*N-1)*(F-2)+2*N+j) = 2
end do
    
!the left and right columns of vertical planar edges on the lattice (not including
!the top-most and bottom-most edges)
do i=2,F-2
    closeneigh(2*N+(i-1)*(3*N-1),2*N+(i-1)*(3*N-1)+N) = 2
    closeneigh(2*N+(i-1)*(3*N-1),2*N+(i-1)*(3*N-1)-(2*N-1)) = 2
    closeneigh(2*N+(i-1)*(3*N-1),2*N+(i-1)*(3*N-1)+N+1) = 2
    closeneigh(2*N+(i-1)*(3*N-1),2*N+(i-1)*(3*N-1)-(2*N-1)+1) = 2

    closeneigh(3*N-1+(i-1)*(3*N-1),3*N-1+(i-1)*(3*N-1)-N) = 2
    closeneigh(3*N-1+(i-1)*(3*N-1),3*N-1+(i-1)*(3*N-1)+(2*N-1)) = 2
    closeneigh(3*N-1+(i-1)*(3*N-1),3*N-1+(i-1)*(3*N-1)-N-1) = 2
    closeneigh(3*N-1+(i-1)*(3*N-1),3*N-1+(i-1)*(3*N-1)+(2*N-1)-1) = 2
end do
    

!finally, do all the remaining edges (i.e. the edges that do not have ghost points on the boundary):

do j=2,N
    do i=2,F-1
        !horizontal egdes

        closeneigh(2+(j-2)*2+(i-1)*(3*N-1),2+(j-2)*2+(i-1)*(3*N-1)+1) = 2
        closeneigh(2+(j-2)*2+(i-1)*(3*N-1),2+(j-2)*2+(i-1)*(3*N-1)-1) = 2
        closeneigh(2+(j-2)*2+(i-1)*(3*N-1),2+(j-2)*2+(i-1)*(3*N-1)-(2*N-1)+N-j) = 1
        closeneigh(2+(j-2)*2+(i-1)*(3*N-1),2+(j-2)*2+(i-1)*(3*N-1)-(2*N-1)+N-j+1) = 1
        closeneigh(2+(j-2)*2+(i-1)*(3*N-1),(3*N-1)*(i-1)+2*N+(j-1)-1) = 1
        closeneigh(2+(j-2)*2+(i-1)*(3*N-1),(3*N-1)*(i-1)+2*N+(j-1)) = 1
    end do
end do

do i=1,F-1
    do j=2,N-1
        !vertical planar edges

        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*(j-1)) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*j) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*(j-1)+(3*N-1)) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*j+(3*N-1)) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),2*N+(j-1)+(i-1)*(3*N-1)-(2*N-1)+(j-1)) = 2
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),2*N+(j-1)+(i-1)*(3*N-1)-(2*N-1)+(j-1)+(3*N-1)) = 2
  
    end do
end do

do i=2,F-1
    do j=2,N-1
        !vertical edges

        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),1+(j-1)*2+(i-1)*(3*N-1)+1) = 2
        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),1+(j-1)*2+(i-1)*(3*N-1)-1) = 2
        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),2*N+(j-2)+(i-2)*(3*N-1)+1) = 2
        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),2*N+(j-2)+(i-2)*(3*N-1)+1+(3*N-1)) = 2
  
    end do
end do

!!Create a 8 x num matrix of neighboring edges. Column number corresponds to edge number, and the row entries represent 
!!neighboring edges
!Create 2 4 x num matrices of neighboring edges. Column number corresponds to edge number, and the row entries represent 
!neighboring edges

do j=1,num
    countc=0
    do i=1,num
        if(closeneigh(j,i)==1) then
            countc=countc+1
            neighborc(countc,j)=i
        elseif(closeneigh(j,i)==2) then
            countc=countc+1
            neighborc(countc,j)=i
            countc=countc+1
            neighborc(countc,j)=i
        elseif(closeneigh(j,i)==4) then
            countc=countc+1
            neighborc(countc,j)=i
            countc=countc+1
            neighborc(countc,j)=i
            countc=countc+1
            neighborc(countc,j)=i
            countc=countc+1
            neighborc(countc,j)=i
        end if
    end do
end do


! read in the data from the micro model, which we obtained from /micro.f90
! READ IN VECTORS FROM MATLAB 



open(200,file=ADJUSTL('data/' // expCode // '/tPAleave' // inFileCode))
do i=1,101
    read(200,*)CDFtPA(i)
end do
close(200)
write(*,*)'read tPAleave.dat'

open(300,file=ADJUSTL('data/' // expCode // '/tsectPA' // inFileCode))
do i=1,101
    read(300,*)tsec1(i)
end do
close(300)
write(*,*)'read tsectPA.dat'



!lysismat_PLG2_tPA01_Q2.dat is a matrix with column corresponding to bin number (1-100) and with entries
!equal to the lysis times obtained in that bin. an entry of 6000 means lysis didn't happen.
!lysismat(:,1)=the first column, i.e. the lysis times for the first 100 (or 500 if we did 50,000 micro runs) tPA leaving times
OPEN(unit=201,FILE=ADJUSTL('data/' // expCode // '/lysismat' // inFileCode))
do i=1,nummicro  !100 if only did 10,000 micro runs, 500 if did 50,000
    READ(201,*)(lysismat(i,ii),ii=1,100)
enddo
close(201)

!lenlysisvect_PLG2_tPA01_Q2.dat saves the first row entry in each column of lysismat_PLG2_tPA01_Q2.dat that lysis
!did not occur, i.e. the first entry there's a 6000
OPEN(unit=202,FILE=ADJUSTL('data/' // expCode // '/lenlysisvect' // inFileCode))
do i=1,100
    READ(202,*)lenlysismat(i)
end do
close(202)


lastmove=0

!The edges that don't contain fibrin are those with edge number <= to enoFB
!enoFB=(3*N-1)*(Ffree-1)

write(*,*)'enoFB=',enoFB

    !!!!!COMMENTED OUT BELOW ON 5/16/16 BECAUSE I DON'T USE THIS DATA IN ANY POST-PROCESSING
    !write(degnextfile,'(57a)') 'degnext_tPA425_PLG2_tPA01_into_and_along_Q2.dat'
    !write(Venextfile,'(59a)') 'Vedgenext_tPA425_PLG2_tPA01_into_and_along_Q2.dat'
    !write(Vbdnextfile,'(57a)') 'Vbdnext_tPA425_PLG2_tPA01_into_and_along_Q2.dat'
    !write(cbindfile,'(57a)') 'numbind_tPA425_PLG2_tPA01_into_and_along_Q2.dat'
    !write(cindfile,'(57a)') 'numindbind_tPA425_PLG2_tPA01_into_and_along_Q2.dat'
    !write(bind1file,'(57a)') 'bind_tPA425_PLG2_tPA01_into_and_along_Q2.dat'
    open(degunit,file=ADJUSTL('data/' // expCode // '/deg' // outFileCode),form=filetype)
    open(Nunit,file=ADJUSTL('data/' // expCode // '/Nsave' // outFileCode),form=filetype)
    open(tunit,file=ADJUSTL('data/' // expCode // '/tsave' // outFileCode),form=filetype)
    open(moveunit,file=ADJUSTL('data/' // expCode // '/move' // outFileCode),form=filetype)
    open(lastmoveunit,file=ADJUSTL('data/' // expCode // '/lastmove' // outFileCode),form=filetype)
    open(plotunit,file=ADJUSTL('data/' // expCode // '/plot' // outFileCode),form=filetype)
    open(mfptunit,file=ADJUSTL('data/' // expCode // '/mfpt' // outFileCode),form=filetype)

!! BRAD 2023-01-21:
    open(t_degrade_unit,file=ADJUSTL('data/' // expCode // '/f_deg_time' // outFileCode),form=filetype)
    open(m_location_unit,file=ADJUSTL('data/' // expCode // '/m_loc' // outFileCode),form=filetype)
    open(m_bound_unit,file=ADJUSTL('data/' // expCode // '/m_bound' // outFileCode),form=filetype)
!   open(m_bind_time_unit,file=ADJUSTL('data/' // expCode // '/m_bind_t' // outFileCode),form=filetype)


    !!!!!COMMENTED OUT BELOW ON 5/16/16 BECAUSE I DON'T USE THIS DATA IN ANY POST-PROCESSING
    !open(degnextunit,file=degnextfile,form=filetype)
    !open(Venextunit,file=Venextfile,form=filetype)
    !open(Vbdnextunit,file=Vbdnextfile,form=filetype)
    !open(cbindunit,file=cbindfile,form=filetype)
    !open(cindunit,file=cindfile,form=filetype)
    !open(bind1unit,file=bind1file,form=filetype)


!initialize variables for MFPT calculation out here b/c we only do this on the first run
!yesfpt=0 !initialize yesfpt to 0, and change individual entries to 1's when that tPA molecule hits the back row of the clot
!mfpt=0

!!!!!! DO "STATS" RUNS OF THE MACRO MODEL
stats_loop: do istat=1,stats

    write(*,*)' run number=',istat

    count=0
    t=0
    count2=0
    countcolr4=0
    count66=0
    countbind=0
    countindep=0
    Nsave=0
    cNsave=0
    bind1=0

!! BRAD 2023-01-04:
    total_binds = 0
    total_regular_moves = 0
    total_restricted_moves = 0
    degraded_fibers = 0
    reached_back_row = 0
    last_degrade_time = 0
    yesfpt=0 !initialize yesfpt to 0, and change individual entries to 1's when that tPA molecule hits the back row of the clot
    mfpt=0
    all_fibers_degraded = .False.
    most_molecules_passed = .False.

    !Initialize vectors to 0
    degrade  =0.0d+00         !vector of degradation state of each edge. 0=not degraded, -t=degraded at time t
    t_degrade=0.0d+00         !vector of the degradation times of each edge
    t_leave  =0.0d+00         !vector of the tPA leaving times for each molecule
    degnext  =0
    Vedgenext = 0
    Vboundnext = 0
    totmove = 0
    time2plot = 0
    bind = 0.0d+00             !vector of the binding times for each tPA

    degrade(1:enoFB) = -1.0  !set the undegradable fibers' degradation status to -1


    write(*,*)' q=',q
    write(*,*)' delx=',delx
    write(*,*)' tstep=',tstep
    write(*,*)' tf=',tf
    write(*,*)' num_t=', num_t
    write(*,*)' kon=',kon
    write(*,*)' bs=',bs
    write(*,*)' avgwait=',avgwait
    write(*,*)' frac_forced=',frac_forced

    !Initial distribution and boundedness of tPA molecules

    !V is matrix of edge loaction and state of boundedness for each tPA. first column is edge
    !molecule is on, second column is 0 if unbound, 1 if bound. start with all molecules unbound

    do i=1,M
        V(2,i)=0
    end do

    do ii=1,enoFB
        init_state(ii) = dble(ii) / dble(enoFB) !make a vector that's the length of one row of lattice and scale so
                                                       !that prob. of being on any edge is equal
    end do

    call vurcw1(rvect,M)

    !use the random numbers to decide where we start the tPAs
    do i=1,M
        if (0.le.rvect(i).and.rvect(i).le.init_state(1)) V(1,i)=1 !init_entry(i)=1
        do j=1,enoFB
            if (init_state(j).lt.rvect(i).and.rvect(i).le.init_state(j+1)) then
                !init_entry(i)=j+1
                V(1,i)=j+1
            end if
        end do
        !V(i,1) = init_entry(i) !molecule i starts on site init_entry(i) so I only allow tPA to start on an edge
                          !that's in the first row of my lattice
    end do
    !    write(*,*)' V=',V  !for debugging 3/31/10



    write(degunit) degrade(:)
    write(tunit) t

!! BRAD 2023-01-21:
    write(t_degrade_unit) t_degrade(:)
    write(m_location_unit) V(1,:)
    write(m_bound_unit) V(2,:)
    Nsave = 10

    write(*,*)' save as deg',outFileCode


    Vedgenext(1,:)=V(1,:)
    Vboundnext(1,:)=V(2,:)
    degnext(1,:)=degrade(:)
    tsave(1) = t



    !Now do the main part of the code - looping over time


!! BRAD 2023-01-06:
    CALL CPU_TIME ( time_begin )
    time_loop: do

        count=count+1
        t = count*tstep
!! BRAD 2023-02-02:
        if(tf>0.and.t>tf) exit time_loop  !only do the big "do" loop while t<tf

!! BRAD 2023-01-10:
!            if(mod(count,100000)==0) then
!            end if

        !at the beginning of each time step, check to see if any of the fibers should be degraded at this time.
        !Degrade fiber before moving and binding/unbinding tPA:
        do i=enoFB+1,num
            if(t_degrade(i).lt.t.and.t_degrade(i).gt.0.and.degrade(i).eq.0) then
                !i.e. if degradation time is smaller than t AND bigger than 0 AND the edge hasn't already been degraded
                degrade(i)=-t
!! BRAD 2023-01-13:
                degraded_fibers = degraded_fibers + 1
                last_degrade_time = t
                if (degraded_fibers==num-enoFB) all_fibers_degraded = .True.

                !if there were any tPAs bound to this edge, have them unbind, and reset their leaving times to 0:
                do j=1,M
!! BRAD 2023-01-31: There was a small chance that a molecule could arrive at a fiber just as it was degrading
!!                  Even though it was not bound, it would still be classified for "restricted movement"
!!                  Fixed 'passerby molecule' bug
                    if(V(1,j)==i.and.V(2,j)==1) then
                        V(2,j)=0
                        t_leave(j)=0
                        !also find the new binding time for this molecule !FOLLOWING LINE ADDED 4/21/2011:
                        bind(j)=0  !because the molecule can never rebind to a degraded edge
                    end if
                end do

            end if
        end do


        !!!!!!!!!!!!!!!! Now do the binding/unbinding and moving part of the algorithm !!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do j=1,M
            if(V(2,j)==1) then   !if the molecule is bound
                if(t_leave(j).le.t) then
                    !if the time to tPA leaving is smaller than current time, have the molecule unbind
                    V(2,j)=0

                    !also find the new binding time for this molecule !FOLLOWING ADDED 4/18/2011:
                    r1=urcw1()
                    bind(j)=t-log(r1)/(kon*bs)-tstep/2   !subtract half a time step so that we round to nearest timestep
                    !else if time to tPA leaving is bigger than current time, keep the molecule bound and continue.
                    !we don't need to change anything in this case
                end if
            end if !end if(V(2,j)==1 statement. The above unbinds tPA with leaving time < current time

            if(V(2,j)==0) then   !if the molecule is unbound
                !first check if it will move or not. if it does not move (i.e. r.le.1-q), check if it can bind.
                !if it can bind, bind it, if not leave it unbound at the same location. if it does move, check if it can bind.
                !if it can't bind, move it. if it can bind, pick a random number and see if r>(t-bind(j))/tstep.
                !if it is bigger, move it. if it isn't bigger, bind it.
                r=urcw1()
                z=V(1,j)

                if (r.le.(1-q)) then   !i.e. if we stay on current edge
                    !check if molecule j can bind
                    if(bind(j).lt.t.and.bind(j).gt.0.and.degrade(V(1,j)).eq.0) then
                        !i.e. if binding time is smaller than t AND bigger than 0  AND the edge hasn't already been degraded
                        V(2,j)=1  !then the molecule binds
                        bind(j)=0 !reset the binding time to 0
                        r3=urcw1()
                        countbind=countbind+1

!! BRAD 2023-01-04:
                        total_binds = total_binds + 1

                        !put a 1 in entry equal to edge number. each second I'll sum up the number of 1 entries, which
                        !will tell me the number of independent bindings (not necessarily successful ones) at each second
                        bind1(V(1,j))=1

                        !find the time that tPA will unbind:
                        colr2=0
                        do i=1,101
                            if(CDFtPA(i).ge.r3)then
                                colr2 = i    !find the first place on CDF that's bigger than r3
                                exit
                            end if
                        enddo

                        if(colr2==0.or.colr2==1) write(*,*)'PROBLEM: colr2 should not equal 0 since CDFtPA goes between 0 and 1 exactly', colr2

                        percent2 = (CDFtPA(colr2)-r3)/(CDFtPA(colr2)-CDFtPA(colr2-1))
                        ttPA = (tsec1(colr2)-(tsec1(colr2)-tsec1(colr2-1))*percent2)
                        !end if

                        t_leave(j) = t + ttPA - tstep/2 !time tPA leaves is current time plus leaving time drawn from distribution
                                                        !minus half a time step so we round to nearest timestep


                        r4=urcw1()
                        r400=ceiling(r4*nummicro)+1


                        !if r400=501 or 101, redefine it to be the 500th or 100th bin, since we don't have 501/101 entries
                        if(r400==nummicro+1) then
                            r400=nummicro
                        end if !for if(r400==nummicro+1) loop


                        if(r400.le.lenlysismat(colr2-1)) then !only have lysis if the random number puts us in a bin that's < the first place we have a 6000, i.e. the first place lysis doesn't happen
                            if(r400==lenlysismat(colr2-1)) then
                                !percent4 = (CDFlys(r400)-r4)/CDFlys(r400)
                                rmicro = lysismat(r400-1,colr2-1)!tseclys(r400)-tseclys(r400)*percent4
                            else
                                percent4 = r400-1-r4*nummicro
                                rmicro = (lysismat(r400,colr2-1)-(lysismat(r400,colr2-1)-lysismat(r400-1,colr2-1))*percent4)
                            end if

                            if(t_degrade(V(1,j))==0) then              !if no tPA has landed on this edge before
                                t_degrade(V(1,j)) = t + rmicro - tstep/2 !time at which degradation occurs is current time plus
                                                                        !the cutting time obtained from the lysis time function
                                                                        !minus half a time step so we round to nearest time step
                                countindep=countindep+1 !save the number of independent binding events
                            else                        !if tPA has previously landed on this edge and dictated a degradation time,
                                t_degrade(V(1,j)) = min(t_degrade(V(1,j)),(t+rmicro-tstep/2)) !choose the smallest time, because
                                                                                                !that's what will happen first
                            end if

                        end if !for if(r400.le.lenlysismat) loop



                        !if the molecule canNOT bind this timestep, it stays unbound on the current edge so we don't have to
                        !change anything. It remains V(1,j)=V(1,j), V(2,j)=0.

                    end if !(end bind(j)...statement)'


                else                     !if r>1-q, i.e. if molecule j has the possibility to move
                    !check if molecule j can bind
                    if(bind(j).lt.t.and.bind(j).gt.0.and.degrade(V(1,j)).eq.0) then
                        !i.e. if binding time is smaller than t AND bigger than 0  AND the edge hasn't already been degraded.
                        !then the molecule could bind. To determine whether it binds or moves, draw a random number, r2.
                        !if r2.gt.(t-bind(j))/tstep, then movement happened before binding, so move the molecule rather than bind it

                        r2=urcw1()
                        if(r2.gt.(t-bind(j))/tstep) then   !if r2 is such that movement happened before binding, move the molecule
                                                           !and calculate the new binding time associated with the new edge
                            if ((1-q).lt.r.and.r.le.((1-q)+q/8)) then
                                V(1,j) = neighborc(1,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round

                            elseif (((1-q)+q/8).lt.r.and.r.le.((1-q)+2*q/8)) then
                                V(1,j) = neighborc(2,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round

                            elseif (((1-q)+2*q/8).lt.r.and.r.le.((1-q)+3*q/8)) then
                                V(1,j) = neighborc(3,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round

                            elseif (((1-q)+3*q/8).lt.r.and.r.le.((1-q)+4*q/8)) then
                                V(1,j) = neighborc(4,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round

                            elseif (((1-q)+4*q/8).lt.r.and.r.le.((1-q)+5*q/8)) then
                                V(1,j) = neighborc(5,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round

                            elseif (((1-q)+5*q/8).lt.r.and.r.le.((1-q)+6*q/8)) then
                                V(1,j) = neighborc(6,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round

                            elseif (((1-q)+6*q/8).lt.r.and.r.le.((1-q)+7*q/8)) then
                                V(1,j) = neighborc(7,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round

                            elseif (((1-q)+7*q/8).lt.r.and.r.le.((1-q)+8*q/8)) then
                                V(1,j) = neighborc(8,z)
                                r1=urcw1()
                                bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round
                            end if !(for diffusion part)

!! BRAD 2023-01-13:
                            total_regular_moves = total_regular_moves + 1

                        else     !i.e. if r2 is less than or equal to (t-bind(j))/tstep, have the molecule bind
                            V(2,j)=1  !then the molecule binds
                            bind(j)=0 !reset the binding time to 0
                            r3=urcw1()

!! BRAD 2023-01-04:
                            total_binds = total_binds + 1

                            !put a 1 in entry equal to edge number. each second I'll sum up the number of 1 entries, which
                            !will tell me the number of independent bindings (not necessarily successful ones) at each second
                            bind1(V(1,j))=1

                            !find the time that tPA will unbind:
                            colr2=0
                            do i=1,101
                                if(CDFtPA(i).ge.r3)then
                                    colr2 = i    !find the first place on CDF that's bigger than r3
                                    exit
                                end if
                            end do

                            if(colr2==0.or.colr2==1) write(*,*)'PROBLEM: colr2 should not equal 0 since CDFtPA goes between 0 and 1 exactly', colr2

                            percent2 = (CDFtPA(colr2)-r3)/(CDFtPA(colr2)-CDFtPA(colr2-1))
                            ttPA = (tsec1(colr2)-(tsec1(colr2)-tsec1(colr2-1))*percent2)
                            !end if

                            t_leave(j) = t + ttPA - tstep/2 !time tPA leaves is current time plus leaving time drawn from distribution
                                                    !minus half a time step so we round to nearest timestep


                            !Using the tPA leaving time, find the lysis time by using the lysis time distribution for the given ttPA

                            r4=urcw1()
                            r400=ceiling(r4*nummicro)+1

                            !if r400=501, redefine it to be the 500th bin, since we don't have 501 entries
                            if(r400==nummicro+1) then
                                r400=nummicro
                            end if !for if(r400==nummicro+1) loop

                            if(r400.le.lenlysismat(colr2-1)) then !only have lysis if the random number puts us in a bin that's < the first place we have a 6000, i.e. the first place lysis doesn't happen
                                if(r400==lenlysismat(colr2-1)) then
                                    !percent4 = (CDFlys(r400)-r4)/CDFlys(r400)
                                    rmicro = lysismat(r400-1,colr2-1)!tseclys(r400)-tseclys(r400)*percent4
                                else
                                    percent4 = r400-1-r4*nummicro
                                    rmicro = (lysismat(r400,colr2-1)-(lysismat(r400,colr2-1)-lysismat(r400-1,colr2-1))*percent4)
                                end if

                                if(t_degrade(V(1,j))==0) then              !if no tPA has landed on this edge before
                                    t_degrade(V(1,j)) = t + rmicro - tstep/2 !time at which degradation occurs is current time plus
                                                                             !the cutting time obtained from the lysis time function
                                                                             !minus half a time step so we round to nearest time step
                                    countindep=countindep+1 !save the number of independent binding events
                                else               !if tPA has previously landed on this edge and dictated a degradation time,
                                    t_degrade(V(1,j)) = min(t_degrade(V(1,j)),(t+rmicro-tstep/2)) !choose the smallest time, because
                                                                                                  !that's what will happen first
                                end if

                            end if !for if(r400.le.lenlysismat) loop

                        end if  !end r2 statement


                    else         !if molecule j canNOT bind this timestep, just have it move

                        if ((1-q).lt.r.and.r.le.((1-q)+q/8)) then
                            V(1,j) = neighborc(1,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                             !which tPA will bind, minus half a time step so we round

                        elseif (((1-q)+q/8).lt.r.and.r.le.((1-q)+2*q/8)) then
                            V(1,j) = neighborc(2,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                             !which tPA will bind, minus half a time step so we round

                        elseif (((1-q)+2*q/8).lt.r.and.r.le.((1-q)+3*q/8)) then
                            V(1,j) = neighborc(3,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                             !which tPA will bind, minus half a time step so we round

                        elseif (((1-q)+3*q/8).lt.r.and.r.le.((1-q)+4*q/8)) then
                            V(1,j) = neighborc(4,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                             !which tPA will bind, minus half a time step so we round

                        elseif (((1-q)+4*q/8).lt.r.and.r.le.((1-q)+5*q/8)) then
                            V(1,j) = neighborc(5,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                             !which tPA will bind, minus half a time step so we round

                        elseif (((1-q)+5*q/8).lt.r.and.r.le.((1-q)+6*q/8)) then
                            V(1,j) = neighborc(6,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                             !which tPA will bind, minus half a time step so we round

                        elseif (((1-q)+6*q/8).lt.r.and.r.le.((1-q)+7*q/8)) then
                            V(1,j) = neighborc(7,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                             !which tPA will bind, minus half a time step so we round

                        elseif (((1-q)+7*q/8).lt.r.and.r.le.((1-q)+8*q/8)) then
                            V(1,j) = neighborc(8,z)
                            r1=urcw1()
                            bind(j)=t-log(r1)/(kon*bs)-tstep/2 !random time chosen from exponential distribution at
                                                                 !which tPA will bind, minus half a time step so we round
                        end if !(for diffusion part)

!! BRAD 2023-01-04:
                        total_regular_moves = total_regular_moves + 1
                    end if !end bind(j) statement
                end if !end r statement

            end if !end V(2,j)=0 statement

            !for the first run only, at the end of each step, check to see if the molecule hit a fiber in the back row
            !      if(istat==1) then
            if (V(1,j)>=backrow.and.yesfpt(j)==0) then !if the molecule is on the back row and it hasn't made it there before
!! BRAD 2023-01-13:
                    reached_back_row = reached_back_row + 1
                    if (reached_back_row >= 0.95*M) most_molecules_passed = .True.

                mfpt(j)=t
                yesfpt(j)=1 !set entry to 1 so we don't track this molecule any more
            end if !end if(V(1,j)>=backrow....) loop
            !      end if !end if(istat=1) loop
        end do !for j=1,M loop


        Ninteger=int(t)


        !if(Ninteger>Nsave) then !if the current time is the 1st past a new a whole second, e.g. t=3.001, save degrade and V
        !write(degunit)    degrade(1:num)
        !write(tunit)  t
        !Nsave=Nsave+1
        !Vedgenext(Nsave+1,:) = V(1,:)
        !Vboundnext(Nsave+1,:) = V(2,:)
        !degnext(Nsave+1,:) = degrade(1:num)
        !tsave(Nsave+1) = t
        !countbindV(istat,Nsave)=countbind
        !countindepV(istat,Nsave)=countindep
        !bind1V(istat,Nsave)=sum(bind1)
        !end if

        if(Ninteger>=Nsave) then !if the current time is the 1st past a new a 10 seconds, e.g. t=10.001, save degrade and V
            write(degunit)    degrade(1:num)
            write(tunit)  t
!! BRAD 2023-01-21:
            write(t_degrade_unit) t_degrade(:)
            write(m_location_unit) V(1,:)
            write(m_bound_unit) V(2,:)

            Nsave=Nsave+10
            cNsave=cNsave+1
            Vedgenext(cNsave+1,:) = V(1,:)
            Vboundnext(cNsave+1,:) = V(2,:)
            degnext(cNsave+1,:) = degrade(1:num)
            tsave(cNsave+1) = t
            !!!!!COMMENTED OUT BELOW ON 5/16/16 BECAUSE I DON'T USE THIS DATA IN ANY POST-PROCESSING
            !countbindV(istat,cNsave)=countbind
            !countindepV(istat,cNsave)=countindep
            !bind1V(istat,cNsave)=sum(bind1)

!! BRAD 2023-02-02:
            rounded_time = real(t)
            degraded_percent = real(degraded_fibers)/(num-enoFB)*100
            reached_back_row_percent = real(reached_back_row)/M*100
            write(*,'(A,F7.2,A,I5,A,F5.1,A,I5,A,F5.1,A)')'After ',rounded_time,' sec, ',&
            degraded_fibers,' fibers are degraded (',degraded_percent,'% of total) and ',&
            reached_back_row,' molecules have reached the back row (',&
            reached_back_row_percent,'% of total).'

            if (all_fibers_degraded.and.most_molecules_passed) exit time_loop
        end if


    end do time_loop !for time loop

!! BRAD 2023-01-06:
    CALL CPU_TIME ( time_end )

    write(*,*)'Processing time: ', time_end - time_begin, ' sec'
    write(*,*)'Total Binds: ',total_binds
    write(*,*)'Total Regular Moves: ',total_regular_moves
    write(*,*)'Total Restricted Moves: ',total_restricted_moves
    write(*,*)'Molecules that reached back row: ',reached_back_row
    write(*,*)'Last fiber degraded at: ',last_degrade_time,' sec'

!! BRAD 2023-02-02:
    write(mfptunit) mfpt(:)

    Nsavevect(istat)=cNsave !CHANGED TO CNSAVE FROM NSAVE 12/17/14
    front=0
    degold=0

    !NOW PROCESS THE DATA WE OBTAINED FROM THE ABOVE RUN
    nplt=Nsavevect(istat)+2 !+1 because I saved once at the beginning, and another +1 to account for the final time point
    write(*,*)'nplt=',nplt
    write(*,*)'r4=',r4

    !  do i=2,nplt
    !
    !      degold=degnext(i,:)
    !      ind=0
    !      place=0
    !   do j=1,N
    !       do k=1,F-1
    !           ind(k) = (3*N-1)*(k-1) + 2*N + j-1 !ind is a vector containing the vertical planar edge numbers above node j
    !           place(k) = degold(ind(k))  !place(k) is the degradation state of each edge above node j
    !       enddo
    !       call findfirstreal(place,F-1,0,zero1)  !find the first undegraded vertical edge above node j
    !       front(i-1,j) = zero1
    !   enddo
    !  enddo
    !
    !
    !  !the columns of front correspond to the nodes on the x-axis, and the
    !  !rows are successive time steps of varying length. The entries are the
    !  !y-position of the first undegraded edge at each x location
    !
    !  !use this information to calculate the front speed. How? First calculate
    !  !the speed at each x-location (so I'll have N speeds. Then average them or
    !  !something)
    !
    !  !all x-locations start with y=1, so find the first place they deviate from
    !  !there:
    !
    !  firstdeg=0
    !  deglast=0
    !
    !  do i=1,N
    !      call findfirstineq(front(:,i),tf,1,fdeg)
    !      if(fdeg==0) then
    !          firstdeg(i)=1
    !      else
    !          firstdeg(i)=fdeg
    !      end if
    !  enddo
    !
    !  do i=1,N
    !      call findfirst(front(2:tf,i),tf-1,0,first0)
    !      deglast(i)=first0
    !  enddo
    !
    !  !so firstdeg saves the row # (i.e. time) at which the front first moves and
    !  !deglast is the time at which total degradation in a single row occurred.
    !  !if deglast=0, then total degradation did not occur
    !
    !  move=0
    !  move(1,:)=firstdeg
    !
    !  do j=2,N
    !      do i=1,N
    !          if(move(j-1,i)==0) then
    !              temp=0
    !          else
    !              call findfirstineq(front(:,i),tf,front(move(j-1,i),i),temp)
    !          end if
    !          move(j,i)=temp
    !      enddo
    !  enddo
    !
    !  !so now "move" saves the saved-time-step at which the front moves for each x location
    !
    !  do i=1,N
    !      call findintineq(move(:,i),N,0,lasti)
    !      lastmove(i,istat)=lasti
    !  enddo
    !
    !  plotstuff=0
    !  plotstuff2=0
    !
    !  do j=1,N
    !      do i=1,lastmove(j,istat)
    !          plotstuff(j,i)=front(move(i,j),j)
    !          plotstuff2(j,i)=(plotstuff(j,i)-1)*dist
    !      enddo
    !  enddo
    !
    !  !now plotstuff has in each row the successive y-positions of x-location
    !  !corresponding to row number, and plotstuff2 has in each row the successive
    !  !y-positions (in microns, instead of node #) of x-location corresponding to
    !  !row number
    !
    !  !In order to plot this and finish the calculations in Matlab, I need to save plotstuff2, lastmove, and move
    !
    !  write(moveunit) move(:,:)
    !  write(plotunit) plotstuff2(:,:)
    !
    !  if(istat==1)then  !choose how many runs you want to save to make a movie
    !     !!!!!COMMENTED OUT BELOW ON 5/16/16 BECAUSE I DON'T USE THIS DATA IN ANY POST-PROCESSING
    !     !write(degnextunit) degnext(:,:)
    !     !write(Venextunit) Vedgenext(:,:)
    !     !write(Vbdnextunit) Vboundnext(:,:)
    !
    !     !!!DO MORE MOVIE PROCESSING BEFORE GOING TO MATLAB
    !
    !      X1plot=0
    !      Y1plot=0
    !      X2plot=0
    !      Y2plot=0
    !      Xvplot=0
    !      Yvplot=0
    !      freetPA=0
    !      bdtPA=0
    !
    !     countintact2=0
    !     counth=0
    !     countv=0
    !     countpv=0
    !     do i=enoFB,num
    !        if(degnext(1,i)==0) then
    !          countintact2=countintact2+1
    !          intact2(countintact2)=i   !finds the undegraded edge numbers
    !        end if
    !     enddo
    !
    !     lenintact2=countintact2
    !
    !     !Assume you have a grid with N nodes in each row of the lattice, and you
    !     !assign the nodes a number simply by counting nodes, starting at the bottom
    !     !left and moving right, then going to the left of row two and moving right,
    !     !etc. in column K, endpts has the node numbers corresponding to the
    !     !endpts of fiber (i.e. edge) K. K=1,...,num
    !
    !     !for vertical edges:
    !     do i=1,F
    !         do j=1,N
    !             endpts(1,(3*N-1)*(i-1)+1+2*(j-1)) = (i-1)*N + j
    !             endpts(2,(3*N-1)*(i-1)+1+2*(j-1)) = (i-1)*N + j
    !         end do
    !     end do
    !
    !     !for horizontal edges
    !     do i=1,F
    !         do j=1,N-1
    !             endpts(1,(3*N-1)*(i-1)+2+2*(j-1)) = j + N*(i-1)
    !             endpts(2,(3*N-1)*(i-1)+2+2*(j-1)) = j + N*(i-1) + 1
    !         end do
    !     end do
    !
    !     !for planar veritcal edges
    !     do i=1,F-1
    !         do j=1,N
    !             endpts(1,(3*N-1)*(i-1)+2*N+(j-1)) = j + N*(i-1)
    !             endpts(2,(3*N-1)*(i-1)+2*N+(j-1)) = j + N*(i-1) + N
    !         end do
    !     end do
    !
    !
    !     do jj=1,lenintact2
    !         !horizontal edges
    !         do iplt=Ffree,F
    !           do j=1,N-1
    !             if(intact2(jj)==2*(j-1)+2+(3*N-1)*(iplt-1)) then    !if we have a horizontal edge
    !                yplace=endpts(1,intact2(jj))/N+1  !find the y value at which the horizontal edge occurs
    !                x2=nint((real(endpts(2,intact2(jj)))/real(N)-floor(real(endpts(2,intact2(jj)))/real(N)))*N)
    !                if(x2==0) x2=N    !if it says the RHS endpoint is 0, force it to actually be N (because otherwise
    !                                  !is says we should plot from N-1 to 0)
    !                counth=counth+1
    !                X1plot(1,counth)=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)
    !                X1plot(2,counth)=x2
    !                Y1plot(1,counth)=yplace
    !                Y1plot(2,counth)=yplace
    !             end if
    !           enddo
    !        enddo
    !        !vertical (planar) edges
    !        do j=Ffree,F-1
    !          do k=1,N
    !             if(intact2(jj)==(3*N-1)*(j-1)+2*N+(k-1)) then   !if we have a vertical (planar) edge
    !               xplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)
    !                                                                !find the x value at which the vertical edge occurs
    !               y1=endpts(1,intact2(jj))/N+1  !find the bottom endpoint of the vertical edge
    !               y2=endpts(2,intact2(jj))/N+1
    !               if(xplace==0) then
    !                 xplace=N
    !                 y1=endpts(1,intact2(jj))/N
    !                 y2=endpts(2,intact2(jj))/N
    !               end if
    !               countpv=countpv+1
    !               X2plot(1,countpv)=xplace
    !               X2plot(2,countpv)=xplace
    !               Y2plot(1,countpv)=y1
    !               Y2plot(2,countpv)=y2
    !            end if
    !         end do
    !       end do
    !       !vertical edges
    !       do i=Ffree,F
    !         do j=1,N
    !            if(intact2(jj)==(3*N-1)*(i-1)+1+2*(j-1)) then  !if we have a vertical edge
    !              yvplace=(endpts(1,intact2(jj))/N)+1  !find the y value at which the vertical edge occurs
    !              xvplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)
    !                                                                      !find the x value at which the vertical edge occurs
    !
    !              if(xvplace==0) then
    !                xvplace=N
    !                yvplace=endpts(1,intact2(jj))/N
    !              end if
    !              countv=countv+1
    !              Xvplot(countv)=xvplace
    !              Yvplot(countv)=yvplace
    !            end if
    !         end do
    !      end do
    !    end do  !for jj loop
    !
    !            open(x1unit,file=ADJUSTL('data/' // expCode // '/X1plot' // outFileCode),form=filetype)
    !            open(x2unit,file=ADJUSTL('data/' // expCode // '/X2plot' // outFileCode),form=filetype)
    !            open(y1unit,file=ADJUSTL('data/' // expCode // '/Y1plot' // outFileCode),form=filetype)
    !            open(y2unit,file=ADJUSTL('data/' // expCode // '/Y2plot' // outFileCode),form=filetype)
    !            open(xvunit,file=ADJUSTL('data/' // expCode // '/Xvplot' // outFileCode),form=filetype)
    !            open(yvunit,file=ADJUSTL('data/' // expCode // '/Yvplot' // outFileCode),form=filetype)
    !
    !write(x1unit) X1plot
    !write(x2unit) X2plot
    !write(y1unit) Y1plot
    !write(y2unit) Y2plot
    !write(xvunit) Xvplot
    !write(yvunit) Yvplot
    !
    !!!Now do location and boundedness of tPA
    !     do i=1,F
    !        do j=1,N-1
    !           !horizontal edges
    !           do jplt=1,M
    !              if(Vedgenext(1,jplt)==2*(j-1)+2+(3*N-1)*(i-1)) then
    !                  Vedgeplace=(endpts(1,Vedgenext(1,jplt))/N)+1
    !                  Vx=nint((real(endpts(2,Vedgenext(1,jplt)))/real(N)-floor(real(endpts(2,Vedgenext(1,jplt)))/real(N)))*N)
    !                  if(Vx==0) then
    !                     Vx=N
    !                  end if
    !                  if(Vboundnext(1,jplt)==0) then   !if unbound, will plot in black
    !                     bdtPA(1,jplt)=-1
    !                     bdtPA(2,jplt)=-1
    !                     freetPA(1,jplt)=Vx-0.5d+00
    !                     freetPA(2,jplt)=Vedgeplace
    !                  elseif(Vboundnext(1,jplt)==1) then  !if bound, will plot in green
    !                     bdtPA(1,jplt)=Vx-0.5d+00
    !                     bdtPA(2,jplt)=Vedgeplace
    !                     freetPA(1,jplt)=-1
    !                     freetPA(2,jplt)=-1
    !                  end if
    !              end if
    !           end do
    !       end do
    !     end do
    !     do j=1,F-1
    !        do k=1,N
    !          !vertical (planar) edges
    !           do kplt=1,M
    !              if(Vedgenext(1,kplt)==(3*N-1)*(j-1)+2*N+(k-1)) then
    !                 xVedgeplace=nint((real(endpts(1,Vedgenext(1,kplt)))/real(N)-floor(real(endpts(1,Vedgenext(1,kplt)))/real(N)))*N)  !find the x value at which the vertical edge occurs
    !                 Vy1=endpts(1,Vedgenext(1,kplt))/N+1  !find the bottom endpoint of the vertical edge
    !                 if(xVedgeplace==0) then
    !                    xVedgeplace=N
    !                    Vy1=endpts(1,Vedgenext(1,kplt))/N
    !                 end if
    !                 if(Vboundnext(1,kplt)==0) then   !if unbound, plot in black
    !                    bdtPA(1,kplt)=-1
    !                    bdtPA(2,kplt)=-1
    !                    freetPA(1,kplt)=xVedgeplace
    !                    freetPA(2,kplt)=Vy1+0.5d+00
    !                 elseif(Vboundnext(1,kplt)==1) then  !if bound, plot in green
    !                    bdtPA(1,kplt)=xVedgeplace
    !                    bdtPA(2,kplt)=Vy1+0.5d+00
    !                    freetPA(1,kplt)=-1
    !                    freetPA(2,kplt)=-1
    !                 end if
    !               end if
    !            end do
    !        end do
    !     end do
    !     do i=1,F
    !        do j=1,N
    !         !vertical edges
    !          do kjplt=1,M
    !             if(Vedgenext(1,kjplt)==(3*N-1)*(i-1)+1+2*(j-1)) then
    !                vertplace=1+(j-1)   !find the x value at which the vertical edge occurs
    !                Vyvert=endpts(1,Vedgenext(1,kjplt))/N+1  !find the bottom endpoint of the vertical edge
    !                if(vertplace==N) then
    !                  !vertplace=N;
    !                   Vyvert=endpts(1,Vedgenext(1,kjplt))/N
    !                end if
    !                if(Vboundnext(1,kjplt)==0) then  !if unbound, plot in black
    !                   bdtPA(1,kjplt)=-1
    !                   bdtPA(2,kjplt)=-1
    !                   freetPA(1,kjplt)=vertplace
    !                   freetPA(2,kjplt)=Vyvert
    !                elseif(Vboundnext(1,kjplt)==1) then   !if bound, plot in green
    !                   bdtPA(1,kjplt)=vertplace
    !                   bdtPA(2,kjplt)=Vyvert
    !                   freetPA(1,kjplt)=-1
    !                   freetPA(2,kjplt)=-1
    !                end if
    !              end if
    !           end do
    !        end do
    !     end do
    !
    !
    !            open(tPAbdunit,file=ADJUSTL('data/' // expCode // '/tPAbd' // outFileCode),form=filetype)
    !            open(tPAfreeunit,file=ADJUSTL('data/' // expCode // '/tPAfree' // outFileCode),form=filetype)
    !
    !write(tPAbdunit) bdtPA
    !write(tPAfreeunit) freetPA
    !
    !
    !!now save different timestep so I can make a matlab movie
    !   do imod=2,nplt
    !      if(mod(imod,6)==0) then  !this means we plot approximately every minute, since we saved data every 10 seconds
    !
    !         X1plot=0
    !         Y1plot=0
    !         X2plot=0
    !         Y2plot=0
    !         Xvplot=0
    !         Yvplot=0
    !         freetPA=0
    !         bdtPA=0
    !
    !         intact2=0
    !         countintact2=0
    !         counth=0
    !         countv=0
    !         countpv=0
    !         do i=enoFB,num  !could do this as"do i=enoFB,num" since we know the first enoFB edges won't equal 0
    !           if(degnext(imod,i)==0) then
    !             countintact2=countintact2+1
    !             intact2(countintact2)=i   !finds the undegraded edge numbers
    !           end if
    !         enddo
    !
    !         lenintact2=countintact2
    !
    !      do jj=1,lenintact2
    !         !horizontal edges
    !         do iplt=Ffree,F
    !           do j=1,N-1
    !             if(intact2(jj)==2*(j-1)+2+(3*N-1)*(iplt-1)) then    !if we have a horizontal edge
    !                yplace=endpts(1,intact2(jj))/N+1  !find the y value at which the horizontal edge occurs
    !                x2=nint((real(endpts(2,intact2(jj)))/real(N)-floor(real(endpts(2,intact2(jj)))/real(N)))*N)
    !                         if(x2==0) x2=N    !if it says the RHS endpoint is 0, force it to actually be N (because otherwise
    !                                           !is says we should plot from N-1 to 0)
    !                         counth=counth+1
    !                         X1plot(1,counth)=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)
    !                         X1plot(2,counth)=x2
    !                         Y1plot(1,counth)=yplace
    !                         Y1plot(2,counth)=yplace
    !             end if
    !           enddo
    !        enddo
    !        !vertical (planar) edges
    !        do j=Ffree,F-1
    !          do k=1,N
    !             if(intact2(jj)==(3*N-1)*(j-1)+2*N+(k-1)) then   !if we have a vertical (planar) edge
    !               xplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)
    !                                                                !find the x value at which the vertical edge occurs
    !               y1=endpts(1,intact2(jj))/N+1  !find the bottom endpoint of the vertical edge
    !               y2=endpts(2,intact2(jj))/N+1
    !               if(xplace==0) then
    !                 xplace=N
    !                 y1=endpts(1,intact2(jj))/N
    !                 y2=endpts(2,intact2(jj))/N
    !               end if
    !               countpv=countpv+1
    !               X2plot(1,countpv)=xplace
    !               X2plot(2,countpv)=xplace
    !               Y2plot(1,countpv)=y1
    !               Y2plot(2,countpv)=y2
    !            end if
    !         end do
    !       end do
    !       !vertical edges
    !       do i=Ffree,F
    !         do j=1,N
    !            if(intact2(jj)==(3*N-1)*(i-1)+1+2*(j-1)) then  !if we have a vertical edge
    !              yvplace=(endpts(1,intact2(jj))/N)+1  !find the y value at which the vertical edge occurs
    !              xvplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)
    !                                                                      !find the x value at which the vertical edge occurs
    !              if(xvplace==0) then
    !                xvplace=N
    !                yvplace=endpts(1,intact2(jj))/N
    !              end if
    !              countv=countv+1
    !              Xvplot(countv)=xvplace
    !              Yvplot(countv)=yvplace
    !            end if
    !         end do
    !      end do
    !    end do  !for jj loop
    !
    !  write(x1unit) X1plot
    !  write(x2unit) X2plot
    !  write(y1unit) Y1plot
    !  write(y2unit) Y2plot
    !  write(xvunit) Xvplot
    !  write(yvunit) Yvplot
    !
    !!!Now do location and boundedness of tPA
    !     do i=1,F
    !        do j=1,N-1
    !           !horizontal edges
    !           do jplt=1,M
    !              if(Vedgenext(imod,jplt)==2*(j-1)+2+(3*N-1)*(i-1)) then
    !                  Vedgeplace=(endpts(1,Vedgenext(imod,jplt))/N)+1
    !                  Vx=nint((real(endpts(2,Vedgenext(imod,jplt)))/real(N)-floor(real(endpts(2,Vedgenext(imod,jplt)))/real(N)))*N)
    !                  if(Vx==0) then
    !                     Vx=N
    !                  end if
    !                  if(Vboundnext(imod,jplt)==0) then   !if unbound, will plot in black
    !                     bdtPA(1,jplt)=-1
    !                     bdtPA(2,jplt)=-1
    !                     freetPA(1,jplt)=Vx-0.5d+00
    !                     freetPA(2,jplt)=Vedgeplace
    !                  elseif(Vboundnext(imod,jplt)==1) then  !if bound, will plot in green
    !                     bdtPA(1,jplt)=Vx-0.5d+00
    !                     bdtPA(2,jplt)=Vedgeplace
    !                     freetPA(1,jplt)=-1
    !                     freetPA(2,jplt)=-1
    !                  end if
    !              end if
    !           end do
    !       end do
    !     end do
    !     do j=1,F-1
    !        do k=1,N
    !          !vertical (planar) edges
    !           do kplt=1,M
    !              if(Vedgenext(imod,kplt)==(3*N-1)*(j-1)+2*N+(k-1)) then
    !                 xVedgeplace=nint((real(endpts(1,Vedgenext(imod,kplt)))/real(N)-floor(real(endpts(1,Vedgenext(imod,kplt)))/real(N)))*N)  !find the x value at which the vertical edge occurs
    !                 Vy1=endpts(1,Vedgenext(imod,kplt))/N+1  !find the bottom endpoint of the vertical edge
    !                 if(xVedgeplace==0) then
    !                    xVedgeplace=N
    !                    Vy1=endpts(1,Vedgenext(imod,kplt))/N
    !                 end if
    !                 if(Vboundnext(imod,kplt)==0) then   !if unbound, plot in black
    !                    bdtPA(1,kplt)=-1
    !                    bdtPA(2,kplt)=-1
    !                    freetPA(1,kplt)=xVedgeplace
    !                    freetPA(2,kplt)=Vy1+0.5d+00
    !                 elseif(Vboundnext(imod,kplt)==1) then  !if bound, plot in green
    !                    bdtPA(1,kplt)=xVedgeplace
    !                    bdtPA(2,kplt)=Vy1+0.5d+00
    !                    freetPA(1,kplt)=-1
    !                    freetPA(2,kplt)=-1
    !                 end if
    !               end if
    !            end do
    !        end do
    !     end do
    !     do i=1,F
    !        do j=1,N
    !         !vertical edges
    !          do kjplt=1,M
    !             if(Vedgenext(imod,kjplt)==(3*N-1)*(i-1)+1+2*(j-1)) then
    !                vertplace=1+(j-1)   !find the x value at which the vertical edge occurs
    !                Vyvert=endpts(1,Vedgenext(imod,kjplt))/N+1  !find the bottom endpoint of the vertical edge
    !                if(vertplace==N) then
    !                  !vertplace=N;
    !                   Vyvert=endpts(1,Vedgenext(imod,kjplt))/N
    !                end if
    !                if(Vboundnext(imod,kjplt)==0) then  !if unbound, plot in black
    !                   bdtPA(1,kjplt)=-1
    !                   bdtPA(2,kjplt)=-1
    !                   freetPA(1,kjplt)=vertplace
    !                   freetPA(2,kjplt)=Vyvert
    !                elseif(Vboundnext(imod,kjplt)==1) then   !if bound, plot in green
    !                   bdtPA(1,kjplt)=vertplace
    !                   bdtPA(2,kjplt)=Vyvert
    !                   freetPA(1,kjplt)=-1
    !                   freetPA(2,kjplt)=-1
    !                end if
    !              end if
    !           end do
    !        end do
    !     end do
    !
    !  write(tPAbdunit) bdtPA
    !  write(tPAfreeunit) freetPA
    !
    !   end if !for if mod(imod,60) loop
    !  end do !for imod loop
    !
    !close(x1unit)
    !close(x2unit)
    !close(y1unit)
    !close(y2unit)
    !close(xvunit)
    !close(yvunit)
    !close(tPAbdunit)
    !close(tPAfreeunit)
    !
    !
    !!!!!! END ADDED STUFF FOR MOVIE
    !  end if
    !



    countindepV(istat,tf)=countindep
end do stats_loop !for stats loop

write(*,*)'Nsavevect=',Nsavevect(:)

!!!!!COMMENTED OUT BELOW ON 5/16/16 BECAUSE I DON'T USE THIS DATA IN ANY POST-PROCESSING
!write(cbindunit) countbindV
!write(cindunit) countindepV
!write(bind1unit) bind1V
write(Nunit) Nsavevect(:)
write(lastmoveunit) lastmove(:,:)
!write(mfptunit) mfpt(:)

close(degunit)
close(Nunit)
close(tunit)
close(moveunit)
close(lastmoveunit)
close(plotunit)
close(mfptunit)

!! BRAD 2023-01-21:
close(t_degrade_unit)
close(m_location_unit)
close(m_bound_unit)
!        close(m_bind_time_unit)

!!!!!COMMENTED OUT BELOW ON 5/16/16 BECAUSE I DON'T USE THIS DATA IN ANY POST-PROCESSING
!close(degnextunit)
!close(Venextunit)
!close(Vbdnextunit)
!close(cbindunit)
!close(cindunit)
!close(bind1unit)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE CONCAT - concatenates a vector and a single number
!
subroutine concat(g,size2,h,sizeg,gnew)

        integer  :: sizeg,size2
        integer  :: h
        integer, dimension(size2)   :: g,gnew

gnew(1:sizeg)=g(1:sizeg)
gnew(sizeg+1)=h

end subroutine concat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE VURCW1 - creates a vector of random numbers
!
subroutine vurcw1(v, n)
      double precision v(*)
      integer k, n

      do k = 1, n
          v(k) = urcw1()
      end do

      end subroutine vurcw1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE FINDFIRSTREAL - finds the first entry equal to a specified number in a double precision array
!
subroutine findfirstreal(g,sizeg,numspec,loca)
        integer  :: sizeg
        integer  :: i,numspec,loca
        double precision, dimension(sizeg)   :: g

        loca=0
        do i=1,sizeg
              if(g(i)==numspec) then
                  loca=i
                  exit
              end if

        enddo
 

end subroutine findfirstreal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE FINDFIRST - finds the first entry equal to a specified number
!
subroutine findfirst(g,sizeg,numspec,loca)
        integer  :: sizeg
        integer  :: i,numspec,loca
        integer, dimension(sizeg)   :: g

        loca=0
        do i=1,sizeg
              if(g(i)==numspec) then
                  loca=i
                  exit
              end if

        enddo
 

end subroutine findfirst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE FIND10 - finds the first multiple of 10 greater than a given number
!
subroutine find10(numint,ans)
        double precision  :: numint
        integer  :: i,ans,intnumint

        ans=0
        intnumint=ceiling(numint/10)  
        do i=1,intnumint
              if(numint<i*10) then
                  ans=i*10
                  exit
              end if

        enddo
 
        if(ans==0) write(*,*) 'Problem: did not find correct multiple of 10'
 

end subroutine find10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE FINDFIRSTINEQ - finds the first entry greater than a specified number
!
subroutine findfirstineq(g,sizeg,numspec,loca)
        integer  :: sizeg
        integer  :: i,numspec,loca
        integer, dimension(sizeg)   :: g

        loca=0
        do i=1,sizeg
              if(g(i)>numspec) then
                  loca=i
                  exit
              end if

        enddo
 

end subroutine findfirstineq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE FINDINTINEQ - finds number of entries in vector g greater than intspec
!
subroutine findintineq(g,sizeg,intspec,sizegnew)
        integer  :: sizeg,sizegnew
        integer  :: i,countg,intspec
        integer, dimension(sizeg)   :: g


        countg=0

        do i=1,sizeg
              if(g(i)>intspec) then
                  countg=countg+1
                 ! gnew(countg) = i
              end if

        enddo

        sizegnew=countg
 

end subroutine findintineq



end program macrolysis

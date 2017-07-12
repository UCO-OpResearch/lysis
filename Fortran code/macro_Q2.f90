program macrolysis
use ncurses      ! load interface to ncurses(3c) C library
!mplicit none
   
!Runs the macroscale model in a clot with 72.7 nm diameter fibers and pore size. 1.0135 uM. FB conc. = 8.8 uM. THIS CODE ACCOUNTS FOR MICRO RUNS IN WHICH 50,000 OR 10,000 INDEPENDENT SIMULATIONS WERE DONE. CHANGE LINE 16 (nummicro=) to 500 or 100 depending on if 50,000 or 10,000 micro runs were completed
implicit none

integer,parameter  :: N=30 !93  !# of lattice nodes in one row in the horizontal direction
integer,parameter  :: F=11 !121 !71 !81  !# of lattice nodes in one column in the vertical direction
integer,parameter  :: Ffree=4 !29 !3 !13 !1st node in vertical direction containing fibers. so if Ffree=10, then rows 1-9
                               !have no fibers, there's one more row of fiber-free planar veritcal edges, and then
                               !the row starting with the Ffree-th (e.g. 10th) vertical node is a full row of fibers 
integer,parameter  :: stats=1 !10
integer,parameter  :: num=(2*N-1)*F+N*(F-1)
integer,parameter  :: M=500 !43074 !total number of tPA molecules: 43074 is Colin's [tPA]=0.6 nM; 86148 is Colin's [tPA]=1.2 nM; 21588 is Colin's [tPA]=0.3 nM
integer,parameter  :: tf=10*60!15*60 !final time in sec
integer,parameter  :: enoFB=(3*N-1)*(Ffree-1) !the last edge number without fibrin
integer,parameter  :: nummicro=283 !if the number of microscale runs was 50,000, take nummicro=500; if it was 10,000, take nummicro=100
integer  :: i, istat
integer  :: j
integer  :: k
integer  :: ii 
integer  :: Nsave, Ninteger, nplt, cNsave
integer  :: count
integer  :: count2, countc, countcolr4 
integer  :: count66 
integer  :: Ntot2
integer  :: colr2, colr4 
integer  :: z
integer  :: numPLi
double precision     :: t
double precision     :: q
double precision     :: delx
double precision     :: Diff
double precision     :: tstep
double precision     :: num_t
double precision     :: dist
double precision     :: kon
double precision     :: bs
double precision     :: t_bind
double precision     :: percent2, percent4
double precision     :: rmicro, ttPA
double precision     :: r, r1, r3, r2, r4


character(50) :: filetype,formatted
character(70) :: filename1
character(80) :: filename2
character(90) :: filename3
character(95) :: filename4
character(75) :: filename6


integer*1, dimension(num,num)  :: closeneigh
integer, dimension(2,num)    :: endpts
integer, dimension(8,num)    :: neighborc
integer, dimension(2,M)   :: V
double precision, dimension(enoFB)      :: init_state
double precision, dimension(2)         :: p
double precision, dimension(2)         :: pfit
double precision, dimension(101)       :: CDFtPA, CDFlys
double precision, dimension(101)       :: tsec1, tseclys
double precision, dimension(num)      :: degrade
double precision, dimension(num)      :: t_degrade
double precision, dimension(M)        :: t_leave
double precision, dimension(M)        :: bind
integer, dimension(stats)             :: Nsavevect 
integer                              :: name1, name2

logical       :: isBinary =  .True.      ! flag for binary output
integer       :: degunit = 20
integer       :: Vunit = 21
integer       :: V2unit = 22
integer       :: Nunit = 23
integer       :: tunit = 24
character(50) :: degfile       ! degradation vector
character(50) :: Vfile
character(50) :: V2file
character(50) :: Nfile
character(50) :: tfile

!stuff for the random # generator. I need to use the file kiss.o when I compile in order for this to work, and kiss.o
!is obtained by compiling the file kiss.c by doing "cc -c kiss.c".
external :: mscw, kiss32, urcw1
integer :: kiss32, mscw, seed, state(4), old_state(4), ui
double precision :: uf, urcw1

double precision, dimension(M)         :: rvect
double precision, dimension(tf+1,num) :: degnext
integer, dimension(tf+1,M) :: Vedgenext
integer, dimension(tf+1,M) :: Vboundnext
integer, dimension(F-1)  :: ind
double precision, dimension(F-1)  :: place
double precision, dimension(num)  :: degold
double precision, dimension(tf+1) :: tsave
integer  :: zero1
integer, dimension(tf,N)  :: front
integer, dimension(N)  :: firstdeg
integer, dimension(N)  :: deglast
integer, dimension(N,stats)  :: lastmove
integer  :: fdeg
integer  :: first0
integer, dimension(N,N)  :: moved
integer  :: temp
integer  :: lasti
integer, dimension(N,N)  :: plotstuff, totmove,time2plot
double precision, dimension(N,N)  :: plotstuff2
integer       :: moveunit = 25
integer       :: lastmoveunit = 26
integer       :: plotunit = 27
integer       :: degnextunit = 28
integer       :: Venextunit = 29
integer       :: Vbdnextunit = 30
character(50) :: movefile    
character(50) :: lastmovefile 
character(50) :: plotfile   
character(50) :: degnextfile  
character(50) :: Venextfile  
character(50) :: Vbdnextfile  

integer, dimension(num)  :: intact2
integer  :: countintact2, lenintact2, counth, countv, countpv
integer  :: jj, iplt, yplace, x2, xplace, y1, y2, yvplace, xvplace, imod
integer  :: jplt, kplt, kjplt, vertplace, Vyvert, Vy1, xVedgeplace, Vedgeplace, Vx 
integer, dimension(2,F*(N-1))  :: X1plot, Y1plot
integer, dimension(2,N*(F-1))  :: X2plot, Y2plot
integer, dimension(N*F)  :: Xvplot, Yvplot
double precision, dimension(2,M)  :: bdtPA, freetPA
double precision, dimension(nummicro,100) :: lysismat !(100,100) if only did 10,000 micro runs, (500,100) if did 50,000
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
character(50)  :: x1file
character(50)  :: x2file
character(50)  :: y1file
character(50)  :: y2file
character(50)  :: xvfile
character(50)  :: yvfile
character(50)  :: tPAbdfile
character(50)  :: tPAfreefile
character(50)  :: cbindfile
character(50)  :: cindfile
character(50)  :: bind1file
integer :: countbind, countindep
integer, dimension(stats,tf)  :: countbindV, countindepV, bind1V
integer, dimension(num) :: bind1

!Stuff for Brad Paynter
integer :: ierr, ikey
!integer*8 :: stdscreen
integer :: iostat

integer, parameter :: RIGHT = 1
integer, parameter :: LEFT = -1
integer, parameter :: UP = 2
integer, parameter :: DOWN = -2
integer, parameter :: OUT = 3
integer, parameter :: IN = -3

integer, parameter :: gridSize = 1
integer*2, parameter :: black2 = 1
integer*2, parameter :: red2 = 2
integer*2, parameter :: green2 = 3
integer*2, parameter :: yellow2 = 4
integer*2, parameter :: blue2 = 5
integer*2, parameter :: white2 = 6
integer*4, parameter :: black4 = 1
integer*4, parameter :: red4 = 2
integer*4, parameter :: green4 = 3
integer*4, parameter :: yellow4 = 4
integer*4, parameter :: blue4 = 5
integer*4, parameter :: white4 = 6

if( isBinary ) then
	filetype = 'unformatted' !if you compile with gfortran or f95
	!filetype = 'binary'      !if you compile with ifort
else
	filetype = 'formatted'
end if
write(*,*)' filetype=',filetype


write(*,*)' N=',N
write(*,*)' F=',F
write(*,*)' Ffree=',Ffree
write(*,*)' num=',num
write(*,*)' M=',M
write(*,*)' obtained using code macro_Q2.f90'

kon = 1.0d-02 !1.0d+00             !tPA binding rate. units of inverse (micromolar*sec)

ui = kiss32()
uf = urcw1()

!seed = mscw()
seed= 912309035
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

enddo

!the left and right columns of vertical edges (not including the top-most
!and bottom-most edges)
do i=2,F-1
   closeneigh((3*N-1)*(i-1)+1,(3*N-1)*(i-1)+1-N) = 2
   closeneigh((3*N-1)*(i-1)+1,(3*N-1)*(i-1)+1+(2*N-1)) = 2
   closeneigh((3*N-1)*(i-1)+1,(3*N-1)*(i-1)+1+1) = 4
   
   closeneigh((3*N-1)*(i-1)+1+2*(N-1),(3*N-1)*(i-1)+1+2*(N-1)-(2*N-1)) = 2
   closeneigh((3*N-1)*(i-1)+1+2*(N-1),(3*N-1)*(i-1)+1+2*(N-1)+N) = 2
   closeneigh((3*N-1)*(i-1)+1+2*(N-1),(3*N-1)*(i-1)+1+2*(N-1)-1) = 4
enddo


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
enddo
    
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
enddo
    

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
    enddo
enddo

do i=1,F-1
    do j=2,N-1
        !vertical planar edges

        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*(j-1)) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*j) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*(j-1)+(3*N-1)) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),(i-1)*(3*N-1)+2*j+(3*N-1)) = 1
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),2*N+(j-1)+(i-1)*(3*N-1)-(2*N-1)+(j-1)) = 2
        closeneigh(2*N+(j-1)+(i-1)*(3*N-1),2*N+(j-1)+(i-1)*(3*N-1)-(2*N-1)+(j-1)+(3*N-1)) = 2
  
    enddo
enddo

do i=2,F-1
    do j=2,N-1
        !vertical edges

        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),1+(j-1)*2+(i-1)*(3*N-1)+1) = 2
        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),1+(j-1)*2+(i-1)*(3*N-1)-1) = 2
        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),2*N+(j-2)+(i-2)*(3*N-1)+1) = 2
        closeneigh(1+(j-1)*2+(i-1)*(3*N-1),2*N+(j-2)+(i-2)*(3*N-1)+1+(3*N-1)) = 2
  
    enddo
enddo

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
   enddo
enddo


! read in the data from the micro model, which we obtained from /micro.f90
! READ IN VECTORS FROM MATLAB 

  

  !write(filename2,'(a44)')'tPAleavePLG135_Q2.dat'
  open(200,file='tPAleavePLG135_Q2.dat')
  do i=1,101
     read(200,*, IOSTAT=iostat)CDFtPA(i)
    if( iostat < 0 )then
     write(*,*) 'Warning: File contains less than 101 entries'
     exit
    else if( iostat > 0 )then
     write(*,*) 'Error: error reading file'
     stop
    end if
  end do
  close(200)
  write(*,*)'read tPAleavePLG135_Q2.dat'

  !write(filename3,'(a43)')'tsectPAPLG135_Q2.dat'
  open(300,file='tsectPAPLG135_Q2.dat')
  do i=1,101
     read(300,*)tsec1(i)
  end do
  close(300)
  write(*,*)'read tsectPAPLG135_Q2.dat'



!lysismat_PLG135_Q2.dat is a matrix with column corresponding to bin number (1-100) and with entries
!equal to the lysis times obtained in that bin. an entry of 6000 means lysis didn't happen.
!lysismat(:,1)=the first column, i.e. the lysis times for the first 100 (or 500 if we did 50,000 micro runs) tPA leaving times
    OPEN(unit=1,FILE='lysismat_PLG135_Q2.dat')
    do i=1,nummicro  !100 if only did 10,000 micro runs, 500 if did 50,000
       READ(1,*)(lysismat(i,ii),ii=1,100)
    enddo
    close(1)

!lenlysisvect_PLG135_Q2.dat saves the first row entry in each column of lysismat_PLG135_Q2.dat that lysis
!did not occur, i.e. the first entry there's a 6000
    OPEN(unit=2,FILE='lenlysisvect_PLG135_Q2.dat')
    do i=1,100
        READ(2,*)lenlysismat(i)
    end do
    close(2)
  

  lastmove=0

!The edges that don't contain fibrin are those with edge number <= to enoFB
!enoFB=(3*N-1)*(Ffree-1)

write(*,*)'enoFB=',enoFB

!!!!!! DO "STATS" RUNS OF THE MACRO MODEL
do istat=1,stats

        write(*,*)' run number=',istat

q=0.5d+00!0.2d+00      !0.2             !q is the probability of moving. Make sure it is small enough that we've converged
delx= 1.0135d-04 !10**(-4)             !pore size (distance between nodes), measured in centimeters
Diff= 5.0d-07 !5*10**(-7)           !diffusion coefficient, measured in cm**2/s
tstep=q*delx**2/(12*Diff)  !4/6/11 CHANGED THIS TO (12*Diff) FROM (8*Diff). SEE WRITTEN NOTES 4/6/11 FOR WHY
num_t=tf/tstep            !number of timesteps
bs = 4.27d+02              !concentration of binding sites in micromolar
dist=1.0862d+00 !microns because distance between nodes is 1.0135 micron and diameter of 1 fiber is 0.0727 micron
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


!Initial distribution and boundedness of tPA molecules

!V is matrix of edge loaction and state of boundedness for each tPA. first column is edge 
!molecule is on, second column is 0 if unbound, 1 if bound. start with all molecules unbound

do i=1,M
   V(2,i)=0
enddo

do ii=1,enoFB
   init_state(ii) = dble(ii) / dble(enoFB) !make a vector that's the length of one row of lattice and scale so  
                                                       !that prob. of being on any edge is equal
enddo

call vurcw1(rvect,M)

!use the random numbers to decide where we start the tPAs
do i=1,M
   if (0.le.rvect(i).and.rvect(i).le.init_state(1)) V(1,i)=1 !init_entry(i)=1
   do j=1,enoFB
      if (init_state(j).lt.rvect(i).and.rvect(i).le.init_state(j+1)) then
         !init_entry(i)=j+1
         V(1,i)=j+1
      end if
   enddo
   !V(i,1) = init_entry(i) !molecule i starts on site init_entry(i) so I only allow tPA to start on an edge 
                          !that's in the first row of my lattice
enddo
!    write(*,*)' V=',V  !for debugging 3/31/10


!write(degfile,'(43a)'  ) 'deg_tPA425_PLG35_Q2_test2.dat'
!write(Nfile,'(45a)' ) 'Nsave_tPA425_PLG35_Q2_test2.dat'
!write(tfile,'(45a)') 'tsave_tPA425_PLG35_Q2_test2.dat'
!write(movefile,'(44a)') 'move_tPA425_PLG35_Q2_test2.dat'
!write(lastmovefile,'(48a)') 'lastmove_tPA425_PLG35_Q2_test2.dat'
!write(plotfile,'(44a)') 'plot_tPA425_PLG35_Q2_test2.dat'
!write(degnextfile,'(47a)') 'degnext_tPA425_PLG35_Q2_test2.dat'
!write(Venextfile,'(49a)') 'Vedgenext_tPA425_PLG35_Q2_test2.dat'
!write(Vbdnextfile,'(47a)') 'Vbdnext_tPA425_PLG35_Q2_test2.dat'
!write(cbindfile,'(47a)') 'numbind_tPA425_PLG35_Q2_test2.dat'
!write(cindfile,'(47a)') 'numindbind_tPA425_PLG35_Q2_test2.dat'
!write(bind1file,'(47a)') 'bind_tPA425_PLG35_Q2_test2.dat'
!open(degunit,file=degfile,form=filetype)
!open(Nunit,file=Nfile,form=filetype)
!open(tunit,file=tfile,form=filetype)
!open(moveunit,file=movefile,form=filetype)
!open(lastmoveunit,file=lastmovefile,form=filetype)
!open(plotunit,file=plotfile,form=filetype)
!open(degnextunit,file=degnextfile,form=filetype)
!open(Venextunit,file=Venextfile,form=filetype)
!open(Vbdnextunit,file=Vbdnextfile,form=filetype)
!open(cbindunit,file=cbindfile,form=filetype)
!open(cindunit,file=cindfile,form=filetype)
!open(bind1unit,file=bind1file,form=filetype)

!write(degunit) degrade(:)
!write(tunit) t

write(*,*)' save as deg_tPA425_PLG35_Q2_test2.dat'


Vedgenext(1,:)=V(1,:)
Vboundnext(1,:)=V(2,:)
degnext(1,:)=degrade(:)
tsave(1) = t

call initializeCurses()

!Now do the main part of the code - looping over time


  do

    count=count+1
    t = count*tstep
    if(t>tf) exit  !only do the big "do" loop while t<tf
    

    !if(mod(count,400000)==0) write(*,*)' t=',t
    if(mod(count,ceiling(0.1/tstep))==0) call printGrid(count)

    !at the beginning of each time step, check to see if any of the fibers should be degraded at this time. 
    !Degrade fiber before moving and binding/unbinding tPA:
    do i=enoFB+1,num
       if(t_degrade(i).lt.t.and.t_degrade(i).gt.0.and.degrade(i).eq.0) then
        !i.e. if degradation time is smaller than t AND bigger than 0 AND the edge hasn't already been degraded
            degrade(i)=-t

            !if there were any tPAs bound to this edge, have them unbind, and reset their leaving times to 0:
            do j=1,M
               if(V(1,j)==i) then
                  V(2,j)=0
                  t_leave(j)=0
                  !also find the new binding time for this molecule !FOLLOWING LINE ADDED 4/21/2011:
                  bind(j)=0  !because the molecule can never rebind to a degraded edge
               end if
            enddo
 
        end if
    enddo


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


                !if r400=501, redefine it to be the 500th bin, since we don't have 501 entries
                if(r400==nummicro+1) then
                   r400=nummicro
                end if !for if(r400==501) loop


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

                else     !i.e. if r2 is less than or equal to (t-bind(j))/tstep, have the molecule bind
                  V(2,j)=1  !then the molecule binds
                  bind(j)=0 !reset the binding time to 0
                  r3=urcw1() 

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

                end if !end bind(j) statement
              end if !end r statement        
  
       end if !end V(2,j)=0 statement
    enddo !for j=1,M loop


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

if(Ninteger>Nsave) then !if the current time is the 1st past a new a 10 seconds, e.g. t=10.001, save degrade and V
!write(degunit)    degrade(1:num)
!write(tunit)  t
Nsave=Nsave+10
cNsave=cNsave+1
Vedgenext(cNsave+1,:) = V(1,:) 
Vboundnext(cNsave+1,:) = V(2,:)
degnext(cNsave+1,:) = degrade(1:num)
tsave(cNsave+1) = t
countbindV(istat,cNsave)=countbind
countindepV(istat,cNsave)=countindep
bind1V(istat,cNsave)=sum(bind1)
end if


enddo !for time loop

!   Nsavevect(istat)=cNsave !CHANGED TO CNSAVE FROM NSAVE 12/17/14
!   front=0
!   degold=0

!   !NOW PROCESS THE DATA WE OBTAINED FROM THE ABOVE RUN
!   nplt=Nsavevect(istat)+2 !+1 because I saved once at the beginning, and another +1 to account for the final time point
! write(*,*)'nplt=',nplt
! write(*,*)'r4=',r4

!   do i=2,nplt

!       degold=degnext(i,:)  	! Store the Degredation state of all fibers
!       ind=0
!       place=0
!    do j=1,N 			! For each column
!        do k=1,F-1		! For each row
!            ind(k) = (3*N-1)*(k-1) + 2*N + j-1 !ind is a vector containing the vertical planar edge numbers above node j
!            place(k) = degold(ind(k))  !place(k) is the degradation state of each edge above node j
!        enddo
!        call findfirstreal(place,F-1,0,zero1)  !find the first undegraded vertical edge above node j 
!        front(i-1,j) = zero1	! Store undegraded vertical edge location
!    enddo
!   enddo


!   !the columns of front correspond to the nodes on the x-axis, and the
!   !rows are successive time steps of varying length. The entries are the
!   !y-position of the first undegraded edge at each x location

!   !use this information to calculate the front speed. How? First calculate
!   !the speed at each x-location (so I'll have N speeds. Then average them or
!   !something)

!   !all x-locations start with y=1, so find the first place they deviate from
!   !there:

!   firstdeg=0
!   deglast=0

!   do i=1,N 		! For each column
!       call findfirstineq(front(:,i),tf,1,fdeg) ! Find the first time that the front stops being 1
!       if(fdeg==0) then		! If its always 1
!           firstdeg(i)=1
!       else
!           firstdeg(i)=fdeg	! Store the time the front moves back from 1
!       end if
!   enddo

!   do i=1,N		! For each column
!       call findfirst(front(2:tf,i),tf-1,0,first0)	! Find the first time all fibers are degraded
!       deglast(i)=first0
!   enddo

!   !so firstdeg saves the row # (i.e. time) at which the front first moves and
!   !deglast is the time at which total degradation in a single row occurred.
!   !if deglast=0, then total degradation did not occur

!   moved=0
!   moved(1,:)=firstdeg

!   do j=2,N	! For each column
!       do i=1,N	! For each column
!           if(moved(j-1,i)==0) then  
!               temp=0
!           else
!               call findfirstineq(front(:,i),tf,front(moved(j-1,i),i),temp) ! Find the first time the front changes from the last value we found
!           end if
!           moved(j,i)=temp
!       enddo
!   enddo

!   !so now "moved" saves the saved-time-step at which the front moves for each x location

!   do i=1,N
!       call findintineq(moved(:,i),N,0,lasti) ! Find the last time the front moved in column i
!       lastmove(i,istat)=lasti
!   enddo

!   plotstuff=0
!   plotstuff2=0

!   do j=1,N
!       do i=1,lastmove(j,istat)
!           plotstuff(j,i)=front(moved(i,j),j) 		! Collapse down the moved() matrix
!           plotstuff2(j,i)=(plotstuff(j,i)-1)*dist	! Convert from coordinates to distance
!       enddo
!   enddo
    
!   !now plotstuff has in each row the successive y-positions of x-location
!   !corresponding to row number, and plotstuff2 has in each row the successive
!   !y-positions (in microns, instead of node #) of x-location corresponding to
!   !row number

!   !In order to plot this and finish the calculations in Matlab, I need to save plotstuff2, lastmove, and moved

!   !write(moveunit) moved(:,:)
!   !write(plotunit) plotstuff2(:,:)

!   if(istat==1)then  !choose how many runs you want to save to make a movie
!      !write(degnextunit) degnext(:,:)
!      !write(Venextunit) Vedgenext(:,:)
!      !write(Vbdnextunit) Vboundnext(:,:)

!      !!!DO MORE MOVIE PROCESSING BEFORE GOING TO MATLAB
      
!       X1plot=0
!       Y1plot=0
!       X2plot=0
!       Y2plot=0
!       Xvplot=0
!       Yvplot=0
!       freetPA=0
!       bdtPA=0

!      countintact2=0
!      counth=0
!      countv=0
!      countpv=0
!      do i=enoFB,num
!         if(degnext(1,i)==0) then
!           countintact2=countintact2+1
!           intact2(countintact2)=i   !finds the undegraded edge numbers
!         end if
!      enddo

!      lenintact2=countintact2

!      !Assume you have a grid with N nodes in each row of the lattice, and you
!      !assign the nodes a number simply by counting nodes, starting at the bottom
!      !left and moving right, then going to the left of row two and moving right,
!      !etc. in column K, endpts has the node numbers corresponding to the
!      !endpts of fiber (i.e. edge) K. K=1,...,num

!      !for vertical edges:
!      do i=1,F
!          do j=1,N
!              endpts(1,(3*N-1)*(i-1)+1+2*(j-1)) = (i-1)*N + j
!              endpts(2,(3*N-1)*(i-1)+1+2*(j-1)) = (i-1)*N + j
!          end do
!      end do

!      !for horizontal edges
!      do i=1,F
!          do j=1,N-1
!              endpts(1,(3*N-1)*(i-1)+2+2*(j-1)) = j + N*(i-1)
!              endpts(2,(3*N-1)*(i-1)+2+2*(j-1)) = j + N*(i-1) + 1
!          end do
!      end do

!      !for planar veritcal edges
!      do i=1,F-1
!          do j=1,N
!              endpts(1,(3*N-1)*(i-1)+2*N+(j-1)) = j + N*(i-1)
!              endpts(2,(3*N-1)*(i-1)+2*N+(j-1)) = j + N*(i-1) + N
!          end do
!      end do


!      do jj=1,lenintact2
!          !horizontal edges
!          do iplt=Ffree,F
!            do j=1,N-1
!              if(intact2(jj)==2*(j-1)+2+(3*N-1)*(iplt-1)) then    !if we have a horizontal edge
!                 yplace=endpts(1,intact2(jj))/N+1  !find the y value at which the horizontal edge occurs
!                 x2=nint((real(endpts(2,intact2(jj)))/real(N)-floor(real(endpts(2,intact2(jj)))/real(N)))*N)
!                 if(x2==0) x2=N    !if it says the RHS endpoint is 0, force it to actually be N (because otherwise 
!                                   !is says we should plot from N-1 to 0)
!                 counth=counth+1
!                 X1plot(1,counth)=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N) 
!                 X1plot(2,counth)=x2
!                 Y1plot(1,counth)=yplace
!                 Y1plot(2,counth)=yplace
!              end if
!            enddo
!         enddo
!         !vertical (planar) edges
!         do j=Ffree,F-1
!           do k=1,N
!              if(intact2(jj)==(3*N-1)*(j-1)+2*N+(k-1)) then   !if we have a vertical (planar) edge
!                xplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N) 
!                                                                 !find the x value at which the vertical edge occurs
!                y1=endpts(1,intact2(jj))/N+1  !find the bottom endpoint of the vertical edge
!                y2=endpts(2,intact2(jj))/N+1
!                if(xplace==0) then
!                  xplace=N
!                  y1=endpts(1,intact2(jj))/N
!                  y2=endpts(2,intact2(jj))/N
!                end if
!                countpv=countpv+1
!                X2plot(1,countpv)=xplace
!                X2plot(2,countpv)=xplace
!                Y2plot(1,countpv)=y1
!                Y2plot(2,countpv)=y2
!             end if
!          end do
!        end do
!        !vertical edges
!        do i=Ffree,F
!          do j=1,N
!             if(intact2(jj)==(3*N-1)*(i-1)+1+2*(j-1)) then  !if we have a vertical edge
!               yvplace=(endpts(1,intact2(jj))/N)+1  !find the y value at which the vertical edge occurs
!               xvplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)  
!                                                                       !find the x value at which the vertical edge occurs
              
!               if(xvplace==0) then
!                 xvplace=N
!                 yvplace=endpts(1,intact2(jj))/N
!               end if
!               countv=countv+1
!               Xvplot(countv)=xvplace
!               Yvplot(countv)=yvplace
!             end if
!          end do
!       end do
!     end do  !for jj loop

! !write(x1file,'(46a)'  ) 'X1plot_tPA425_PLG35_Q2_test2.dat'
! !open(x1unit,file=x1file,form=filetype)
! !write(x2file,'(46a)'  ) 'X2plot_tPA425_PLG35_Q2_test2.dat'
! !open(x2unit,file=x2file,form=filetype)
! !write(y1file,'(46a)'  ) 'Y1plot_tPA425_PLG35_Q2_test2.dat'
! !open(y1unit,file=y1file,form=filetype)
! !write(y2file,'(46a)'  ) 'Y2plot_tPA425_PLG35_Q2_test2.dat'
! !open(y2unit,file=y2file,form=filetype)
! !write(xvfile,'(46a)'  ) 'Xvplot_tPA425_PLG35_Q2_test2.dat'
! !open(xvunit,file=xvfile,form=filetype)
! !write(yvfile,'(46a)'  ) 'Yvplot_tPA425_PLG35_Q2_test2.dat'
! !open(yvunit,file=yvfile,form=filetype)

! !write(x1unit) X1plot
! !write(x2unit) X2plot
! !write(y1unit) Y1plot
! !write(y2unit) Y2plot
! !write(xvunit) Xvplot
! !write(yvunit) Yvplot

! !!Now do location and boundedness of tPA
!      do i=1,F
!         do j=1,N-1
!            !horizontal edges
!            do jplt=1,M
!               if(Vedgenext(1,jplt)==2*(j-1)+2+(3*N-1)*(i-1)) then
!                   Vedgeplace=(endpts(1,Vedgenext(1,jplt))/N)+1
!                   Vx=nint((real(endpts(2,Vedgenext(1,jplt)))/real(N)-floor(real(endpts(2,Vedgenext(1,jplt)))/real(N)))*N)
!                   if(Vx==0) then
!                      Vx=N
!                   end if
!                   if(Vboundnext(1,jplt)==0) then   !if unbound, will plot in black
!                      bdtPA(1,jplt)=-1
!                      bdtPA(2,jplt)=-1
!                      freetPA(1,jplt)=Vx-0.5d+00
!                      freetPA(2,jplt)=Vedgeplace
!                   elseif(Vboundnext(1,jplt)==1) then  !if bound, will plot in green
!                      bdtPA(1,jplt)=Vx-0.5d+00
!                      bdtPA(2,jplt)=Vedgeplace
!                      freetPA(1,jplt)=-1
!                      freetPA(2,jplt)=-1
!                   end if
!               end if
!            end do
!        end do
!      end do
!      do j=1,F-1
!         do k=1,N
!           !vertical (planar) edges
!            do kplt=1,M
!               if(Vedgenext(1,kplt)==(3*N-1)*(j-1)+2*N+(k-1)) then
!                  xVedgeplace=nint((real(endpts(1,Vedgenext(1,kplt)))/real(N)-floor(real(endpts(1,Vedgenext(1,kplt)))/real(N)))*N)  !find the x value at which the vertical edge occurs
!                  Vy1=endpts(1,Vedgenext(1,kplt))/N+1  !find the bottom endpoint of the vertical edge
!                  if(xVedgeplace==0) then
!                     xVedgeplace=N
!                     Vy1=endpts(1,Vedgenext(1,kplt))/N
!                  end if
!                  if(Vboundnext(1,kplt)==0) then   !if unbound, plot in black
!                     bdtPA(1,kplt)=-1
!                     bdtPA(2,kplt)=-1
!                     freetPA(1,kplt)=xVedgeplace
!                     freetPA(2,kplt)=Vy1+0.5d+00
!                  elseif(Vboundnext(1,kplt)==1) then  !if bound, plot in green
!                     bdtPA(1,kplt)=xVedgeplace
!                     bdtPA(2,kplt)=Vy1+0.5d+00
!                     freetPA(1,kplt)=-1
!                     freetPA(2,kplt)=-1
!                  end if
!                end if
!             end do
!         end do
!      end do
!      do i=1,F
!         do j=1,N
!          !vertical edges
!           do kjplt=1,M
!              if(Vedgenext(1,kjplt)==(3*N-1)*(i-1)+1+2*(j-1)) then
!                 vertplace=1+(j-1)   !find the x value at which the vertical edge occurs
!                 Vyvert=endpts(1,Vedgenext(1,kjplt))/N+1  !find the bottom endpoint of the vertical edge
!                 if(vertplace==N) then
!                   !vertplace=N;
!                    Vyvert=endpts(1,Vedgenext(1,kjplt))/N
!                 end if
!                 if(Vboundnext(1,kjplt)==0) then  !if unbound, plot in black
!                    bdtPA(1,kjplt)=-1
!                    bdtPA(2,kjplt)=-1
!                    freetPA(1,kjplt)=vertplace
!                    freetPA(2,kjplt)=Vyvert
!                 elseif(Vboundnext(1,kjplt)==1) then   !if bound, plot in green
!                    bdtPA(1,kjplt)=vertplace
!                    bdtPA(2,kjplt)=Vyvert
!                    freetPA(1,kjplt)=-1
!                    freetPA(2,kjplt)=-1
!                 end if
!               end if
!            end do
!         end do
!      end do


! !write(tPAbdfile,'(46a)'  ) 'tPAbd_tPA425_PLG35_Q2_test2.dat'
! !open(tPAbdunit,file=tPAbdfile,form=filetype)
! !write(tPAfreefile,'(47a)'  ) 'tPAfree_tPA425_PLG35_Q2_test2.dat'
! !open(tPAfreeunit,file=tPAfreefile,form=filetype)

! !write(tPAbdunit) bdtPA
! !write(tPAfreeunit) freetPA


! !now save different timestep so I can make a matlab movie
!    do imod=2,nplt
!       if(mod(imod,60)==0) then  !this means we plot approximately every minute

!          X1plot=0
!          Y1plot=0
!          X2plot=0
!          Y2plot=0
!          Xvplot=0
!          Yvplot=0
!          freetPA=0
!          bdtPA=0         

!          intact2=0
!          countintact2=0
!          counth=0
!          countv=0
!          countpv=0
!          do i=enoFB,num  !could do this as"do i=enoFB,num" since we know the first enoFB edges won't equal 0
!            if(degnext(imod,i)==0) then
!              countintact2=countintact2+1
!              intact2(countintact2)=i   !finds the undegraded edge numbers
!            end if
!          enddo

!          lenintact2=countintact2 
         
!       do jj=1,lenintact2
!          !horizontal edges
!          do iplt=Ffree,F
!            do j=1,N-1
!              if(intact2(jj)==2*(j-1)+2+(3*N-1)*(iplt-1)) then    !if we have a horizontal edge
!                 yplace=endpts(1,intact2(jj))/N+1  !find the y value at which the horizontal edge occurs
!                 x2=nint((real(endpts(2,intact2(jj)))/real(N)-floor(real(endpts(2,intact2(jj)))/real(N)))*N)  
!                          if(x2==0) x2=N    !if it says the RHS endpoint is 0, force it to actually be N (because otherwise 
!                                            !is says we should plot from N-1 to 0)
!                          counth=counth+1
!                          X1plot(1,counth)=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)
!                          X1plot(2,counth)=x2
!                          Y1plot(1,counth)=yplace
!                          Y1plot(2,counth)=yplace
!              end if
!            enddo
!         enddo
!         !vertical (planar) edges
!         do j=Ffree,F-1
!           do k=1,N
!              if(intact2(jj)==(3*N-1)*(j-1)+2*N+(k-1)) then   !if we have a vertical (planar) edge
!                xplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N) 
!                                                                 !find the x value at which the vertical edge occurs
!                y1=endpts(1,intact2(jj))/N+1  !find the bottom endpoint of the vertical edge
!                y2=endpts(2,intact2(jj))/N+1
!                if(xplace==0) then
!                  xplace=N
!                  y1=endpts(1,intact2(jj))/N
!                  y2=endpts(2,intact2(jj))/N
!                end if
!                countpv=countpv+1
!                X2plot(1,countpv)=xplace
!                X2plot(2,countpv)=xplace
!                Y2plot(1,countpv)=y1
!                Y2plot(2,countpv)=y2
!             end if
!          end do
!        end do
!        !vertical edges
!        do i=Ffree,F
!          do j=1,N
!             if(intact2(jj)==(3*N-1)*(i-1)+1+2*(j-1)) then  !if we have a vertical edge
!               yvplace=(endpts(1,intact2(jj))/N)+1  !find the y value at which the vertical edge occurs
!               xvplace=nint((real(endpts(1,intact2(jj)))/real(N)-floor(real(endpts(1,intact2(jj)))/real(N)))*N)  
!                                                                       !find the x value at which the vertical edge occurs
!               if(xvplace==0) then
!                 xvplace=N
!                 yvplace=endpts(1,intact2(jj))/N
!               end if
!               countv=countv+1
!               Xvplot(countv)=xvplace
!               Yvplot(countv)=yvplace
!             end if
!          end do
!       end do
!     end do  !for jj loop
  
!   !write(x1unit) X1plot
!   !write(x2unit) X2plot
!   !write(y1unit) Y1plot
!   !write(y2unit) Y2plot
!   !write(xvunit) Xvplot
!   !write(yvunit) Yvplot

! !!Now do location and boundedness of tPA
!      do i=1,F
!         do j=1,N-1
!            !horizontal edges
!            do jplt=1,M
!               if(Vedgenext(imod,jplt)==2*(j-1)+2+(3*N-1)*(i-1)) then
!                   Vedgeplace=(endpts(1,Vedgenext(imod,jplt))/N)+1
!                   Vx=nint((real(endpts(2,Vedgenext(imod,jplt)))/real(N)-floor(real(endpts(2,Vedgenext(imod,jplt)))/real(N)))*N)
!                   if(Vx==0) then
!                      Vx=N
!                   end if
!                   if(Vboundnext(imod,jplt)==0) then   !if unbound, will plot in black
!                      bdtPA(1,jplt)=-1
!                      bdtPA(2,jplt)=-1
!                      freetPA(1,jplt)=Vx-0.5d+00
!                      freetPA(2,jplt)=Vedgeplace
!                   elseif(Vboundnext(imod,jplt)==1) then  !if bound, will plot in green
!                      bdtPA(1,jplt)=Vx-0.5d+00
!                      bdtPA(2,jplt)=Vedgeplace
!                      freetPA(1,jplt)=-1
!                      freetPA(2,jplt)=-1
!                   end if
!               end if
!            end do
!        end do
!      end do
!      do j=1,F-1
!         do k=1,N
!           !vertical (planar) edges
!            do kplt=1,M
!               if(Vedgenext(imod,kplt)==(3*N-1)*(j-1)+2*N+(k-1)) then
!                  xVedgeplace=nint((real(endpts(1,Vedgenext(imod,kplt)))/real(N)-floor(real(endpts(1,Vedgenext(imod,kplt)))/real(N)))*N)  !find the x value at which the vertical edge occurs
!                  Vy1=endpts(1,Vedgenext(imod,kplt))/N+1  !find the bottom endpoint of the vertical edge
!                  if(xVedgeplace==0) then
!                     xVedgeplace=N
!                     Vy1=endpts(1,Vedgenext(imod,kplt))/N
!                  end if
!                  if(Vboundnext(imod,kplt)==0) then   !if unbound, plot in black
!                     bdtPA(1,kplt)=-1
!                     bdtPA(2,kplt)=-1
!                     freetPA(1,kplt)=xVedgeplace
!                     freetPA(2,kplt)=Vy1+0.5d+00
!                  elseif(Vboundnext(imod,kplt)==1) then  !if bound, plot in green
!                     bdtPA(1,kplt)=xVedgeplace
!                     bdtPA(2,kplt)=Vy1+0.5d+00
!                     freetPA(1,kplt)=-1
!                     freetPA(2,kplt)=-1
!                  end if
!                end if
!             end do
!         end do
!      end do
!      do i=1,F
!         do j=1,N
!          !vertical edges
!           do kjplt=1,M
!              if(Vedgenext(imod,kjplt)==(3*N-1)*(i-1)+1+2*(j-1)) then
!                 vertplace=1+(j-1)   !find the x value at which the vertical edge occurs
!                 Vyvert=endpts(1,Vedgenext(imod,kjplt))/N+1  !find the bottom endpoint of the vertical edge
!                 if(vertplace==N) then
!                   !vertplace=N;
!                    Vyvert=endpts(1,Vedgenext(imod,kjplt))/N
!                 end if
!                 if(Vboundnext(imod,kjplt)==0) then  !if unbound, plot in black
!                    bdtPA(1,kjplt)=-1
!                    bdtPA(2,kjplt)=-1
!                    freetPA(1,kjplt)=vertplace
!                    freetPA(2,kjplt)=Vyvert
!                 elseif(Vboundnext(imod,kjplt)==1) then   !if bound, plot in green
!                    bdtPA(1,kjplt)=vertplace
!                    bdtPA(2,kjplt)=Vyvert
!                    freetPA(1,kjplt)=-1
!                    freetPA(2,kjplt)=-1
!                 end if
!               end if
!            end do
!         end do
!      end do
  
!   !write(tPAbdunit) bdtPA
!   !write(tPAfreeunit) freetPA

!    end if !for if mod(imod,60) loop 
!   end do !for imod loop

! !close(x1unit)
! !close(x2unit)
! !close(y1unit)
! !close(y2unit)
! !close(xvunit)
! !close(yvunit)
! !close(tPAbdunit)
! !close(tPAfreeunit)


! !!!!! END ADDED STUFF FOR MOVIE
!   end if


 

!  countindepV(istat,tf)=countindep
 enddo  !for stats loop

! write(*,*)'Nsavevect=',Nsavevect(:)

!write(cbindunit) countbindV
!write(cindunit) countindepV
!write(bind1unit) bind1V
!write(Nunit) Nsavevect(:)
!write(lastmoveunit) lastmove(:,:)

!close(degunit)
!close(Nunit)
!close(tunit)
!close(moveunit)
!close(lastmoveunit)
!close(plotunit)
!close(degnextunit)
!close(Venextunit)
!close(Vbdnextunit)
!close(cbindunit)
!close(cindunit)
!close(bind1unit)

ierr = getch()
ierr = endwin()

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!

subroutine initializeCurses()
  stdscr = initscr()
  ierr = noecho()
  ierr = curs_set(0)
  ierr = start_color()
  ierr = init_pair(black2, COLOR_BLACK, COLOR_BLACK)
  ierr = init_pair(red2, COLOR_RED, COLOR_BLACK)
  ierr = init_pair(green2, COLOR_GREEN, COLOR_BLACK)
  ierr = init_pair(yellow2, COLOR_YELLOW, COLOR_BLACK)
  ierr = init_pair(blue2, COLOR_BLUE, COLOR_BLACK)
  ierr = init_pair(white2, COLOR_WHITE, COLOR_BLACK)
end subroutine initializeCurses


subroutine printGrid(t)
    integer :: t
    integer(C_INT) :: maxX,maxY
    integer(C_INT) :: x,y
    real :: currentTime


    currentTime = t*tstep
    call getmaxyx(stdscr, maxY, maxX)
    maxX = N*(gridSize*2 + 2) + 10
    ierr = attron(COLOR_PAIR(red4))
    ierr = mvprintw(10, maxX, "|"//C_NULL_CHAR)
    ierr = attron(COLOR_PAIR(blue4))
    ierr = mvprintw(12, maxX, "|"//C_NULL_CHAR)
    ierr = attron(COLOR_PAIR(green4))
    ierr = mvprintw(14, maxX, "*"//C_NULL_CHAR)
    ierr = attron(COLOR_PAIR(yellow4))
    ierr = mvprintw(16, maxX, "*"//C_NULL_CHAR)
    ierr = attron(COLOR_PAIR(white4))
    ierr = mvprintw(0, 0, "Time: %f4.3 secs"//C_NULL_CHAR, currentTime)
    !ierr = r_mvprintw(0, 10, C_NULL_CHAR, currentTime)
    !ierr = mvprintw(0, 20, "secs"//C_NULL_CHAR)
    ierr = mvprintw(10, maxX+2, "Degraded Fiber"//C_NULL_CHAR)
    ierr = mvprintw(12, maxX+2, "Non-degraded Fiber"//C_NULL_CHAR)
    ierr = mvprintw(14, maxX+2, "Unbound Molecule"//C_NULL_CHAR)
    ierr = mvprintw(16, maxX+2, "Bound Molecule"//C_NULL_CHAR)
    !mvprintw(0, 0, ("Time: " + to_string(t * timeStep) + " sec");
    do x=1,N
      do y=1,F
        !call printNode(x,y);
      enddo
    enddo
    do j=1,M
      !call printMolecule(j)
    enddo
    ierr = refresh()
    ierr = getch()

end subroutine printGrid

!
!Returns the (single-dimensional) index of the fiber 
!leaving node (nodeX, nodeY) in direction 'fiber'.
!If an impossible request is made, the method returns '-1'.
!NOTE: All of these indices are 0-indexed.
!
INTEGER function fiberIndex(nodeX, nodeY, direction)
    integer :: nodeX, nodeY, direction
    if ((nodeX < 0) || (nodeY < 0))                        
        write(*,*)"Index must be non-negative."
    else if (nodeX >= nodesInRow)                          
        write(*,*)"Index out of bounds."
    else if (nodeY >= nodesInColumn)                       
        write(*,*)"Index out of bounds."
    else
      ! This section only calculates the fiber index for UP, RIGHT, and OUT. 
      ! If one of the other directions are requested, it moves to a neighboring node
      ! and converts the request.
      ! i.e., if DOWN is requested, we recursively call the method again 
      ! asking for the UP fiber from the node below.
        switch (direction) {                            // Determine the direction of the fiber requested
            case DOWN:
                return (nodeY == 0) ? -1 : fiberIndex(nodeX, nodeY-1, UP);      // Ensure that we are not at the bottom edge of the grid, then recurse
            case LEFT:                                                     
                return (nodeX == 0) ? -1 : fiberIndex(nodeX-1, nodeY, RIGHT);   // Ensure that we are not on the left edge of the grid, then recurse
            case IN:
                return fiberIndex(nodeX, nodeY, OUT);                 // Every node has an IN fiber, so just recurse
            case UP:
                return (nodeY == nodesInColumn - 1) ? -1 :              // Ensure that we are not on the top edge of the grid
                    nodeY * fullRow +                           // The index is all fibers from all rows below, 
                     xzRow +                              // plus the x- and z-fibers from the current row,
                     nodeX;                               // plus the y-fibers for all nodes to the left in this row.
            case RIGHT:
                return (nodeX == nodesInRow - 1) ? -1 :               // Ensure that we are not on the right edge of the grid
                    nodeY * fullRow +                           // The index is all fibers from all rows below,
                     nodeX * 2 +                            // plus all the x- and z-fibers for all nodes to the left in this row,
                     1;                                 // plus the z-fiber for the current node
            case OUT:
                return nodeY * fullRow +                      // The index is all fibers from all rows below,
                 nodeX * 2;                             // plus all the x- and z-fibers for all nodes to the left in this row.
            default:
                return -1;
        }
}

/*
 * Returns the x-coordinate of the node at the bottom/right endpoint
 * of a fiber given the (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short nodeX(unsigned short fiberIndex) {
    if (fiberIndex < 0)                               // Ensures that the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers)                       // Ensures that the index is in bounds
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;                // Find the fiber's position in its row 
    if (rowPosition >= xzRow)                             // If this fiber is in the second (y-fiber) part of its row
      return rowPosition - xzRow;                           // Just subtract the number of x- and z-fibers to get the node index
    else                                      // Else if this fiber is in the first (x-, z-fiber) part of its row
        return (unsigned short)(rowPosition / 2);                   // Its node index is half its fiber index (rounded/cast) down
}

/*
 * Returns the y-coordinate of the node at the bottom/right endpoint
 * of a fiber given the (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short nodeY(unsigned short fiberIndex) {
    if (fiberIndex < 0)                               // Ensures that the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers)                       // Ensures that the index is in bounds
        throw invalid_argument("Index out of bounds.");
    return (unsigned short)(fiberIndex / fullRow);                  // Determine the number of full rows preceding this fiber
}

/*
 * Returns the direction of a fiber given the 
 * (single-dimensional) index of a fiber.
 * NOTE: All of these indices are 0-indexed.
 */
unsigned short fiberDirection(unsigned short fiberIndex) {
    if (fiberIndex < 0)                               // Ensures that the index is non-negative
        throw invalid_argument("Index must be non-negative.");
    else if (fiberIndex >= totalFibers)                       // Ensures that the index is in bounds
        throw invalid_argument("Index out of bounds.");
    unsigned short rowPosition = fiberIndex % fullRow;                // Find the fiber's position in its row 
    if (rowPosition >= xzRow)                             // If this fiber is in the second (y-fiber) part of its row
      return UP;                                    // It must be an UP fiber
    else if (rowPosition % 2 == 0)                          // Else, the even fibers are OUT and the odd ones are RIGHT
      return OUT;
    else
      return RIGHT; 
}

end program macrolysis

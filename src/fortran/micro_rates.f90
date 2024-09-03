program micromodel

    !! BRAD 2024-01-13: This code has been modified in the following ways:
    !!                  - Data folder is relative to git repository root
    !!                  - Data is stored in subfolders based on expCode
    !!                  - Data file codes are now set globally from the (in/out)FileCode variables

    implicit none
    character(15) :: expCode = '2024-01-13-0710'
    character(6)  :: inFileCode = 'Q2.dat'
    character(40)   :: outFileCode = 'PLG2_tPA01_Q2.dat'
    !!!! This code is the microscale model with lots of opportunities for changing the rate constants and initial concentrations
    !!!! Lines 19-25 allow you to set the various dissociation constants, binding rates, and the concentration of free PLG
    !!!! This code treats degradation and exposure in the gillespie algorithm, rather than separately with
    !!!! a degradation timer. It also allows PLi to degrade any exposed doublet at the same binding location, and
    !!!! tPA to convert any PLG on the same binding location to PLi.


    !! Q0      46 nm diameter fibers,     4 nodes per row,  fibrin concentration 710.78 uM, binding site concentration 355.390 uM
    !! Q1      57.4 nm diameter fibers,   5 nodes per row, 
    !! Q2      72.7 nm diameter fibers,   7 nodes per row, fibrin concentration 871.17 uM, binding site concentration 426.848 uM
    !! Q3      81.3 nm diameter fibers,   8 nodes per row
    !! TF-v    105.1 nm diameter fibers,  5 nodes per row, fibrin concentration 212.75 uM, binding site concentration 136.159 uM
    !! TF-vii  105.1 nm diameter fibers,  7 nodes per row, fibrin concentration 416.84 uM, binding site concentration 204.238 uM
    !! TF-x    105.1 nm diameter fibers, 10 nodes per row, fibrin concentration 850.69 uM, binding site concentration 306.357 uM
    !! TB-xi   123.0 nm diameter fibers, 11 nodes per row, fibrin concentration 751.54 uM, binding site concentration 248.531 uM
    !! Q4      145.4 nm diameter fibers, 13 nodes per row, fibrin concentration 751.16 uM, binding site concentration 213.0 uM
    !!

    integer  :: nodes = 7 !total number of nodes in one row of the lattice. This is the only difference between thin and thick runs, so it is the only change that must be made. 5 for Q1 (57.4 nm), 7 for Q2 (72.7 nm), 8 for Q3 (81.3 nm), 13 for Q4 (145.4 nm)
    double precision  :: radius = 72.7/2/1000 ! fiber bundle radius in microns !! MAKE SURE you enter this as a decimal, even if it ends in .0
    integer, parameter :: Nplginit = 1 !number of exposed doublets initially - i.e. intact doublets - at each spatial location
    integer, parameter  :: Ninit = 5*Nplginit !number of cryptic doublets at each spatial location
    integer, parameter  :: Ntot = Ninit + Nplginit !total number of doublets at each spatial location
    integer  :: runs = 50000 !How many independent simulations to run of the microscale model
    integer, parameter  :: sV = 26  !number of reactions, i.e. size of V
    integer          :: stats

    !Define the Kd's and on rates that will be used in the given run. These will be used to help define the off rates later.
    double precision  :: KdtPAyesplg = 0.02 !0.02 !units uM, tPA Kd in presence of PLG
    double precision  :: KdtPAnoplg = 0.36 !0.36 !units uM, tPA Kd in absence of PLG
    double precision  :: KdPLGintact = 38 !10 !units uM, PLG Kd to intact fibrin !38 in original model
    double precision  :: KdPLGnicked = 2.2 !1 !units uM, PLG Kd to nicked fibrin !2.2 in original model
    double precision  :: ktPAon = 0.1 !0.1 !units 1/(uM*s), tPA binding rate to fibrin !0.01 in original model
    double precision  :: kplgon = 0.1 !units 1/(uM*s), PLG bindind rate to fibrin
    double precision  :: freeplg = 2 !1.5 !units uM, concentration of free PLG

    double precision  :: kdeg = 5 !units 1/s, plasmin-mediated degradation rate
    double precision  :: kplioff = 57.6 !units 1/s, PLi unbinding rate
    double precision  :: kapcat = 0.1 !units 1/s, tPA-mediated rate of conversion of PLG to PLi
    double precision  :: kncat = 5 !units 1/s, PLi-mediated rate of exposure of new binding sites

    double precision  :: kplgoff  !units 1/s, off rate for intact FB
    double precision  :: kplgoffnick !units 1/(s*uM), off rate for nicked FB
    double precision  :: kaoff12 !units 1/s, off rate in presence of PLG
    double precision  :: kaoff10 !units 1/s, off rate in absense of PLG

    !double precision  :: prob_N02
    !double precision  :: prob_N00
    !double precision  :: prob_N22

    !! BRAD 2024-01-14: Change default to zero later
    integer :: seed = 0 ! 981681759

    !double precision, dimension(11,2)  :: param  !matrix that holds all the various parameter values we can use

    integer, dimension(:, :), allocatable  :: state, statetemp !matrices to save the state of each doublet. There are 6 doublets at each node, and nodes^2 total nodes. recall: in fortran, columns are listed 1st, rows 2nd
    double precision, dimension(:), allocatable  :: init_state, cumsuminit
    integer  :: init_entry
    integer, dimension(:, :), allocatable  :: Lat !matrix of the connections between nodes. For instance, if node 1 is not a direct neighbor of node 6, then there would be a 0 in the (1,6) and (6,1) entries of Lat
    !! BRAD 2024-01-13:
    ! integer, dimension(nodes**2,nodes**2)  :: Lat_temp !matrix of the connections between nodes. For instance, if node 1 is not a direct neighbor of node 6, then there would be a 0 in the (1,6) and (6,1) entries of Lat

    double precision  :: r, rin, randdegexp
    integer  :: Tdoublets
    double precision  :: Tdoublets2
    double precision  :: kcat

    integer, dimension(9, sV) :: V !stoichiometric matrix. each row corresponds to a different reaction and each column represents a different state
    integer, dimension(9*sV)  :: V_tmp !temporary stoichiometric matrix (it gets reshaped into V). each row corresponds to a different reaction and each column represents a different state
    double precision  :: c2, c3, c5, c6, c8, c8a, c9
    double precision  :: c11, c12, c13, c14, c15, c20, c25
    integer, dimension(9, 9)  :: q
    double precision  :: t, tfinal
    double precision  :: blahtest, kunbind
    integer  :: count, countPLi, countstate, countstate2, countstate3, countstate122
    integer  :: countlat, sizenonzero_connectind, countstate4, countstate5, countstate4a, countstate5a
    integer  :: countstate6, countstate7, countstate8, countstate9, countstate10, countlat2, countlat3
    integer  :: countstate11, countstate12, countstate13, countstate14, countstate15, countstate10a, countstate11a
    integer  :: countstate11b, cccc
    integer  :: plotnumber
    double precision, dimension(1000000)  :: tvals !vector of times at which reactions happened. These vectors are set to 1000000 because I don't know how many timesteps will be taken during the Gillespie Algorithm
    integer, dimension(1000000)  :: Anumber
    double precision, dimension(300)  :: PLideg, tau2
    integer, dimension(300)  :: ones
    integer, dimension(2, 300)  :: PLilocation, PLilocationt
    integer, dimension(2)      :: PLilocsmall
    integer, dimension(300)  :: PLiwaste
    integer, dimension(1000000)  :: tPA, tPA_loc
    double precision, dimension(1000000)  :: PLG
    integer, dimension(1000000)  :: PLi
    integer, dimension(1000000)  :: per_lysis
    double precision, dimension(1000000)  :: timevect
    double precision, dimension(1000000)  :: percent_degrade, percent_degrade2
    integer, dimension(1000000)  :: degrade_tot
    integer, dimension(1000000)  :: undegrade_tot
    integer, dimension(1000000)  :: undegrade_new
    integer, dimension(1000000)  :: PLi_loc
    integer  :: temp

    ! integer, dimension(nodes**2,nodes**2)  :: num_degrade !I don't think this gets used anymore
    ! integer, dimension(nodes**2,nodes**2)  :: num_undegrade !I don't think this gets used anymore
    double precision  :: tnext
    double precision  :: D
    double precision  :: vol3
    double precision  :: bs, bs0
    integer, dimension(300)  :: row, col, rowPLi, colPLi
    integer  :: A
    double precision, dimension(12000)  :: tau !these vectors are likely WAY bigger than necessary. They only need to be as big as the total number of "active" sites in the fiber, i.e., doublets that have tPA and/or plasmin on them, or are degraded. Likely this would never be bigger than 2*6*nodes^2... Could define an integer parameter like "sizeGillespievector" and make all of these have dimension "sizeGillespievector". Then we could play around with the actual value in a simpler way than manually having to change all of these each time
    double precision, dimension(12000, sV)  :: aa
    double precision, dimension(sV)  :: cumsumaa, asumtemp
    double precision, dimension(12000, 9)  :: qstate
    double precision, dimension(12000)  :: asum
    integer, dimension(12000)  :: bigJ
    integer  :: newsite, newsite2
    integer  :: newdoublet, newdoublet2
    double precision, dimension(5)  :: P, cumP
    double precision, dimension(2)  :: Pnsmall, cumPnsmall, Pismall, cumPismall
    integer, dimension(2)  :: currentPLib
    integer  :: currentPLi
    integer  :: countcurPLi
    integer  :: i, ii, j, jj, iPLi2, iPLi !, ipar
    double precision, dimension(12300)  :: tautot
    integer    :: place, loc, loc2b
    double precision  :: entry
    integer  :: loc2, reaction2, loca, reaction3
    integer, dimension(6)  :: connect_ind, connect_ind2temp
    integer, dimension(5, 6)  :: connect_ind2
    double precision  :: Fi, Fn, Ftot   !I want these to be integers, but I eventually divide them, so it doesn't work...
    double precision  :: prob_N02i
    double precision  :: prob_N00i
    double precision  :: prob_N22i
    double precision  :: prob_N02n
    double precision  :: prob_N00n
    double precision  :: prob_N22n

    double precision, dimension(6)  :: tot_sitesPLi
    double precision, dimension(6)  :: prob_sitePLi, cumprob
    integer, dimension(6)  :: sites_N03
    integer, dimension(6)  :: sites_N10
    integer, dimension(6)  :: sites_Nplg
    integer, dimension(6)  :: sites_Nplgn
    integer, dimension(6)  :: sites_Nplgi
    integer, dimension(6)  :: sites_degrade
    integer, dimension(5)  :: Nplgn

    double precision  :: alpha
    integer  :: statePLi
    integer  :: numunexp
    integer  :: degexp
    integer  :: NtoNplg
    integer  :: NplgtoN12
    integer  :: N03toN13
    integer  :: N10toN13
    integer  :: NplgtoN23
    integer  :: N03toN33
    integer  :: degradetoPLi

    double precision  :: exposed
    double precision  :: r0
    double precision  :: r1
    double precision  :: k0
    double precision  :: p_rebind
    double precision  :: rr
    double precision, dimension(2, 1000000)  :: Prob_rebind
    integer :: rrebind, crebind, newrow, newcol, indnewsite
    integer, dimension(6)  :: rtest, ctest, testr, testc

    integer, dimension(6)  :: tot_sitestPA
    integer, dimension(6)  :: prob_sitetPA
    integer, dimension(6)  :: sitestPA_N03
    integer, dimension(6)  :: sitestPA_N10
    integer, dimension(6)  :: sitestPA_Nplg
    integer, dimension(6)  :: sitestPA_Nplgn
    integer, dimension(6)  :: sitestPA_Nplgi
    integer, dimension(5)  :: tPANplgn

    integer, dimension(9)  :: newstate
    integer  :: statetype
    integer, dimension(1000000)  :: reaction
    double precision  :: randPLi, rmpt

    integer, dimension(:), allocatable  :: row2, col2
    integer, dimension(2)  :: rowtPA, coltPA
    double precision, dimension(:), allocatable  :: lysis_time
    double precision, dimension(:), allocatable  :: tPA_time
    integer, dimension(:), allocatable  :: tPAunbind, tPAPLiunbd, ltPA
    integer, dimension(:), allocatable  :: Plasmin
    integer, dimension(:), allocatable  :: max_Plg
    integer, dimension(:), allocatable  :: lysiscomplete
    integer, dimension(:), allocatable  :: countvect, PLGbd, PLGunbd
    double precision, dimension(100) :: PLitime
    integer  :: numN, numP1, numP1b, numP2, numP6, qconvert, colPLG, colPLGb, colPLGfind, ntest
    double precision, dimension(4) :: probnumP, cumprobnumP
    double precision :: numPt
    integer  :: n10, n12, n13, nplgi, nplgn2, n23, n33, n03, findit

    double precision :: pi, rP

    integer  :: Nsave, Ninteger
    double precision, dimension(96)  :: tsave, persave !96 because I save every percent from 0 to 95

    integer :: countfp
    double precision, dimension(:), allocatable :: firstPLi

    integer  :: PLGbind, PLGunbind

    character(40) :: filetype, formatted
    character(70) :: filename1
    logical       :: isBinary = .True.      ! flag for binary output
    integer       :: lysunit = 22
    integer       :: tPAunit = 23
    integer       :: PLiunit = 24
    integer       :: endunit = 25
    integer       :: lysunit2 = 26
    integer       :: tPAunit2 = 27
    integer       :: PLiunit2 = 28
    integer       :: endunit2 = 29
    integer       :: lysunit3 = 30
    integer       :: tPAunit3 = 31
    integer       :: PLiunit3 = 32
    integer       :: endunit3 = 33
    integer       :: lysunit4 = 34
    integer       :: tPAunit4 = 35
    integer       :: PLiunit4 = 36
    integer       :: endunit4 = 37
    integer       :: lysunit5 = 38
    integer       :: tPAunit5 = 39
    integer       :: PLiunit5 = 40
    integer       :: endunit5 = 41
    integer       :: lysunit6 = 42
    integer       :: tPAunit6 = 43
    integer       :: PLiunit6 = 44
    integer       :: endunit6 = 45
    integer       :: perdegunit = 20
    integer       :: tunit = 21
    integer       :: t2unit = 19
    integer       :: perunit = 18
    integer       :: ptunit = 17
    integer       :: countunit = 16
    integer       :: tvectunit = 46
    integer       :: sunit = 47
    integer       :: plgunit = 48
    integer       :: ctunit = 49
    integer       :: plgunbdunit = 50
    integer       :: plgbdunit = 51
    integer       :: plitimeunit = 52
    integer       :: tPAunbdunit = 53
    integer       :: tPAPLiunit = 54
    integer       :: s2unit = 55
    integer       :: fpunit = 56

    character(80) :: lysfile
    character(80) :: tPAfile
    character(80) :: PLifile
    character(80) :: lysfile2
    character(80) :: tPAfile2
    character(80) :: PLifile2
    character(80) :: lysfile3
    character(80) :: tPAfile3
    character(80) :: PLifile3
    character(80) :: lysfile4
    character(80) :: tPAfile4
    character(80) :: PLifile4
    character(80) :: lysfile5
    character(80) :: tPAfile5
    character(80) :: PLifile5
    character(80) :: endfile5
    character(80) :: lysfile6
    character(80) :: tPAfile6
    character(80) :: PLifile6
    character(80) :: endfile6
    character(80) :: perdegfile
    character(80) :: tfile
    character(80) :: endfile
    character(80) :: endfile2
    character(80) :: endfile3
    character(80) :: endfile4
    character(80) :: t2file
    character(80) :: perfile
    character(80) :: countfile
    character(80) :: ptfile
    character(80) :: tvectfile
    character(80) :: sfile
    character(80) :: plgfile
    character(80) :: ctfile
    character(80) :: plgbdfile
    character(80) :: plgunbdfile
    character(80) :: plitimefile
    character(80) :: tPAunbdfile
    character(80) :: tPAPLifile
    character(80) :: s2file
    character(80) :: fpfile

    !stuff for the random # generator. I need to use the file kiss.o when I compile in order for this to work, and kiss.o
    !is obtained by compiling the file kiss.c by doing "cc -c kiss.c".
    external :: mscw, kiss32, urcw1
    integer :: kiss32, mscw, stater(4), old_stater(4), ui
    double precision :: uf, urcw1

    !! BRAD 2024-01-13:

    integer :: cmd_count, param_i, param_len, param_val_len, cmd_status, io_status
    character(80) :: param_name, param_value

    cmd_count = command_argument_count()
    write (*, *) 'number of command arguments = ', cmd_count

    param_i = 0
    do while (param_i < cmd_count)
        param_i = param_i + 1
        call get_command_argument(param_i, param_name, param_len, cmd_status)
        if (cmd_status /= 0) then
            write (*, *) ' get_command_argument failed: status = ', cmd_status, ' arg = ', param_i
            stop
        end if
        write (*, *) 'command arg ', param_i, ' = ', param_name(1:param_len)
        param_i = param_i + 1
        call get_command_argument(param_i, param_value, param_val_len, cmd_status)
        if (cmd_status /= 0) then
            write (*, *) ' get_command_argument failed: status = ', cmd_status, ' arg = ', param_i
            stop
        end if
        write (*, *) 'command arg ', param_i, ' = ', param_value(1:param_val_len)
        select case (param_name(3:param_len))
        case ('expCode')
            expCode = param_value(1:param_val_len)
            write (*, *) 'Setting expCode = ', expCode
        case ('inFileCode')
            inFileCode = param_value(1:param_val_len)
            write (*, *) 'Setting inFileCode = ', inFileCode
        case ('outFileCode')
            outFileCode = param_value(1:param_val_len)
            write (*, *) 'Setting outFileCode = ', outFileCode
        case ('nodes')
            read (param_value, *, iostat=io_status) nodes
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting nodes = ', nodes
        case ('runs')
            read (param_value, *, iostat=io_status) runs
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting runs = ', runs
        case ('radius')
            read (param_value, *, iostat=io_status) radius
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting radius = ', radius
        case ('KdtPAyesplg')
            read (param_value, *, iostat=io_status) KdtPAyesplg
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting KdtPAyesplg = ', KdtPAyesplg
        case ('KdtPAnoplg')
            read (param_value, *, iostat=io_status) KdtPAnoplg
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting KdtPAnoplg = ', KdtPAnoplg
        case ('KdPLGintact')
            read (param_value, *, iostat=io_status) KdPLGintact
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting KdPLGintact = ', KdPLGintact
        case ('KdPLGnicked')
            read (param_value, *, iostat=io_status) KdPLGnicked
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting KdPLGnicked = ', KdPLGnicked
        case ('ktPAon')
            read (param_value, *, iostat=io_status) ktPAon
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting ktPAon = ', ktPAon
        case ('kplgon')
            read (param_value, *, iostat=io_status) kplgon
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting kplgon = ', kplgon
        case ('freeplg')
            read (param_value, *, iostat=io_status) freeplg
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting freeplg = ', freeplg
        case ('kdeg')
            read (param_value, *, iostat=io_status) kdeg
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting kdeg = ', kdeg
        case ('kplioff')
            read (param_value, *, iostat=io_status) kplioff
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting kplioff = ', kplioff
        case ('kapcat')
            read (param_value, *, iostat=io_status) kapcat
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting kapcat = ', kapcat
        case ('kncat')
            read (param_value, *, iostat=io_status) kncat
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting kncat = ', kncat
        case ('seed')
            read (param_value, *, iostat=io_status) seed
            if (io_status /= 0) then
                write (*, *) 'String conversion error'
                stop
            end if
            write (*, *) 'Setting seed = ', seed
        case default
            write (*, *) 'Unrecognized parameter'
            stop
        end select
    end do

    write (*, *) 'command line processed'

    allocate (state(Ntot, nodes**2))
    allocate (statetemp(Ntot, nodes**2))
    allocate (init_state(nodes))
    allocate (cumsuminit(nodes))
    allocate (Lat(nodes**2, nodes**2))
    allocate (row2(nodes**2*Ntot))
    allocate (col2(nodes**2*Ntot))
    allocate (lysis_time(runs))
    allocate (tPA_time(runs))
    allocate (tPAunbind(runs))
    allocate (tPAPLiunbd(runs))
    allocate (ltPA(runs))
    allocate (Plasmin(runs))
    allocate (max_Plg(runs))
    allocate (lysiscomplete(runs))
    allocate (countvect(runs))
    allocate (PLGbd(runs))
    allocate (PLGunbd(runs))
    allocate (firstPLi(runs))

    ! Calculate dependent parameters
    kplgoff = kplgon*KdPLGintact !units 1/s, off rate for intact FB
    kplgoffnick = kplgon*KdPLGnicked !units 1/(s*uM), off rate for nicked FB
    kaoff12 = ktPAon*KdtPAyesplg    !units 1/s, off rate in presence of PLG
    kaoff10 = ktPAon*KdtPAnoplg     !units 1/s, off rate in absense of PLG

    !! BRAD 2024-01-13:
    do i = 1, nodes
        do j = 1, nodes
            Lat((i - 1)*nodes + j, (i - 1)*nodes + j) = 1 ! Put a 1 along the diagonal for (i, j) to itself
            if (i /= 1) Lat((i - 1)*nodes + j, (i - 2)*nodes + j) = 1 ! UP ! Put a 1 from (i, j) to (i-1, j)
            if (j /= 1) Lat((i - 1)*nodes + j, (i - 1)*nodes + j - 1) = 1 ! RIGHT ! Put a 1 from (i, j) to (i, j-1)
            if (i /= nodes) Lat((i - 1)*nodes + j, (i)*nodes + j) = 1 ! DOWN ! Put a 1 from (i, j) to (i+1, j)
            if (j /= nodes) Lat((i - 1)*nodes + j, (i - 1)*nodes + j + 1) = 1 ! LEFT ! Put a 1 from (i, j) to (i+1, j)
        end do
    end do

    !! END BRAD

    if (isBinary) then
        !filetype = 'unformatted' !if you compile with gfortran or f95
        filetype = 'binary'      !if you compile with ifort
    else
        filetype = 'formatted'
    end if
    write (*, *) ' filetype=', filetype

    !the next part (up until the definition of pi) is related to the random numer generator
    ui = kiss32()

    uf = urcw1()

    !! BRAD 2024-01-13:
    if (seed == 0) seed = mscw()

    ! seed = mscw() !randomly generate seed
    ! seed=-34041038
    write (*, *), ' seed=', seed

    stater(1) = 129281
    stater(2) = 362436069
    stater(3) = 123456789
    stater(4) = seed
    call set_kiss32(stater)

    call get_kiss32(stater)
    !call vurcw1(rvect,M)

    Tdoublets = Ntot*nodes**2   !total number of doublets in system
    pi = ACOS(-1.0)
    ones = 1

    write (*, *), 'nodes=', nodes
    write (*, *), 'KdtPAnoplg=', KdtPAnoplg
    write (*, *), 'KdtPAyesplg=', KdtPAyesplg
    write (*, *), 'KdPLGnicked=', KdPLGnicked
    write (*, *), 'KdPLGintact=', KdPLGintact

    !parameters we use to find reaction rates

    !columns of param array represent, in order: kdeg, freeplg, kplgoff, kplgon, kplgoffnick, kplioff, kaoff12, kaoff10,
    !ktPAon, kapcat, kncat
    !case 1 (baseline parameter value)
    !param(1,1) = 5                     !kdeg    units 1/s, plasmin-mediated degradation rate
    !param(2,1) = freeplg                   !freeplg  constant concentration of free plasminogen, units uM
    !param(3,1) = kplgon*KdPLGintact    !kplgoff  units 1/s, off rate for intact FB
    !param(4,1) = kplgon                !kplgon   units 1/(s*uM)
    !param(5,1) = kplgon*KdPLGnicked    !kplgoffnick  units 1/(s*uM), off rate for nicked FB
    !param(6,1) = 57.6                  !kplioff  units 1/s, PLi unbinding rate
    !param(7,1) = ktPAon*KdtPAyesplg    !kaoff12  units 1/s, off rate in presence of PLG
    !param(8,1) = ktPAon*KdtPAnoplg     !kaoff10  units 1/s, off rate in absense of PLG
    !param(9,1) = ktPAon                !ktPAon   units 1/(s*uM)
    !param(10,1) = 0.1                  !kapcat  units 1/s, tPA-mediated rate of conversion of PLG to PLi
    !param(11,1) = 5                    !kncat  units 1/s, PLi-mediated rate of exposure of new binding sites

    !case 2 (baseline parameter values with [plg]=0.92 uM).
    !param(1,2) = 5     !kdeg    units 1/s, degradation rate
    !param(2,2) = 0.92      !freeplg  constant concentration of free plasminogen, in muM
    !param(3,2) = 3.8    !kplgoff  units 1/s
    !param(4,2) = 0.1    !kplgon   units 1/(s*muM)
    !param(5,2) = 1.7273 !kplgonnick  units 1/(s*muM), so Kd for nicked FB is 2.2 muM
    !param(6,2) = 57.6     !kplioff  units 1/s
    !param(7,2) = 0.0002 !kaoff12  units 1/s
    !param(8,2) = 0.0036 !kaoff10  units 1/s
    !param(9,2) = 0.01   !ktPAon   units 1/(s*muM)
    !param(10,2) = 0.1    !kapcat  units 1/s
    !param(11,2) = 5     !kncat  units 1/s

    !stoichiometric matrix. each row corresponds to a different reaction and each column represents a different state:
    ![Nplg,N12,N10,N13,N03,N23,N33,degradation,degradation plus 1 PLi]

    V_tmp = (/(/0, -1, 1, 0, 0, 0, 0, 0, 0/), &  !, N12 goes to  N10
              (/0, 1, -1, 0, 0, 0, 0, 0, 0/), &  !N10 goes to  N12
              (/0, 0, 1, -1, 0, 0, 0, 0, 0/), &  !N13 goes to  N10
              (/0, 0, 0, -1, 1, 0, 0, 0, 0/), &  !N13 goes to  N03
              (/0, 0, 0, 0, -1, 1, 0, 0, 0/), &  !N03 goes to  N23
              (/0, 0, 0, 0, 1, -1, 0, 0, 0/), &  !N23 goes to  N03
              (/0, 0, 0, 0, 1, 0, -1, 0, 0/), &  !N33 goes to  N03
              (/1, -1, 0, 0, 0, 0, 0, 0, 0/), &  !N12 goes to  Nplg
              (/1, 0, -1, 0, 0, 0, 0, 0, 0/), &  !N10 goes to  Nplg
              (/1, 0, 0, 0, -1, 0, 0, 0, 0/), &  !N03 goes to  Nplg
              (/1, 0, 0, 0, 0, -1, 0, 0, 0/), &  !N23 goes to  Nplg
              !(/0, -1, 0, 1, 0, 0, 0, 0, 0/),&  !tPA converts PLG to PLi only on the same doublet. i.e. N12 goes to  N13
              (/0, 0, 0, 0, 0, 0, 0, 0, 0/), &  !tPA converts PLG to PLi on any doublet at same binding location
              (/0, 0, 0, 0, 0, 0, 0, 1, -1/), &  !degraded doublet with one PLi bound becomes a degraded doublet with no PLi bound
              (/0, 0, 0, -1, 0, 0, 0, 1, 0/), &  !N13 degrades
              (/0, 0, 0, 0, -1, 0, 0, 1, 0/), &  !N03 degrades
              (/0, 0, 0, 0, 0, -1, 0, 1, 0/), &  !N23 degrades
              (/0, 0, 0, 0, 0, 0, -1, 1, 0/), &  !N33 degrades
              (/0, 0, 0, 0, 0, 0, 0, 0, 0/), &   !N13 converts N to Nplg
              (/0, 0, 0, 0, 0, 0, 0, 0, 0/), &   !N03 converts N to Nplg
              (/0, 0, 0, 0, 0, 0, 0, 0, 0/), &   !N23 converts N to Nplg
              (/0, 0, 0, 0, 0, 0, 0, 0, 0/), &   !N33 converts N to Nplg
              (/0, 0, 0, 0, 0, 0, 0, 0, 0/), &   !degraded with PLi converts N to Nplg
              (/0, 0, -1, 0, 0, 0, 0, 1, 0/), &  !N10 degrades
              (/0, -1, 0, 0, 0, 0, 0, 1, 0/), &  !N12 degrades
              (/-1, 0, 0, 0, 0, 0, 0, 1, 0/), &  !intact Nplg degrades
              (/-1, 0, 0, 0, 0, 0, 0, 1, 0/)/)  !nicked Nplg degrades

    V = reshape(V_tmp, (/9, sV/))

    q = 0 !q is a diagonal matrix with 1's down the diagonal
    do i = 1, 9
        q(i, i) = 1
    end do

    lysiscomplete = 0
    firstPLi = 0
    countfp = 0

    write (*, *) 'runs=', runs

    !below starts the big loop where we do 50000 independent simulations ("do stats=1,runs")
    do stats = 1, runs
        if (MODULO(stats, 1000) == 0) write (*, *) ' stats=', stats

        Nsave = 0
        persave = 0
        tsave = 0

        if (stats <= runs) then   !change back to 10000 when I change runs
            !ipar=1 !says if stats<=50000, use the "case 1" parameters defined around line 357

            !elseif(stats<=20000.and.stats.gt.10000) then
            !ipar=16
        end if

        if (stats == 1) then
            write (lysfile, '(58a)') 'data/'//expCode//'/lysis_'//outFileCode
            write (tPAfile, '(61a)') 'data/'//expCode//'/tPA_time_'//outFileCode
            write (PLifile, '(56a)') 'data/'//expCode//'/PLi_'//outFileCode
            write (endfile, '(64a)') 'data/'//expCode//'/lyscomplete_'//outFileCode
            !write(plgfile,'(56a)' ) 'data/' // expCode // '/PLG_' // outFileCode
            !write(ctfile,'(58a)' ) 'data/' // expCode // '/count_' // outFileCode
            !write(plgbdfile,'(27a)') 'data/' // expCode // '/PLGunbindPLG2_tPA01_Q2.dat'
            !write(plgunbdfile,'(29a)') 'data/' // expCode // '/PLGbindPLG2_tPA01_Q2.dat'
            !write(plitimefile,'(63a)') 'data/' // expCode // '/PLitime_' // outFileCode
            write (tPAPLifile, '(63a)') 'data/'//expCode//'/tPAPLiunbd_'//outFileCode
            !write(sfile,'(28a)' ) 'data/' // expCode // '/statetPAPLG2_tPA01_Q2.dat'
            !write(profile,'(23a)' ) 'data/' // expCode // '/ersavePLG2_tPA01_Q2.dat'
            !write(t2file,'(23a)' ) 'data/' // expCode // '/tsavePLG2_tPA01_Q2.dat'
            write (tPAunbdfile, '(61a)') 'data/'//expCode//'/tPAunbind_'//outFileCode
            write (s2file, '(67a)') 'data/'//expCode//'/lasttPA_'//outFileCode
            write (fpfile, '(77a)') 'data/'//expCode//'/firstPLi_'//outFileCode
            open (lysunit, file=lysfile, form=filetype)
            open (tPAunit, file=tPAfile, form=filetype)
            open (PLiunit, file=PLifile, form=filetype)
            open (endunit, file=endfile, form=filetype)
            !open(plgunit,file=plgfile,form=filetype)
            !open(ctunit,file=ctfile,form=filetype)
            !open(plgbdunit,file=plgbdfile,form=filetype)
            !open(plgunbdunit,file=plgunbdfile,form=filetype)
            !open(plitimeunit,file=plitimefile,form=filetype)
            open (tPAPLiunit, file=tPAPLifile, form=filetype)
            !open(sunlit,file=sfile,form=filetype)
            !open(permit,file=perfile,form=filetype)
            !open(t2unit,file=t2file,form=filetype)
            open (tPAunbdunit, file=tPAunbdfile, form=filetype)
            open (s2unit, file=s2file, form=filetype)
            open (fpunit, file=fpfile, form=filetype)

            write (*, *) ' kncat=', kncat
            write (*, *) ' kapcat=', kapcat
            write (*, *) ' ktPAon=', ktPAon
            write (*, *) ' kaoff10=', kaoff10
            write (*, *) ' kaoff12=', kaoff12
            write (*, *) ' kplioff=', kplioff
            write (*, *) ' kplgoffnick=', kplgoffnick
            write (*, *) ' kplgon=', kplgon
            write (*, *) ' kplgoff=', kplgoff
            write (*, *) ' freeplg=', freeplg
            write (*, *) ' kdeg=', kdeg

        end if

        ! set the reaction rates to be the appropriate parameters.
        ! I think we likely don't need this step in future versions of the code;
        ! we could directly use the parameters,
        ! e.g. "kplgon*KdPLGintact" to represent the unbinding rate of PLG
        c2 = kplgoff
        c3 = kplgon*freeplg
        c5 = kplioff
        c6 = kaoff12
        c8 = kplgon*freeplg
        c8a = kplgoffnick!*freeplg
        c9 = kplgoff
        c11 = kplioff
        c12 = kaoff12
        c13 = kaoff10
        c14 = kplioff
        c15 = kplioff
        c20 = kapcat
        c25 = kplioff
        ! kdeg = param(1,ipar)
        ! kncat = param(11,ipar)

        !define the plg probabilities (calculated in thesis, or see "PLG_prob_calc.mw").
        !These are the probabilities of an Nplg doublet being N02, N00, or N22. i.e., having 1 PLG bound, no PLG bound, or 2 PLG bound
        !probability for fully intact fibrin
        prob_N02i = 1/(1 + (0.5*KdPLGintact/freeplg) + (0.5*freeplg/KdPLGintact))
        prob_N00i = 0.5*KdPLGintact/(freeplg*(1 + (0.5*KdPLGintact/freeplg) + (0.5*freeplg/KdPLGintact)))
        prob_N22i = 0.5*freeplg/(KdPLGintact*(1 + (0.5*KdPLGintact/freeplg) + (0.5*freeplg/KdPLGintact)))

        !probability for nicked fibrin
        prob_N02n = 1/(1 + (0.5*KdPLGnicked/freeplg) + (0.5*freeplg/KdPLGnicked))
        prob_N00n = 0.5*KdPLGnicked/(freeplg*(1 + (0.5*KdPLGnicked/freeplg) + (0.5*freeplg/KdPLGnicked)))
        prob_N22n = 0.5*freeplg/(KdPLGnicked*(1 + (0.5*KdPLGnicked/freeplg) + (0.5*freeplg/KdPLGnicked)))

        !write(*,*)'prob_N02i=',prob_N02i
        !write(*,*)'prob_N00i=',prob_N00i
        !write(*,*)'prob_N22i=',prob_N22i
        !write(*,*)'prob_N02n=',prob_N02n
        !write(*,*)'prob_N00n=',prob_N00n
        !write(*,*)'prob_N22n=',prob_N22n
        !write(*,*)'sumNs=',prob_N02i+prob_N00i+prob_N22i+prob_N02n+prob_N22n+prob_N00n

        state = 0
        state(1:Nplginit, 1:nodes**2) = 1 !set the initial state of all exposed doublets to be Nplg (1)

        !begin by setting one of the exposed doublets to an activated doublet: either N12 (2) or N10 (3).

        init_state = 1.0/nodes   !vector that's length of one row of lattice, scaled so sum(init_state)=1

        rin = urcw1()  !pick a random number and use that to decide where we place the 1st tPA

        call cumsum(init_state, nodes, cumsuminit)
        do i = 1, nodes
            if (cumsuminit(i) >= rin) then
                init_entry = i !find the first entry in the cumsum vector that is >= random number, and have this be the doublet tPA starts on
                exit
            end if
        end do

        !! BRAD   write(*,*)'init_entry=',init_entry

        tPA_loc(1) = init_entry; !tPA starts on node init_entry

        r = urcw1()  ! pick a different random number to determine whether we get N12 or N10 I don't take into account the fact
        ! that the Nplg on that site could be an N22...is this a problem? What I do below just says, "if r isn't
        ! less than prob_N02, make state be a 3" even though there's a small chance it is an N22, and therefore tPA
        ! can't bind...I don't think it's an issue, because in the macroscale model I'm assuming tPA DOES bind,
        ! so it must have a site to bind to (i.e., not be N22)
        if (r <= prob_N02i) then
            state(1, init_entry) = 2 !set the state of the doublet to be N12
        else
            state(1, init_entry) = 3 !set the state of the double to be N10
        end if

        !! BRAD  write(*,*)'state(1,init_entry)=',state(1,init_entry)

        t = 0
        tfinal = 24*60*60   ! final time to run Gillespie algorithm for
        ! (if the algorithm doesn't stop before then because fiber degraded or all tPA and PLi unbound)
        count = 1
        countPLi = 0
        PLGbind = 0
        PLGunbind = 0
        plotnumber = 1
        PLitime = 0
        reaction2 = 0
        reaction3 = 0
        tvals = 0
        reaction = 0
        Anumber = 0
        PLi = 0
        PLideg = 0
        PLilocation = 0
        PLiwaste = 0
        tPA = 0
        tPA_loc = 0
        PLG = 0
        per_lysis = 0
        tau2 = 0
        percent_degrade = 0
        timevect = 0
        Prob_rebind = 0
        degrade_tot = 0
        undegrade_tot = 0
        undegrade_new = 0
        undegrade_tot(1) = Nplginit*nodes**2
        undegrade_new(1) = 0
        PLi_loc = 0
        Anumber(1) = 1       !we always start with just one active site
        temp = 0
        PLi(1) = 0         !we initially start with no PLi in the system
        tPA(1) = 1         !we initially start with 1 tPA

        countstate = 0
        countstate2 = 0

        !find the number of doublets in states 2 and 6, i.e. the number of N12 and N23 doublets. the goal is to count up the amount of PLG
        do i = 1, nodes**2
            do j = 1, Ntot
                if (state(j, i) == 2 .or. state(j, i) == 6) then
                    countstate = countstate + 1
                end if
            end do
        end do

        !find the number of doublets in state 1, i.e. the number of Nplg doublets.
        do i = 1, nodes**2
            do j = 1, Ntot
                if (state(j, i) == 1) then
                    countstate2 = countstate2 + 1
                end if
            end do
        end do

        PLG(1) = countstate + (prob_N02i + 2*prob_N22i)*countstate2 ! the initial amount of PLG
        ! = (amount in states N12 and N23)
        ! + (amount in state Nplg)*(probability the doublet was N02 or N22)
        degrade_tot(1) = 0          !we initially start with no degraded doublets
        !  num_degrade=0
        !  num_undegrade=0 !don't think this gets used anymore
        !  num_undegrade(1,1:nodes**2)=Nplginit !don't think this gets used anymore
        percent_degrade(1) = 0 !vector that saves the percentage of doublets that have been degraded at each time step

        D = 10**(7)             !diffusion coefficient in units of nm^2/s

        vol3 = 0.001*(radius)**2*pi   !volume of 1 "location", i.e. cleavage site, is 1nm=0.001microns times radius squared times pi.

        bs = Nplginit*nodes**2/vol3  !# of binding sites per volume where we take only the volume of the current cleavage site -
        ! not the whole fiber.
        bs = bs/602.2                !converts bs from units of #/volume to units of micromolar

        !Now do the Gillespie algorithm
        do

            if (t > tfinal) exit     !only do the loop while t<=tfinal

            !find where "active" sites (i.e., doublets with tPA and/or plasmin) are located in state vector
            countstate3 = 0
            row = 0
            col = 0

            do j = 1, Ntot
                do i = 1, nodes**2
                    if (state(j, i) > 1) then
                        countstate3 = countstate3 + 1   ! if the state of the doublet is larger than 1 (i.e., any state other than Nplg),
                        ! count it as an active site
                        row(countstate3) = i
                        col(countstate3) = j
                    end if
                end do
            end do
            A = countstate3   !number of activated states

            if (A == 0) then    !if there are no activated sites, stop the algorithm
            !! BRAD         write(*,*)'no tPA - stopped algorithm'
                exit
            end if

            tau = 0       !reset tau to zero each time
            newsite = 0
            newsite2 = 0
            newdoublet = 0
            newdoublet2 = 0
            P = 0
            currentPLi = 0
            countcurPLi = 0

            !loop over all active sites. Every doublet that has tPA or plasmin bound is a doublet on which the next reaction can happen. We'll compile all possible "next" reactions and corresponding reaction times, and then choose the one that happens soonest
            do i = 1, A

                qstate(i, :) = q(state(col(i), row(i)), :) !qstate(i,:) is a vector that has a "1" in the entry corresponding to the state of the active site, and 0's elsewhere
                qconvert = 0
                numN = 0
                numP2 = 0
                numP1 = 0
                numP1b = 0
                numPt = 0
                numP6 = 0
                n12 = 0
                n10 = 0
                n13 = 0
                n03 = 0
                n23 = 0
                n33 = 0
                nplgi = 0
                nplgn2 = 0

                call findint(state(:, row(i)), Ntot, 0, numN) !find the number of N's (i.e., cryptic, unexposed doublets) on the binding location of interest

                !!!uncomment below if tPA can convert any PLG at same binding location
                !if the doublet of interest has a tPA bound (i.e., state is 2, 3, or 4), find the number of PLG molecules on the current binding location
                if (state(col(i), row(i)) >= 2 .and. state(col(i), row(i)) <= 4) then
                    call findint(state(:, row(i)), Ntot, 2, numP2) !find number of N12
                    call findint(state(2:Ntot, row(i)), Ntot - 1, 1, numP1) !find number of exposed (or "nicked") Nplg
                    call findint(state(1, row(i)), 1, 1, numP1b) !find number of intact Nplg
                    call findint(state(:, row(i)), Ntot, 6, numP6) !find number of N23
                    numPt = (2*prob_N22n + prob_N02n)*numP1 + numP2 + numP6 + (2*prob_N22i + prob_N02i)*numP1b
                    if (numPt > 0) then  !if there are PLG molecules at the current binding location AND a tPA molecule
                        qconvert = 1
                    end if
                end if
                !!!end commentable section

                !if the doublet of interest has a PLi bound, find the number of exposed doublets of each type at the current binding location
                if (state(col(i), row(i)) >= 4 .and. state(col(i), row(i)) <= 7 .or. state(col(i), row(i)) == 9) then
                    call findint(state(:, row(i)), Ntot, 2, n12) !find number of N12
                    call findint(state(:, row(i)), Ntot, 3, n10) !find number of N10
                    call findint(state(:, row(i)), Ntot, 4, n13) !find number of N13
                    call findint(state(:, row(i)), Ntot, 5, n03) !find number of N03
                    call findint(state(:, row(i)), Ntot, 6, n23) !find number of N23
                    call findint(state(:, row(i)), Ntot, 7, n33) !find number of N33
                    call findint(state(1, row(i)), 1, 1, nplgi)  !find number of intact Nplg
                    call findint(state(2:Ntot, row(i)), Ntot - 1, 1, nplgn2) !find number of nicked Nplg
                end if

                !q(state(col(i),row(i)),:) is the column of the identity matrix corresponding to the
                !current state of active doublet i. aa is a matrix of propensity fcns.
                !the first entry of aa represents the doublet. So if there are 3 "active" doublets (i.e., doublets with tPA or plasmin bound), then aa(1:3,:). If there's only 1, we have aa(1,:). The second entry of aa represents the reaction number. So in this case, there are 26 different reactions that can occur.

                !WHAT IS qstate(i,:)?!?!?! It's a vector with 0's in all entries except for the entry corresponding to the state of the doublet (which has entry "1"). 2nd entry in qstate represents the type of doublet we have: 1=Nplg, 2=N12, 3=N10, 4=N13, 5=N03, 6=N23, 7=N33, 8=degraded, 9=degraded with 1 plasmin bound. so qstate(i,4) means that the ith active site is of type 4 (N13)

                if (col(i) == 1) then       !if the doublet is in column 1 of state vector (i.e. was initially exposed) use intact fibrin rates
                    aa(i, 1) = c2*qstate(i, 2)   !N12 goes to N10 on intact
                    aa(i, 6) = c2*qstate(i, 6)   !N23 goes to N03 on intact
                elseif (col(i) > 1) then       !if doublet is in columns 2-6 of state vector (i.e. was unexposed) use nicked fibrin rates
                    aa(i, 1) = c8a*qstate(i, 2)   !N12 goes to N10 on nicked
                    aa(i, 6) = c8a*qstate(i, 6)   !N23 goes to N03 on nicked
                end if

                aa(i, 2) = c3*qstate(i, 3)     !N10 goes to N12
                aa(i, 5) = c3*qstate(i, 5)     !N03 goes to N23
                aa(i, 3) = c5*qstate(i, 4)   !N13 goes to N10
                aa(i, 4) = c6*qstate(i, 4)   !N13 goes to N03
                aa(i, 7) = c11*qstate(i, 7)   !N33 goes to N03
                aa(i, 8) = c12*qstate(i, 2)   !N12 goes to Nplg
                aa(i, 9) = c13*qstate(i, 3)   !N10 goes to Nplg
                aa(i, 10) = c14*qstate(i, 5)   !N03 goes to Nplg
                aa(i, 11) = c15*qstate(i, 6)   !N23 goes to Nplg
                !aa(i,12)=c20*qstate(i,2)   !N12 goes to N13
                !comment out below aa(i,12) when tPA only converts PLG on same doublet. Comment out above aa(i,12) when it can convert any PLG on same location
                aa(i, 12) = c20*numPt*qconvert   !tPA converts PLG to PLi
                aa(i, 13) = c25*qstate(i, 9)   !degraded doublet with one PLi goes to degraded doublet w/no PLi
                aa(i, 14) = kdeg*n13    !N13 goes to degraded
                aa(i, 15) = kdeg*n03    !N03 goes to degraded
                aa(i, 16) = kdeg*n23    !N23 goes to degraded
                aa(i, 17) = kdeg*n33    !N33 goes to degraded
                aa(i, 18) = kncat*numN*qstate(i, 4) !N13 converts N to Nplg
                aa(i, 19) = kncat*numN*qstate(i, 5) !N03 converts N to Nplg
                aa(i, 20) = kncat*numN*qstate(i, 6) !N23 converts N to Nplg
                aa(i, 21) = kncat*numN*qstate(i, 7) !N33 converts N to Nplg
                aa(i, 22) = kncat*numN*qstate(i, 9) !degraded with PLi converts N to Nplg
                aa(i, 23) = kdeg*n10    !N10 goes to degraded
                aa(i, 24) = kdeg*n12    !N12 goes to degraded
                aa(i, 25) = kdeg*nplgi  !intact Nplg goes to degraded
                aa(i, 26) = kdeg*nplgn2  !nicked Nplg goes to degraded

                !find the sum of the propensity functions
                asum(i) = sum(aa(i, :))
                asumtemp = aa(i, :)/asum(i) !create new vector of propensities that is scaled by the total sum
                cumsumaa = 0
                call cumsum(asumtemp, sV, cumsumaa) !define cumsumaa to be a vector of the cumulative sums of the asumtemp vector. e.g., if asumtemp=(0.3, 0.6, 0.1) then cumsumaa=(0.3,0.9,1.0)
                aa(i, :) = 0   !reset to 0

                r = urcw1() !pick a random number to choose next time reaction occurs
                rin = urcw1() !pick another random number to choose next reaction that occurs
                tau(i) = log(1/r)/asum(i)       !picks the next time an event will occur
                asum(i) = 0
                do j = 1, sV !loop over all possible reactions and choose the first one for which cumsumaa(j) is bigger than or equal to rin
                    if (cumsumaa(j) >= rin) then
                        bigJ(i) = j          !picks the next state by choosing which reaction happens next. j should be a number between 1 and 26
                        exit
                    end if
                end do
            end do

            !choose the smallest time (out of all possible reaction times found above), and make that the next time step. then change the state associated with that time.
            !"place" is the location in the vector the smallest entry is found, and "entry" is its value

            call minvect(tau, 12000, place)       !find the location of minimum non-zero value in tautot vector
            entry = tau(place)

            !reset the unsused qstates and J's to 0. Reset the used qstate and bigJ to 0 after using them below
            do i = 1, A
                if (i /= place) then
                    qstate(i, :) = 0
                    bigJ(i) = 0
                end if
            end do

            tot_sitesPLi = 0
            prob_sitePLi = 0
            sites_N03 = 0
            sites_N10 = 0
            sites_Nplg = 0
            sites_Nplgn = 0
            sites_Nplgi = 0
            sites_degrade = 0

            alpha = 1    !parameter describing how likely PLi is to move onto a degraded doublet rather than an
            !undegraded one. 0<alpha<=1. alpha=1 means PLi is equally likely to move onto either

            !change the state of the doublet to the new one determined by rxn J
            !proceed as normal with the Gillespie steps

            if (bigJ(place) == 14) then
                !if N13 degrades
                call findfirst(state(:, row(place)), Ntot, 4, findit) !find the first entry that is N13
                if (findit == 0) write (*, *) 'Problem: should have an N13'
                state(findit, row(place)) = -1 !set the state of the doublet to be -1, i.e., degraded

                !I think the next do and if statements are mostly checks that nothing funny is going on. The critical step is defining "loc", the location of the plasmin molecule
                do iPLi = 1, countPLi !for all the plasmin in the system, check to see how much of it is involved in the current reaction
                    if (PLilocation(1, iPLi) == row(place) .and. PLilocation(2, iPLi) == findit) then
                        countcurPLi = countcurPLi + 1 !counts the number of plasmin molecules on the current binding location
                        currentPLib(countcurPLi) = iPLi
                    end if
                end do

                if (countcurPLi == 2) then !finds location of PLi involved in this event
                    write (*, *) 'Problem: should only have one plasmin'
                    write (*, *) 'bigJ(place)=', bigJ(place)
                    write (*, *) 'countcurPLi=', countcurPLi
                elseif (countcurPLi == 1) then
                    loc = currentPLib(1)
                else
                    write (*, *) 'Problem: should be 1 plasmin'
                    write (*, *) 'bigJ(place)=', bigJ(place)
                    write (*, *) 'countcurPLi=', countcurPLi
                end if

                count = count + 1 !update "count", which keeps track of the number of time steps taken
                t = t + entry        !entry=tau(place);

                tvals(count) = t !saves the time associated with time step "count"
                Anumber(count) = A !saves the number of active sites at time step "count"
                reaction(count) = bigJ(place) !saves the reaction that happened at time step "count"
                connect_ind = 0
                countlat = 0
                do j = 1, nodes**2
                    if (Lat(j, PLilocation(1, loc)) > 0) then    !find nodes neighboring current node
                        countlat = countlat + 1
                        connect_ind(countlat) = j
                    end if

                end do
                sizenonzero_connectind = countlat

                !release one tPA and one PLi (because when doublet degrades, tPA and plasmin are released). start with tPA:

                call movetpa(Ntot, nodes, state, prob_N02i, prob_N00i, prob_N02n, prob_N00n, vol3, ktPAon, D, count, p_rebind)
                !! BRAD          write(*,*)'p_rebind=',p_rebind

                reaction3 = reaction3 + 1

                !then move PLi and update the state of the doublets
                call movepli(Ntot, nodes, sizenonzero_connectind, state, connect_ind, alpha, prob_N02i, prob_N00i, &
                             prob_N02n, prob_N00n, loc, PLilocation, PLilocationt, statetemp)
                PLilocation = PLilocationt
                state = statetemp
            elseif (bigJ(place) == 15) then
                !if N03 degrades
                call findfirst(state(:, row(place)), Ntot, 5, findit) !find the first entry that is N03
                if (findit == 0) write (*, *) 'Problem: should have an N03'
                state(findit, row(place)) = -1

                do iPLi = 1, countPLi
                    if (PLilocation(1, iPLi) == row(place) .and. PLilocation(2, iPLi) == findit) then
                        countcurPLi = countcurPLi + 1
                        currentPLib(countcurPLi) = iPLi
                    end if
                end do

                if (countcurPLi == 2) then !finds location of PLi involved in this event
                    write (*, *) 'Problem: should only have one plasmin'
                    write (*, *) 'bigJ(place)=', bigJ(place)
                    write (*, *) 'countcurPLi=', countcurPLi
                elseif (countcurPLi == 1) then
                    loc = currentPLib(1)
                else
                    write (*, *) 'Problem: should be 1 plasmin'
                    write (*, *) 'bigJ(place)=', bigJ(place)
                    write (*, *) 'countcurPLi=', countcurPLi
                end if

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)
                connect_ind = 0
                countlat = 0
                do j = 1, nodes**2
                    if (Lat(j, PLilocation(1, loc)) > 0) then    !find nodes neighboring current node
                        countlat = countlat + 1
                        connect_ind(countlat) = j
                    end if

                end do
                sizenonzero_connectind = countlat

                !release 1 PLi and update state of the doublets

                call movepli(Ntot, nodes, sizenonzero_connectind, state, connect_ind, alpha, prob_N02i, prob_N00i, &
                             prob_N02n, prob_N00n, loc, PLilocation, PLilocationt, statetemp)
                PLilocation = PLilocationt
                state = statetemp

            elseif (bigJ(place) == 16) then
                !if N23 degrades
                call findfirst(state(:, row(place)), Ntot, 6, findit) !find the first entry that is N23
                if (findit == 0) write (*, *) 'Problem: should have an N23'
                state(findit, row(place)) = -1

                do iPLi = 1, countPLi
                    if (PLilocation(1, iPLi) == row(place) .and. PLilocation(2, iPLi) == findit) then
                        countcurPLi = countcurPLi + 1
                        currentPLib(countcurPLi) = iPLi
                    end if
                end do

                if (countcurPLi == 2) then !finds location of PLi involved in this event
                    write (*, *) 'Problem: should only have one plasmin'
                    write (*, *) 'bigJ(place)=', bigJ(place)
                    write (*, *) 'countcurPLi=', countcurPLi
                elseif (countcurPLi == 1) then
                    loc = currentPLib(1)
                else
                    write (*, *) 'Problem: should be 1 plasmin'
                    write (*, *) 'bigJ(place)=', bigJ(place)
                    write (*, *) 'countcurPLi=', countcurPLi
                end if

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)
                connect_ind = 0
                countlat = 0

                do j = 1, nodes**2
                    if (Lat(j, PLilocation(1, loc)) > 0) then    !find nodes neighboring current node
                        countlat = countlat + 1
                        connect_ind(countlat) = j
                    end if

                end do
                sizenonzero_connectind = countlat

                !release 1 PLi and update state of the doublets

                call movepli(Ntot, nodes, sizenonzero_connectind, state, connect_ind, alpha, prob_N02i, prob_N00i, &
                             prob_N02n, prob_N00n, loc, PLilocation, PLilocationt, statetemp)
                PLilocation = PLilocationt
                state = statetemp

            elseif (bigJ(place) == 17) then

                !if N33 degrades
                call findfirst(state(:, row(place)), Ntot, 7, findit) !find the first entry that is N33
                if (findit == 0) write (*, *) 'Problem: should have an N33'
                state(findit, row(place)) = -1

                do iPLi = 1, countPLi
                    if (PLilocation(1, iPLi) == row(place) .and. PLilocation(2, iPLi) == findit) then
                        countcurPLi = countcurPLi + 1
                        currentPLib(countcurPLi) = iPLi
                    end if
                end do

                if (countcurPLi == 2) then !finds location of PLi involved in this event
                    loc = currentPLib(1)
                    loc2 = currentPLib(2)
                elseif (countcurPLi == 1) then
                    write (*, *) 'Problem: should have two plasmins'
                else
                    write (*, *) 'Problem: should be 2 plasmins'
                end if

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)
                connect_ind = 0
                countlat = 0

                do j = 1, nodes**2
                    if (Lat(j, PLilocation(1, loc)) > 0) then    !find nodes neighboring current node
                        countlat = countlat + 1
                        connect_ind(countlat) = j
                    end if

                end do
                sizenonzero_connectind = countlat

                !release 2 PLi and update state
                PLilocsmall(1) = PLilocation(1, loc)
                PLilocsmall(2) = PLilocation(2, loc)
                do iPLi2 = 1, countPLi
                  if (PLilocation(1, iPLi2) == PLilocsmall(1) .and. PLilocation(2, iPLi2) == PLilocsmall(2) .and. iPLi2 /= loc) then
                        loc2b = iPLi2
                        if (loc2b /= loc2) then
                            write (*, *) 'Problem: loc2 and loc2b not the same'
                            write (*, *) 'loc2=', loc2
                            write (*, *) 'loc2b=', loc2b
                        end if
                    end if
                end do

                call move2pli(Ntot, nodes, sizenonzero_connectind, state, connect_ind, alpha, prob_N02i, prob_N00i, prob_N02n, &
                              prob_N00n, loc, loc2, PLilocation, PLilocationt, statetemp)

                PLilocation = PLilocationt
                state = statetemp

            elseif (bigJ(place) == 23) then

                !if N10 degrades

                call findfirst(state(:, row(place)), Ntot, 3, findit) !find the first entry that is N10
                if (findit == 0) write (*, *) 'Problem: should have an N10'
                state(findit, row(place)) = -1

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

                !release tPA

                call movetpa(Ntot, nodes, state, prob_N02i, prob_N00i, prob_N02n, prob_N00n, vol3, ktPAon, D, count, p_rebind)
                !! BRAD          write(*,*)'p_rebind=',p_rebind

                reaction3 = reaction3 + 1

            elseif (bigJ(place) == 24) then

                !if N12 degrades

                call findfirst(state(:, row(place)), Ntot, 2, findit) !find the first entry that is N12
                if (findit == 0) write (*, *) 'Problem: should have an N12'
                state(findit, row(place)) = -1

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

                !release tPA

                call movetpa(Ntot, nodes, state, prob_N02i, prob_N00i, prob_N02n, prob_N00n, vol3, ktPAon, D, count, p_rebind)
                !! BRAD          write(*,*)'p_rebind=',p_rebind

                reaction3 = reaction3 + 1

            elseif (bigJ(place) == 25) then

                !if intact Nplg degrades

                if (state(1, row(place)) == 1) then
                    state(1, row(place)) = -1
                else
                    write (*, *) 'Problem: should have an intact Nplg'
                end if

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

            elseif (bigJ(place) == 26) then

                !if nicked Nplg degrades

                call findfirst(state(2:Ntot, row(place)), Ntot - 1, 1, findit) !find the first entry that is nicked Nplg
                if (findit == 0) write (*, *) 'Problem: should have a nicked Nplg'
                findit = findit + 1  !because we searched entries 2-6, not 1-6.
                state(findit, row(place)) = -1

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

         elseif (bigJ(place) == 18 .or. bigJ(place) == 19 .or. bigJ(place) == 20 .or. bigJ(place) == 21 .or. bigJ(place) == 22) then
                !if exposure of binding site occurs:

                !update the current doublet to its new state
                newstate = qstate(place, :) + V(:, bigJ(place))
                qstate(place, :) = 0   !reset qstate to 0 so we're ready for the next timestep

                do j = 1, 9
                    if (newstate(j) > 0) then
                        statetype = j
                        exit
                    end if
                end do

                if (statetype == state(col(place), row(place))) then       !if the rxn that occurs is one in which N goes to Nplg, then
                    call findfirst(state(:, row(place)), Ntot, 0, NtoNplg)  !find the next 0 entry in row of state matrix
                    state(NtoNplg, row(place)) = 1                         !set the 0 entry to 1, i.e. change N to Nplg
                    if (NtoNplg == 0) write (*, *) 'Problem: there must be Ns to expose if this rxn occurs'
                else
                    write (*, *) 'Problem:this should be exposure so statetype should equal state(col(place),row(place))'
                end if

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

           elseif (bigJ(place) == 3 .or. bigJ(place) == 7 .or. bigJ(place) == 10 .or. bigJ(place) == 11 .or. bigJ(place) == 13) then
                !if a reaction that releases a PLi occurs:

                !update the current doublet to its new state
                newstate = qstate(place, :) + V(:, bigJ(place))
                qstate(place, :) = 0   !reset qstate to 0 so we're ready for the next timestep

                do j = 1, 9
                    if (newstate(j) > 0) then
                        statetype = j
                        exit
                    end if
                end do

                if (statetype == 8) statetype = -1                        !if a doublet degrades, set its state to be -1

                state(col(place), row(place)) = statetype

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

                !figure out which PLi is released

                do iPLi = 1, countPLi
                    if (PLilocation(1, iPLi) == row(place) .and. PLilocation(2, iPLi) == col(place)) then
                        countcurPLi = countcurPLi + 1
                        currentPLib(countcurPLi) = iPLi
                    end if
                end do

                !!figure out how to move tPA and PLi

                connect_ind = 0
                countlat3 = 0
                do j = 1, nodes**2
                    if (Lat(j, row(place)) > 0) then    !find nearest neighbor sites to node tPA unbound from
                        countlat3 = countlat3 + 1
                        connect_ind(countlat3) = j
                    end if

                end do

                sizenonzero_connectind = countlat3

                tot_sitesPLi = 0
                prob_sitePLi = 0
                sites_N03 = 0
                sites_N10 = 0
                sites_Nplg = 0
                sites_Nplgn = 0
                sites_Nplgi = 0
                sites_degrade = 0

                alpha = 1  !parameter describing how likely PLi is to move onto a degraded doublet rather than an undegraded one.
                !0<alpha<=1. alpha=1 means PLi is equally likely to move onto degraded or undegraded.

                !if the reaction involved a PLi unbinding from an N33, then randomly choose which PLi is the one to move:

                if (countcurPLi == 2) then
                    randPLi = urcw1()

                    if (randPLi <= 0.5) then
                        currentPLi = currentPLib(1)
                    else
                        currentPLi = currentPLib(2)
                    end if
                elseif (countcurPLi == 1) then
                    currentPLi = currentPLib(1)
                    if (bigJ(place) == 7) write (*, *) 'Problem: there should be 2 PLis if rxn 7 occurs'
                end if

                call movepli(Ntot, nodes, sizenonzero_connectind, state, connect_ind, alpha, prob_N02i, prob_N00i, prob_N02n, &
                             prob_N00n, currentPLi, PLilocation, PLilocationt, statetemp)
                PLilocation = PLilocationt
                state = statetemp

            elseif (bigJ(place) == 12) then
                !if a reaction that converts PLG to PLi occurs
                !!!below is for tPA only converting PLG on the same doublet to PLi
                !countPLi=countPLi+1                         !keep track of # of PLi
                !PLilocation(1,countPLi)=row(place)          !location of current PLi
                !PLilocation(2,countPLi)=col(place)          !doublet PLi is on at location PLilocation(countPLi,1)
                !PLitime(countPLi)=t
                !
                !!!update the current doublet to its new state
                !newstate=qstate(place,:)+V(:,bigJ(place))
                !qstate(place,:)=0   !reset qstate to 0 so we're ready for the next timestep
                !
                !do j=1,9
                !   if(newstate(j)>0) then
                !      statetype=j
                !      exit
                !   end if
                !enddo
                !
                !if(statetype/=4) then
                !  write(*,*)'Problem: here the only state should be N13'
                !  write(*,*)'statetype=',statetype
                !end if
                !
                !state(col(place),row(place))=statetype
                !
                !count = count + 1
                !t = t + entry        !entry=tau(place);
                !
                !tvals(count) = t
                !Anumber(count)=A
                !reaction(count)=bigJ(place)
                !!!end section for tPA only converting PLG on the same doublet to PLi

                !!!below is for tPA converting any PLG on the same binding location to PLi

                countPLi = countPLi + 1                         !keep track of # of PLi
                PLilocation(1, countPLi) = row(place)          !location of current PLi

                call findint(state(:, row(place)), Ntot, 2, numP2) !find number N12
                call findint(state(2:Ntot, row(place)), Ntot - 1, 1, numP1) !find number newly exposed ("nicked") Nplg
                call findint(state(1, row(place)), 1, 1, numP1b) !find number initially exposed ("intact") Nplg
                call findint(state(:, row(place)), Ntot, 6, numP6) !find number N23
                numPt = (2*prob_N22n + prob_N02n)*numP1 + numP2 + numP6 + (2*prob_N22i + prob_N02i)*numP1b !find total number of PLG

                !figure out which PLG molecule gets converted, i.e. which doublet the PLG molecule is on
                probnumP(1) = (2*prob_N22n + prob_N02n)*dble(numP1)/numPt   !nicked Nplg
                probnumP(2) = dble(numP2)/numPt !N12
                probnumP(3) = dble(numP6)/numPt !N23
                probnumP(4) = (2*prob_N22i + prob_N02i)*dble(numP1b)/numPt  !intact Nplg
                call cumsum(probnumP, 4, cumprobnumP)

                rP = urcw1()
                do j = 1, 4
                    if (rP <= cumprobnumP(j)) then
                        colPLGb = j
                        exit
                    end if
                end do

                if (colPLGb == 3) then
                    colPLGfind = 6
                    call randchoice(state(:, row(place)), Ntot, colPLGfind, colPLG)
                elseif (colPLGb == 4) then
                    colPLG = 1  !if intact doublet, it must be the 1st doublet
                elseif (colPLGb == 2) then
                    colPLGfind = 2
                    call randchoice(state(:, row(place)), Ntot, colPLGfind, colPLG)
                elseif (colPLGb == 1) then
                    colPLGfind = 1
                    call randchoice(state(2:Ntot, row(place)), Ntot - 1, colPLGfind, colPLG)
                    colPLG = colPLG + 1  !because we searched entries 2 to 6, but it registers at 1-5
                end if

                PLilocation(2, countPLi) = colPLG          !doublet PLi is on at location PLilocation(countPLi,1)
                PLitime(countPLi) = t
                if (countPLi == 1) then
                    countfp = countfp + 1
                    !firstPLi(countfp)=t+entry
                    firstPLi(stats) = t + entry
                    !write(*,*) 'countfp=',countfp
                    !write(*,*) 'stats=',stats
                end if

                !!update the current doublet to its new state
                !newstate=qstate(place,:)+V(:,bigJ(place))
                qstate(place, :) = 0   !reset qstate to 0 so we're ready for the next timestep

                !change the PLG to PLi:
                if (state(colPLG, row(place)) == 1) then
                    !if state of doublet is Nplg, new state is either statetype=6 or statetype=5 figure out which
                    rP = urcw1()
                    if (colPLG == 1) then
                        if (rP <= prob_N02i/(prob_N02i + prob_N22i)) then
                            statetype = 5
                        else
                            statetype = 6
                        end if
                    else
                        if (rP <= prob_N02n/(prob_N02n + prob_N22n)) then
                            statetype = 5
                        else
                            statetype = 6
                        end if
                    end if
                elseif (state(colPLG, row(place)) == 2) then
                    statetype = 4 !if state of doublet is N12, new state is N13
                elseif (state(colPLG, row(place)) == 6) then
                    statetype = 7 !if state of doublet is N23, new state is N33
                else
                    write (*, *) 'Problem: should be a doublet with PLG'
                end if

                state(colPLG, row(place)) = statetype

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)
                !!!end section for tPA converting any PLG on the same binding location to PLi

            elseif (bigJ(place) == 4 .or. bigJ(place) == 8 .or. bigJ(place) == 9) then
                !if a rxn that releases 1 tPA occurs:
                reaction2 = reaction2 + 1

                !update the current doublet to its new state
                newstate = qstate(place, :) + V(:, bigJ(place))
                qstate(place, :) = 0   !reset qstate to 0 so we're ready for the next timestep

                do j = 1, 9
                    if (newstate(j) > 0) then
                        statetype = j
                        exit
                    end if
                end do

                if (statetype == 8) then
                    statetype = -1                        !if a doublet degrades, set its state to be -1
                    write (*, *) 'Problem: degradation should not occur here'
                end if

                state(col(place), row(place)) = statetype

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

                !!figure out how to move tPA. I'm not sure I allow tPA to move (if it unbinds, it's gone forever). If I do let it move, I'll need to look more carefully here...

                connect_ind = 0
                countlat3 = 0
                do j = 1, nodes**2
                    if (Lat(j, row(place)) > 0) then    !find nearest neighbor sites to node tPA unbound from
                        countlat3 = countlat3 + 1
                        connect_ind(countlat3) = j
                    end if

                end do

                sizenonzero_connectind = countlat3

                tot_sitesPLi = 0
                prob_sitePLi = 0
                sites_N03 = 0
                sites_N10 = 0
                sites_Nplg = 0
                sites_Nplgn = 0
                sites_Nplgi = 0
                sites_degrade = 0

                !first calculate the rebinding probability
                call movetpa(Ntot, nodes, state, prob_N02i, prob_N00i, prob_N02n, prob_N00n, vol3, ktPAon, D, count, p_rebind)
                !! BRAD             write(*,*)'p_rebind=',p_rebind

            elseif (bigJ(place) == 2 .or. bigJ(place) == 5 .or. bigJ(place) == 1 .or. bigJ(place) == 6) then
                !if a rxn that results in a binding or unbinding of PLG occurs

                !update the current doublet to its new state
                newstate = qstate(place, :) + V(:, bigJ(place))
                qstate(place, :) = 0   !reset qstate to 0 so we're ready for the next timestep

                do j = 1, 9
                    if (newstate(j) > 0) then
                        statetype = j
                        exit
                    end if
                end do

                if (statetype == 8) then
                    statetype = -1                        !if a doublet degrades, set its state to be -1
                    write (*, *) 'Problem: degradation should not occur here'
                end if

                state(col(place), row(place)) = statetype

                count = count + 1
                t = t + entry        !entry=tau(place);

                tvals(count) = t
                Anumber(count) = A
                reaction(count) = bigJ(place)

            end if !for if(bigJ(place)==14....)
            bigJ(place) = 0

            !keep track of how much tPA is in system at each step, and where it is in space

            if (tPA(count - 1) == 1) then
                countstate9 = 0
                do i = 1, Ntot
                    do j = 1, nodes**2
                        if (state(i, j) >= 2 .and. state(i, j) <= 4) then
                            countstate9 = countstate9 + 1
                            rowtPA(countstate9) = j
                        end if
                    end do
                end do

                tPA(count) = countstate9
                if (tPA(count) > 1) write (*, *) 'should only have 1 tPA in system. something is wrong'
                if (tPA(count) == 1) tPA_loc(count) = rowtPA(countstate9)
            end if

            !keep track of how much PLG is in the system at each step
            countstate13 = 0
            do i = 1, Ntot
                do j = 1, nodes**2
                    if (state(i, j) == 2 .or. state(i, j) == 6) countstate13 = countstate13 + 1
                end do
            end do

            countstate14 = 0
            do j = 1, nodes**2
                if (state(1, j) == 1) countstate14 = countstate14 + 1
            end do

            countstate15 = 0
            do i = 2, Ntot
                do j = 1, nodes**2
                    if (state(i, j) == 1) countstate15 = countstate15 + 1
                end do
            end do

            PLG(count) = countstate13 + (prob_N02i + 2*prob_N22i)*countstate14 + (prob_N02n + 2*prob_N22n)*countstate15

            !keep track of how much PLi is in system at each step

            countstate10 = 0
            do i = 1, Ntot
                do j = 1, nodes**2
                    if (state(i, j) >= 4 .and. state(i, j) <= 6) then
                        countstate10 = countstate10 + 1
                        rowPLi(countstate10) = j
                        colPLi(countstate10) = i
                    end if
                end do
            end do

            countstate10a = 0
            do i = 1, Ntot
                do j = 1, nodes**2
                    if (state(i, j) == 9) then
                        countstate10a = countstate10a + 1
                        rowPLi(countstate10a) = j
                        colPLi(countstate10a) = i
                    end if
                end do
            end do

            countstate11 = 0
            do i = 1, Ntot
                do j = 1, nodes**2
                    if (state(i, j) == 7) then
                        countstate11 = countstate11 + 1
                        rowPLi(countstate11) = j
                        colPLi(countstate11) = i
                    end if
                end do
            end do

            PLi(count) = countstate10 + 2*countstate11 + countstate10a

            if (PLi(count) .lt. PLi(count - 1)) then
                write (*, *) 'Problem: should not lose any PLi from system. something is wrong.'
                write (*, *) 'count=', count
                write (*, *) 'PLi(count)=', PLi(count)
                write (*, *) 'reaction(count)=', reaction(count)
                write (*, *) 'PLi(count-1)=', PLi(count - 1)
                write (*, *) 'reaction(count-1)=', reaction(count - 1)
                exit
            end if

            timevect(count) = t

            countstate12 = 0
            do i = 1, Ntot
                do j = 1, nodes**2
                    if (state(i, j) .gt. -1 .and. state(i, j) .lt. 9) then
                        countstate12 = countstate12 + 1 !count up the number of undegraded doublets, i.e, those doublets not in state -1 or 9
                        row2(countstate12) = j
                        col2(countstate12) = i
                    end if
                end do
            end do

            Tdoublets2 = 6.0*nodes**2 !total number of doublets is 6*nodes^2 because there are 6 doublets at each node
            percent_degrade(count) = 1 - countstate12/Tdoublets2 !the percent of doublets that have been degraded is 1 minus the fraction of undegraded doublets
            if (mod(count, runs) == 0) then
                write (*, *) 'percent_degrade=', percent_degrade(count)
                write (*, *) 't=', t
            end if

            !!!!!!!!! Added 5/4/11 to calculate time to reach a given percentage of degraded binding doublets
            !Ninteger=int(percent_degrade(count)*100)
            !
            !
            !if(Ninteger>Nsave) then !if the current percentage is the 1st past a new a whole percent, save the current time
            !Nsave=Nsave+1
            !persave(Nsave+1) = percent_degrade(count)
            !tsave(Nsave+1) = t
            !end if
            !!!!!!!

            if (percent_degrade(count) > 2.0/3.0) then
                lysiscomplete(stats) = 1
                exit    !if more than 2/3 of the total number of doublets have been
                !degraded, stop the algorithm because we consider this to be lysis
            end if

        end do !enddo for big "do" Gillespie loop

        !!!!!uncomment below 2 lines, as well as above 5/4/11 stuff if I want to save percent degraded
        !write(permit) persave
        !write(t2unit) tsave

        lysis_time(stats) = tvals(count)                       !saves the time for lysis to complete
        if (tPA(count) == 1) then
            tPA_time(stats) = tvals(count)
        else
            call findfirst(tPA, 1000000, 0, loca)
            tPA_time(stats) = tvals(loca)                      !saves the time tPA leaves the system
        end if

        Plasmin(stats) = PLi(count)                          !saves the total number of PLi generated
        max_Plg(stats) = maxval(PLG)                            !saves the maximum number of PLG generated
        countvect(stats) = count
        tPAunbind(stats) = reaction2
        tPAPLiunbd(stats) = reaction3
        ltPA(stats) = tPA(count)

        if (stats == runs) then   !change back to 10000 when I change runs
            write (lysunit) lysis_time(1:stats)  !writes the lysis_time to file
            write (tPAunit) tPA_time(1:stats)
            write (PLiunit) Plasmin(1:stats)
            write (endunit) lysiscomplete(1:stats)
            !write(plgbdunit) PLGbd(1:stats)
            !write(plgunbdunit) PLGunbd(1:stats)
            write (tPAunbdunit) tPAunbind(1:stats)
            write (tPAPLiunit) tPAPLiunbd(1:stats)
            write (s2unit) ltPA(1:stats)
            write (fpunit) firstPLi(1:stats)
            close (lysunit)
            close (tPAunit)
            close (PLiunit)
            close (endunit)
            !close(plgbdunit)
            !close(plgunbdunit)
            !close(sunit)
            close (tPAunbdunit)
            !close(perunit)
            close (s2unit)
            !close(t2unit)
            close (tPAPLiunit)
            close (fpunit)
            write (*, *) ' kncat=', kncat
            write (*, *) ' kapcat=', kapcat
            write (*, *) ' ktPAon=', ktPAon
            write (*, *) ' kaoff10=', kaoff10
            write (*, *) ' kaoff12=', kaoff12
            write (*, *) ' kplioff=', kplioff
            write (*, *) ' kplgoffnick=', kplgoffnick
            write (*, *) ' kplgon=', kplgon
            write (*, *) ' kplgoff=', kplgoff
            write (*, *) ' freeplg=', freeplg
            write (*, *) ' kdeg=', kdeg
        end if

    end do   !enddo for stats loop

    !close(plgunit)
    !close(ctunit)
    !close(plitimeunit)

CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE CONCAT - concatenates a vector and a single number
    !
    subroutine concat(g, size2, h, sizeg, gnew)

        integer  :: sizeg, size2
        double precision  :: h
        double precision, dimension(size2)   :: g, gnew

        gnew(1:sizeg) = g(1:sizeg)
        gnew(sizeg + 1) = h

    end subroutine concat

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE HORZCONCAT - concatenates 2 vectors
    !
    subroutine horzconcat(g, sizeg, h, sizeh, gnew)

        integer  :: sizeh, sizeg
        double precision, dimension(sizeh)   :: h
        double precision, dimension(sizeg)   :: g
        double precision, dimension(sizeg + sizeh)  :: gnew

        gnew(1:sizeg) = g(1:sizeg)
        gnew(sizeg + 1:sizeg + sizeh) = h(1:sizeh)

    end subroutine horzconcat

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE VURCW1 - 
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
    !SUBROUTINE FINDZERO - finds all the zero elements in a vector, and returns their index number in a new vector
    !
    subroutine findzero(g, sizeg, sizegnew)!,gnew)

        integer  :: sizeg, sizegnew
        integer  :: i, countg
        double precision, dimension(sizeg)   :: g
        integer, dimension(sizegnew)  :: gnew

        countg = 0

        do i = 1, sizeg
            if (g(i) == 0) then
                countg = countg + 1
                gnew(countg) = i
            end if

        end do

        sizegnew = countg

    end subroutine findzero

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE FINDZEROINT - finds all the zero elements in a vector, and returns their index number in a new vector
    !
    subroutine findzeroint(g, sizeg, sizegnew)!,gnew)

        integer  :: sizeg, sizegnew
        integer  :: i, countg
        integer, dimension(sizeg)   :: g

        countg = 0

        do i = 1, sizeg
            if (g(i) == 0) then
                countg = countg + 1
            end if

        end do

        sizegnew = countg

    end subroutine findzeroint

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE FINDINT - finds number of entries in vector g equal to intspec
    !
    subroutine findint(g, sizeg, intspec, sizegnew)
        integer  :: sizeg, sizegnew
        integer  :: i, countg, intspec
        integer, dimension(sizeg)   :: g

        countg = 0

        do i = 1, sizeg
            if (g(i) == intspec) then
                countg = countg + 1
            end if

        end do

        sizegnew = countg

    end subroutine findint

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE FINDFIRST - finds the first entry equal to a specified number
    !
    subroutine findfirst(g, sizeg, numspec, loca)
        integer  :: sizeg
        integer  :: i, numspec, loca
        integer, dimension(sizeg)   :: g

        loca = 0
        do i = 1, sizeg
            if (g(i) == numspec) then
                loca = i
                exit
            end if

        end do

    end subroutine findfirst

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE RANDCHOICE - finds the entries equal to a specified number, and randomly chooses one of them to be used later
    !
    subroutine randchoice(g, sizeg, numspec, loca)

        integer  :: sizeg

        integer  :: i, numspec, loca, count, sizel

        integer, dimension(sizeg)   :: g
        double precision :: r
        integer, dimension(sizeg) :: location
        double precision, dimension(sizeg) :: prob, cumprob

        count = 0
        prob = 0
        cumprob = 0

        do i = 1, sizeg

            if (g(i) == numspec) then

                count = count + 1
                location(count) = i
            end if

        end do

        sizel = count
        if (sizel == 1) then
            loca = location(count)
        elseif (sizel == 0) then
            loca = 0
            write (*, *) 'Problem:should only access subroutine if at least one entry is this number'
        else
            do i = 1, sizel
                prob(i) = 1.0/sizel  !have it be equally likely to choose any entry
            end do
            call cumsum(prob(1:sizel), sizel, cumprob(1:sizel))

            r = urcw1()

            do i = 1, sizel
                if (r <= cumprob(i)) then
                    loca = location(i)
                    exit
                end if
            end do
        end if

    end subroutine randchoice

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE CUMSUM - returns a vector with entries cumulative sums of original vector
    !
    subroutine cumsum(g, sizeg, gnew)

        integer  :: sizeg
        integer  :: i
        double precision :: currentsum
        double precision, dimension(sizeg)   :: g, gnew

        currentsum = 0

        do i = 1, sizeg
            currentsum = currentsum + g(i)
            gnew(i) = currentsum
        end do

    end subroutine cumsum

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE MINVECT - finds location of minimum entry in vector
    !
    subroutine minvect(g, sizeg, minlocation)

        integer  :: sizeg, minlocation
        double precision, dimension(sizeg)   :: g
        double precision  :: g_temp
        integer  :: i, g_old

        g_old = 1

        do i = 1, sizeg
            g_temp = g(i)
            if (g_temp .lt. g(g_old) .and. g_temp /= 0) then
                g_old = i
            end if
        end do

        minlocation = g_old

    end subroutine minvect

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE MOVETPA - moves one tPA molecule. the only output is "p_rebind"
    !
    subroutine movetpa(Ntot, nodes, state, prob_N02i, prob_N00i, prob_N02n, prob_N00n, vol3, param9, D, count, p_rebind)

        integer :: Ntot, nodes
        integer :: countstate6, countstate7, countstate8, count, j, jj
        integer, dimension(Ntot, nodes**2) :: state
        double precision :: prob_N02i, prob_N00i, prob_N02n, prob_N00n, vol3, exposed, bs0
        double precision :: r0, r1, k0, param9, p_rebind, D
        double precision :: rr

        !calculate the rebinding probability by first calculating the number of exposed binding sites
        !number of N03 doublets
        countstate6 = 0
        do j = 1, Ntot
            do jj = 1, nodes**2
                if (state(j, jj) == 5) countstate6 = countstate6 + 1
            end do
        end do

        !number of intact Nplg
        countstate7 = 0
        do jj = 1, nodes**2
            if (state(1, jj) == 1) countstate7 = countstate7 + 1
        end do

        !number of nicked Nplg
        countstate8 = 0
        do jj = 1, nodes**2
            do j = 2, Ntot
                if (state(j, jj) == 1) countstate8 = countstate8 + 1
            end do
        end do

        exposed = countstate6 + (prob_N02i + 2*prob_N00i)*countstate7 + (prob_N02n + 2*prob_N00n)*countstate8

        bs0 = exposed/vol3/602.2

        !if tPA molecule can escape in 3 dimensions:
        r0 = 0.5    !units nm
        r1 = 22.5   !units nm
        k0 = param9*bs0  !binding rate

        p_rebind = 1 - r1/r0*(1/(1 + (r1 - r0)*sqrt(k0/D)*cosh(sqrt(k0/D)*r0)/sinh(sqrt(k0/D)*r0)))
        !write(*,*)'p_rebind=',p_rebind

        rr = urcw1()

        !!!!! IF REBINDING PROB IS BIGGER THAN ZERO, THEN I NEED TO ADD A WHOLE SECTION HERE PERTAINING TO THE REBINDING.
        !!!!! I STARTED IT BELOW (MARKED WITH "!!") BUT STOPPED BECAUSE I DON'T NEED IT YET. IF I EVENTUALLY NEED IT,
        !!!!! I NEED TO ADD LINES 435-521 FROM MACRO_NEW/FINE/MICRO_CHANGEPLI_NEW.M    (8/13/10)
        if (rr <= p_rebind) then
            write (*, *) 'Problem: random number smallerer than p_rebind. must add rebinding to code'
            !exit
        end if
        ! if(rr<=p_rebind)then
        !         !if random number is less than or equal to p_rebind, have tPA rebind to any available doublet.
        !         !1st have it choose which site to go to

        !         connect_ind2=0    !zeros(length(connect_ind),5);

        !         do i=1,sizenonzero_connectind
        !             countlat2=0
        !             connect_ind2temp=0
        !             do j=nodes**2
        !             if(Lat(j,connect_ind(i))/=0) then      !find sites 2 lengths away from node tPA unbound from
        !                 countlat2=countlat2+1
        !                 connect_ind2temp(countlat2) = j
        !             end if
        !             enddo
        !             connect_ind2(1:countlat2,i)=connect_ind2temp(1:countlat2)
        !     enddo
        !         !look at the matrix connect_ind2 for repeated entries. For any node
        !         !number that appears more than once, reset all the entries with that
        !         !node number = 0 except for one. So we're left with a matrix in which
        !         !the only nonzero entries are the nodes we can go to, and none of
        !         !these numbers are repeated.
        !         do i=1,length(connect_ind2(1,:))
        !             for j=1:5
        !                 [rtest,ctest]=find(connect_ind2==connect_ind2(i,j));
        !                 if length(rtest) > 1 && connect_ind2(i,j) > 0
        !                     for jj=1:length(rtest)-1
        !                         connect_ind2(rtest(jj),ctest(jj))=0;
        !                     end
        !                 end
        !             end
        !         end


    end subroutine movetpa

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE MOVEPLI - moves one PLi molecule. the output is "state" and "PLilocationt"
    !
    subroutine movepli(Ntot, nodes, sizenonzero_connectind, state, connect_ind, alpha, prob_N02i, prob_N00i, prob_N02n, &
                       prob_N00n, loc, PLilocation, PLilocationt, statetemp)

        integer :: sizenonzero_connectind, Ntot, nodes
        integer, dimension(Ntot, nodes**2) :: state, statetemp
        integer, dimension(6)  :: connect_ind, sites_N10, sites_N03, sites_Nplgn, sites_degrade, sites_Nplgi
        integer :: j, temp
        double precision, dimension(6)  :: tot_sitesPLi, prob_sitePLi, cumprob
        double precision :: prob_N02i, prob_N00i, prob_N02n, prob_N00n, alpha
        double precision :: rin
        integer :: newsite2, newdoublet2, loc
        double precision, dimension(5)  :: P, cumP
        integer :: N10toN13, NplgtoN23, N03toN33, degradetoPLi
        integer, dimension(2, 300)  :: PLilocation, PLilocationt
        double precision, dimension(2)  :: Pnsmall, cumPnsmall, Pismall, cumPismall

        !move PLi

        newsite2 = 0
        newdoublet2 = 0
        prob_sitePLi = 0
        cumprob = 0
        tot_sitesPLi = 0
        PLilocationt = PLilocation

        do j = 1, sizenonzero_connectind
            call findint(state(:, connect_ind(j)), Ntot, 3, sites_N10(j))
            call findint(state(:, connect_ind(j)), Ntot, 5, sites_N03(j))
            call findint(state(:, connect_ind(j)), Ntot, -1, sites_degrade(j))  !will assume each degraded doublet has
            !room for 1 PLi to bind
            call findint(state(2:6, connect_ind(j)), 5, 1, sites_Nplgn(j))      ! # of nicked doublets in state 1
            if (state(1, connect_ind(j)) == 1) then
                sites_Nplgi(j) = 1                            ! # of intact doublets in state 1
            else !if(state(1,connect_ind(j))==0)then
                sites_Nplgi(j) = 0
            end if

            tot_sitesPLi(j) = sites_N10(j) + sites_N03(j) + (prob_N02n + 2*prob_N00n)*sites_Nplgn(j) + &
                              (prob_N02i + 2*prob_N00i)*sites_Nplgi(j) + alpha*sites_degrade(j)
            !so tot_situsPLi(j) is the total number of empty binding sites for PLi at location j.
            !it's not an integer because we only know the probability of Nplg being N02 or N00.

        end do

        do j = 1, sizenonzero_connectind
            prob_sitePLi(j) = tot_sitesPLi(j)/sum(tot_sitesPLi)      !prob. PLi goes to site j
        end do

        call cumsum(prob_sitePLi, 6, cumprob)
        rin = urcw1()
        do j = 1, 6
            if (rin <= cumprob(j)) then
                newsite2 = j
                exit
            end if
        end do

        do j = 1, sizenonzero_connectind
            if (newsite2 == j) then
                P(1) = sites_N10(j)/tot_sitesPLi(j)                            !prob. of moving onto N10
                P(2) = (prob_N02n + 2*prob_N00n)*sites_Nplgn(j)/tot_sitesPLi(j)  !prob. of moving onto nicked Nplg
                P(3) = (prob_N02i + 2*prob_N00i)*sites_Nplgi(j)/tot_sitesPLi(j)  !prob. of moving onto intact Nplg
                P(4) = sites_N03(j)/tot_sitesPLi(j)                            !prob of moving onto N03
                P(5) = alpha*sites_degrade(j)/tot_sitesPLi(j)                  !prob. of moving onto degraded doublet
            elseif (newsite2 > sizenonzero_connectind) then
                write (*, *) 'Problem: newsite2'
                exit
            end if
        end do

        call cumsum(P, 5, cumP)
        rin = urcw1()
        do j = 1, 5
            if (rin <= cumP(j)) then      !find the doublet PLi will move to
                newdoublet2 = j
                exit
            end if
        end do

        if (newdoublet2 == 1) then                      !if PLi is going to N10
            call findfirst(state(:, connect_ind(newsite2)), Ntot, 3, N10toN13)  !find the next 3 entry in state vector
            state(N10toN13, connect_ind(newsite2)) = 4                         !change the state to 4, N13
            PLilocationt(1, loc) = connect_ind(newsite2)                        !location of current PLi
            PLilocationt(2, loc) = N10toN13

        elseif (newdoublet2 == 2) then                  !if PLi is going to nicked Nplg
            call findfirst(state(2:6, connect_ind(newsite2)), 5, 1, NplgtoN23)
            NplgtoN23 = 1 + NplgtoN23                    !find the next 1 entry in columns 2-6 of state vector
            if (NplgtoN23 == 1) then
                write (*, *) 'Problem: NplgtoN23'
            end if
            !need the 1+ in above defn. because we're searcing 5 columns, but we want to start with column 2.
            !So if the first 1 is in column 2, NplgtoN23 would equal 1. We need to add 1 to it.
            Pnsmall(1) = prob_N00n/(prob_N02n + prob_N00n)
            Pnsmall(2) = prob_N02n/(prob_N02n + prob_N00n)
            call cumsum(Pnsmall, 2, cumPnsmall)
            rin = urcw1()
            do j = 1, 2
                if (rin <= cumPnsmall(j)) then      !find the doublet PLi will move to
                    temp = j
                    exit
                end if
            end do

            if (temp == 1) then
                state(NplgtoN23, connect_ind(newsite2)) = 5           !change the state to 5, N03
            elseif (temp == 2) then
                state(NplgtoN23, connect_ind(newsite2)) = 6           !change the state to 6, N23
            end if
            PLilocationt(1, loc) = connect_ind(newsite2)               !location of current PLi
            PLilocationt(2, loc) = NplgtoN23

        elseif (newdoublet2 == 3) then                   !if PLi is going to intact Nplg
            NplgtoN23 = 1                               !only place we have intact fibrin is in column 1
            Pismall(1) = prob_N00i/(prob_N02i + prob_N00i)
            Pismall(2) = prob_N02i/(prob_N02i + prob_N00i)
            call cumsum(Pismall, 2, cumPismall)
            rin = urcw1()

            do j = 1, 2
                if (rin <= cumPismall(j)) then         !find the doublet PLi will move to
                    temp = j
                    exit
                end if
            end do

            if (temp == 1) then
                state(NplgtoN23, connect_ind(newsite2)) = 5            !change the state to 5, N03

            elseif (temp == 2) then
                state(NplgtoN23, connect_ind(newsite2)) = 6            !change the state to 6, N23
            end if
            PLilocationt(1, loc) = connect_ind(newsite2)                !location of current PLi
            PLilocationt(2, loc) = NplgtoN23

        elseif (newdoublet2 == 4) then                  !if PLi is going to N03
            call findfirst(state(:, connect_ind(newsite2)), Ntot, 5, N03toN33)  !find the next 5 entry in state vector
            state(N03toN33, connect_ind(newsite2)) = 7                 !change the state to 7, N33
            PLilocationt(1, loc) = connect_ind(newsite2)                !location of current PLi
            PLilocationt(2, loc) = N03toN33

        elseif (newdoublet2 == 5) then                  !if PLi is going to degraded doublet
            call findfirst(state(:, connect_ind(newsite2)), Ntot, -1, degradetoPLi) !find the next -1 entry in state vector
            state(degradetoPLi, connect_ind(newsite2)) = 9        !change the state to 9, degraded doublet with PLi bound
            PLilocationt(1, loc) = connect_ind(newsite2)           !location of current PLi
            PLilocationt(2, loc) = degradetoPLi

        end if

        statetemp = state

    end subroutine movepli

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !SUBROUTINE MOVE2PLI - moves 2 PLi molecules. the output is "state" and "PLilocationt"
    !
    subroutine move2pli(Ntot, nodes, sizenonzero_connectind, state, connect_ind, alpha, prob_N02i, prob_N00i, prob_N02n, &
                        prob_N00n, loc, loc2, PLilocation, PLilocationt, statetemp)

        integer :: sizenonzero_connectind, Ntot, nodes
        integer, dimension(Ntot, nodes**2) :: state, statetemp
        integer, dimension(6)  :: connect_ind, sites_N10, sites_N03, sites_Nplgn, sites_degrade, sites_Nplgi
        integer :: j, temp
        double precision, dimension(6)  :: tot_sitesPLi, prob_sitePLi, cumprob
        double precision :: prob_N02i, prob_N00i, prob_N02n, prob_N00n, alpha
        double precision :: rin
        integer :: newsite, newdoublet, newsite2, newdoublet2, loc, loc2
        double precision, dimension(5)  :: P, cumP
        integer :: N10toN13, NplgtoN23, N03toN33, degradetoPLi
        integer, dimension(2, 300)  :: PLilocation, PLilocationt
        double precision, dimension(2)  :: Pnsmall, cumPnsmall, Pismall, cumPismall

        newsite = 0
        newdoublet = 0
        newsite2 = 0
        newdoublet2 = 0
        prob_sitePLi = 0
        cumprob = 0
        tot_sitesPLi = 0
        PLilocationt = PLilocation

        !move one PLi first
        do j = 1, sizenonzero_connectind

            call findint(state(:, connect_ind(j)), Ntot, 3, sites_N10(j))
            call findint(state(:, connect_ind(j)), Ntot, 5, sites_N03(j))
            call findint(state(:, connect_ind(j)), Ntot, -1, sites_degrade(j))  !will assume each degraded doublet has
            !room for 1 PLi to bind
            call findint(state(2:6, connect_ind(j)), 5, 1, sites_Nplgn(j))      ! # of nicked doublets in state 1
            if (state(1, connect_ind(j)) == 1) then
                sites_Nplgi(j) = 1                            ! # of intact doublets in state 1
            else !if(state(1,connect_ind(j))==0)then
                sites_Nplgi(j) = 0
            end if

        tot_sitesPLi(j) = sites_N10(j) + sites_N03(j) + (prob_N02n+2*prob_N00n)*sites_Nplgn(j) + (prob_N02i+2*prob_N00i)*sites_Nplgi(j) + alpha*sites_degrade(j)
            !so tot_sitesPLi(j) is the total number of empty binding sites for PLi at location j.
            !it's not an integer because we only know the probability of Nplg being N02 or N00.

        end do

        do j = 1, sizenonzero_connectind
            prob_sitePLi(j) = tot_sitesPLi(j)/sum(tot_sitesPLi)         !prob. PLi go to site j
        end do

        !randomly choose which site we go to
        call cumsum(prob_sitePLi, 6, cumprob)
        rin = urcw1()

        do j = 1, 6
            if (rin <= cumprob(j)) then
                newsite = j
                exit
            end if
        end do

        do j = 1, sizenonzero_connectind
            if (newsite == j) then
                P(1) = sites_N10(j)/tot_sitesPLi(j)                            !prob. of moving onto N10
                P(2) = (prob_N02n + 2*prob_N00n)*sites_Nplgn(j)/tot_sitesPLi(j)  !prob. of moving onto nicked Nplg
                P(3) = (prob_N02i + 2*prob_N00i)*sites_Nplgi(j)/tot_sitesPLi(j)  !prob. of moving onto intact Nplg
                P(4) = sites_N03(j)/tot_sitesPLi(j)                            !prob of moving onto N03
                P(5) = alpha*sites_degrade(j)/tot_sitesPLi(j)                  !prob. of moving onto degraded doublet
            elseif (newsite > sizenonzero_connectind) then
                write (*, *) 'Problem: newsite PLi2'
                exit
            end if
        end do

        call cumsum(P, 5, cumP)
        rin = urcw1()

        do j = 1, 5
            if (rin <= cumP(j)) then      !randomly chooses the doublet PLi will move to
                newdoublet = j
                exit
            end if
        end do

        if (newdoublet == 1) then                      !if PLi is going to N10
            call findfirst(state(:, connect_ind(newsite)), Ntot, 3, N10toN13)  !find the next 3 entry in state vector
            state(N10toN13, connect_ind(newsite)) = 4                         !change the state to 4, N13
            PLilocationt(1, loc) = connect_ind(newsite)                        !location of current PLi
            PLilocationt(2, loc) = N10toN13

        elseif (newdoublet == 2) then                  !if PLi is going to nicked Nplg
            call findfirst(state(2:6, connect_ind(newsite)), 5, 1, NplgtoN23)
            NplgtoN23 = 1 + NplgtoN23                    !find the next 1 entry in columns 2-6 of state vector
            if (NplgtoN23 == 1) then
                write (*, *) 'Problem: NplgtoN23'
            end if
            !need the 1+ in above defn. because we're searcing 5 columns, but we want to start with column 2.
            !So if the first 1 is in column 2, NplgtoN23 would equal 1. We need to add 1 to it.
            Pnsmall(1) = prob_N00n/(prob_N02n + prob_N00n)
            Pnsmall(2) = prob_N02n/(prob_N02n + prob_N00n)
            call cumsum(Pnsmall, 2, cumPnsmall)
            rin = urcw1()

            do j = 1, 2
                if (rin <= cumPnsmall(j)) then      !find the doublet PLi will move to
                    temp = j
                    exit
                end if
            end do

            if (temp == 1) then
                state(NplgtoN23, connect_ind(newsite)) = 5           !change the state to 5, N03
            elseif (temp == 2) then
                state(NplgtoN23, connect_ind(newsite)) = 6           !change the state to 6, N23
            end if
            PLilocationt(1, loc) = connect_ind(newsite)               !location of current PLi
            PLilocationt(2, loc) = NplgtoN23

        elseif (newdoublet == 3) then                    !if PLi is going to intact Nplg
            NplgtoN23 = 1                               !only place we have intact fibrin is in column 1
            Pismall(1) = prob_N00i/(prob_N02i + prob_N00i)
            Pismall(2) = prob_N02i/(prob_N02i + prob_N00i)
            call cumsum(Pismall, 2, cumPismall)
            rin = urcw1()

            do j = 1, 2
                if (rin <= cumPismall(j)) then         !find the doublet PLi will move to
                    temp = j
                    exit
                end if
            end do

            if (temp == 1) then
                state(NplgtoN23, connect_ind(newsite)) = 5            !change the state to 5, N03

            elseif (temp == 2) then
                state(NplgtoN23, connect_ind(newsite)) = 6            !change the state to 6, N23
            end if

            PLilocationt(1, loc) = connect_ind(newsite)                !location of current PLi
            PLilocationt(2, loc) = NplgtoN23

        elseif (newdoublet == 4) then                  !if PLi is going to N03
            call findfirst(state(:, connect_ind(newsite)), Ntot, 5, N03toN33)  !find the next 5 entry in state vector
            state(N03toN33, connect_ind(newsite)) = 7                 !change the state to 7, N33
            PLilocationt(1, loc) = connect_ind(newsite)                !location of current PLi
            PLilocationt(2, loc) = N03toN33

        elseif (newdoublet == 5) then                  !if PLi is going to degraded doublet
            call findfirst(state(:, connect_ind(newsite)), Ntot, -1, degradetoPLi) !find the next -1 entry in state vector
            state(degradetoPLi, connect_ind(newsite)) = 9        !change the state to 9, degraded doublet with PLi bound
            PLilocationt(1, loc) = connect_ind(newsite)           !location of current PLi
            PLilocationt(2, loc) = degradetoPLi

        end if

        !then move the other PLi

        do j = 1, sizenonzero_connectind
            call findint(state(:, connect_ind(j)), Ntot, 3, sites_N10(j))
            call findint(state(:, connect_ind(j)), Ntot, 5, sites_N03(j))
            call findint(state(:, connect_ind(j)), Ntot, -1, sites_degrade(j))  !will assume each degraded doublet has
            !room for 1 PLi to bind
            call findint(state(2:6, connect_ind(j)), 5, 1, sites_Nplgn(j))      ! # of nicked doublets in state 1
            if (state(1, connect_ind(j)) == 1) then
                sites_Nplgi(j) = 1                            ! # of intact doublets in state 1
            else
                sites_Nplgi(j) = 0
            end if

        tot_sitesPLi(j) = sites_N10(j) + sites_N03(j) + (prob_N02n+2*prob_N00n)*sites_Nplgn(j) + (prob_N02i+2*prob_N00i)*sites_Nplgi(j) + alpha*sites_degrade(j)
            !so tot_sitesPLi(j) is the total number of empty binding sites for PLi at location j.
            !it's not an integer because we only know the probability of Nplg being N02 or N00.
            !           Nplgn = 0    !reset N

        end do

        do j = 1, sizenonzero_connectind
            prob_sitePLi(j) = tot_sitesPLi(j)/sum(tot_sitesPLi)      !prob. PLi goes to site j
        end do

        call cumsum(prob_sitePLi, 6, cumprob)
        rin = urcw1()

        do j = 1, 6
            if (rin <= cumprob(j)) then
                newsite2 = j
                exit
            end if
        end do

        do j = 1, sizenonzero_connectind
            if (newsite2 == j) then
                P(1) = sites_N10(j)/tot_sitesPLi(j)                            !prob. of moving onto N10
                P(2) = (prob_N02n + 2*prob_N00n)*sites_Nplgn(j)/tot_sitesPLi(j)  !prob. of moving onto nicked Nplg
                P(3) = (prob_N02i + 2*prob_N00i)*sites_Nplgi(j)/tot_sitesPLi(j)  !prob. of moving onto intact Nplg
                P(4) = sites_N03(j)/tot_sitesPLi(j)                            !prob of moving onto N03
                P(5) = alpha*sites_degrade(j)/tot_sitesPLi(j)                  !prob. of moving onto degraded doublet
            elseif (newsite2 > sizenonzero_connectind) then
                write (*, *) 'Problem: newsite2'
                exit
            end if
        end do

        call cumsum(P, 5, cumP)
        rin = urcw1()

        do j = 1, 5
            if (rin <= cumP(j)) then      !find the doublet PLi will move to
                newdoublet2 = j
                exit
            end if
        end do

        if (newdoublet2 == 1) then                      !if PLi is going to N10
            call findfirst(state(:, connect_ind(newsite2)), Ntot, 3, N10toN13)  !find the next 3 entry in state vector
            state(N10toN13, connect_ind(newsite2)) = 4                         !change the state to 4, N13
            PLilocationt(1, loc2) = connect_ind(newsite2)                        !location of current PLi
            PLilocationt(2, loc2) = N10toN13

        elseif (newdoublet2 == 2) then                  !if PLi is going to nicked Nplg
            call findfirst(state(2:6, connect_ind(newsite2)), 5, 1, NplgtoN23)
            NplgtoN23 = 1 + NplgtoN23                    !find the next 1 entry in columns 2-6 of state vector
            if (NplgtoN23 == 1) then
                write (*, *) 'Problem: NplgtoN23'
            end if
            !need the 1+ in above defn. because we're searcing 5 columns, but we want to start with column 2.
            !So if the first 1 is in column 2, NplgtoN23 would equal 1. We need to add 1 to it.
            Pnsmall(1) = prob_N00n/(prob_N02n + prob_N00n)
            Pnsmall(2) = prob_N02n/(prob_N02n + prob_N00n)
            call cumsum(Pnsmall, 2, cumPnsmall)
            rin = urcw1()

            do j = 1, 2
                if (rin <= cumPnsmall(j)) then      !find the doublet PLi will move to
                    temp = j
                    exit
                end if
            end do

            if (temp == 1) then
                state(NplgtoN23, connect_ind(newsite2)) = 5           !change the state to 5, N03
            elseif (temp == 2) then
                state(NplgtoN23, connect_ind(newsite2)) = 6           !change the state to 6, N23
            end if
            PLilocationt(1, loc2) = connect_ind(newsite2)               !location of current PLi
            PLilocationt(2, loc2) = NplgtoN23

        elseif (newdoublet2 == 3) then                   !if PLi is going to intact Nplg
            NplgtoN23 = 1                               !only place we have intact fibrin is in column 1
            Pismall(1) = prob_N00i/(prob_N02i + prob_N00i)
            Pismall(2) = prob_N02i/(prob_N02i + prob_N00i)
            call cumsum(Pismall, 2, cumPismall)
            rin = urcw1()

            do j = 1, 2
                if (rin <= cumPismall(j)) then         !find the doublet PLi will move to
                    temp = j
                    exit
                end if
            end do

            if (temp == 1) then
                state(NplgtoN23, connect_ind(newsite2)) = 5            !change the state to 5, N03

            elseif (temp == 2) then
                state(NplgtoN23, connect_ind(newsite2)) = 6            !change the state to 6, N23
            end if
            PLilocationt(1, loc2) = connect_ind(newsite2)                !location of current PLi
            PLilocationt(2, loc2) = NplgtoN23

        elseif (newdoublet2 == 4) then                  !if PLi is going to N03
            call findfirst(state(:, connect_ind(newsite2)), Ntot, 5, N03toN33)  !find the next 5 entry in state vector
            state(N03toN33, connect_ind(newsite2)) = 7                 !change the state to 7, N33
            PLilocationt(1, loc2) = connect_ind(newsite2)                !location of current PLi
            PLilocationt(2, loc2) = N03toN33

        elseif (newdoublet2 == 5) then                  !if PLi is going to degraded doublet
            call findfirst(state(:, connect_ind(newsite2)), Ntot, -1, degradetoPLi) !find the next -1 entry in state vector
            state(degradetoPLi, connect_ind(newsite2)) = 9        !change the state to 9, degraded doublet with PLi bound
            PLilocationt(1, loc2) = connect_ind(newsite2)           !location of current PLi
            PLilocationt(2, loc2) = degradetoPLi

        end if

        statetemp = state

    end subroutine move2pli

end program micromodel

module SCClientSubs

#if MPI
    use mpi
#endif
use, intrinsic          :: ISO_C_Binding
use, intrinsic          :: ISO_FORTRAN_ENV, only : stdout=>output_unit ! get output unit for stdout, see https://stackoverflow.com/questions/8508590/standard-input-and-output-units-in-fortran-90

implicit none
private

integer                         :: server_comm ! client/server communicator that will be output by MPI_COMM_CONNECT
integer                         :: nc
integer                         :: iT
integer                         :: nT
real                            :: SC_DT       ! timestep for the supercontroller (FAST.farm)
integer                         :: UseSC        
integer                         :: Mod_AmbWind
integer                         :: Nseeds, Iseed
character(256)                  :: MPIinit_filename
logical, parameter              :: verbose = .true.
logical, parameter              :: textout = .true.
integer                         :: output_unit_eff
character(255), save            :: dir_ctrl, dir_fast, dir_farm  
logical                         :: ReadSC = .true.      ! false: bypass reading SC_input.dat and force UseSC=0 (for single-turbine use)

public                          :: SC_MPI, TSC, get_TSC

! User defined type for super controller - deprecated. internal use in SCClientsSubs module only, initialized by SC_Init without going through calling routine (turbine controller)
type :: TSC
    integer                     :: useSC
    integer(C_SIZE_T)           :: iT                           ! Current turbine number, for super controller
    integer(C_SIZE_T)           :: nT                           ! Number of turbines, for super controller
    real(C_FLOAT)               :: SC_DT                        ! FAST.farm time step
    character(255)              :: dir_ctrl, dir_fast           ! directory of this routine (containing dll); root OpenFAST directory (normally two folders up)
end type TSC

contains

!-------------------------------------------------------------------------
! Main super controller routine
subroutine SC_MPI(status, avrSWAP, lfilename, SCinit_filename, ierror)
    integer, intent(inout)                  :: status
    real(C_FLOAT), intent(inout)            :: avrSWAP(:)                       ! The swap array, used to pass data to, and receive data from, the DLL controller.
    integer, intent(in)                     :: lfilename
    character(kind=c_char), intent(in)      :: SCinit_filename(lfilename)       ! The name of the parameter input .IN file
    integer                                 :: ierror
    ! logical, parameter                      :: powerramp=.true.

    ! logical                                 :: initflag
    ! integer                                 :: taskID, ntasks

    if (status==1) status=1+modulo(nint(avrSWAP(2)/avrSWAP(3))-1,nint(SC_DT/avrSWAP(3))) !-1: termination, 0: initialization, 1: one OpenFAST time step ahead FAST.farm time step increment, 2: FAST.farm time step increment, else: ineffective OpenFAST time steps in between FAST.farm steps
    ! if (floor(avrSWAP(2)/SC_DT)==0 .and. status==2) status=3 !At first FAST.farm time step, do not receive from server
    avrSWAP(1) = status ! Update status in avrSWAP for communication with server
    
    ! if (status /=0) print*, iT
    if (UseSC==0) then
        avrSWAP(13)=12e6 ! make sure the commanded power is higher than rated power
        ! if(powerramp) avrSWAP(13)=avrSWAP(13)-ABS(FLOOR(avrSWAP(2)/100.0)*1.0e6-avrSWAP(13))  ! Power ramp
        if(status/=0) goto 10
    endif

    ! Check rank and size, this should be 0 and 1 as we are not running multiple instances of the same program (using MPI_Exec.exe), but parallel threads with OpenMP. 
    ! call MPI_INITIALIZED(initflag,ierror)
    ! taskiD=0
    ! ntasks=0
    ! if(initflag)  then   
    !     call MPI_COMM_RANK(MPI_COMM_WORLD, taskID, ierror)
    !     call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierror)
    ! endif    
    ! print*, taskID, ntasks

    if(verbose .and. status/=0) write(output_unit_eff,*) 'Client Main: turbine nr:', iT, '     status:', status, 'timestep:', nint(avrSWAP(2)/avrSWAP(3))

#if MPI
    
    select case (status) 
    case(0) !Initialize
        
        call SC_init(lfilename, SCinit_filename, avrSWAP(3), ierror)
        ! print*, iT
        if (UseSC==0) goto 10
        nc=size(avrSWAP)
        server_comm=0
        if(iT==1) call MPIClient_init(ierror)

    case(1:) !Step
        ! IF(status==1) call MPIClient_send(avrSWAP,ierror) ! Send to server at each farm-level timestep
        ! IF(status==2) call MPIClient_receive(avrSWAP,ierror) ! Receive command from server at each farm-level timestep, with one turbine-level timestep delay to let the supercontroller gather and process info from all turbines
        ! IF(modulo(status,2)==1) call MPIClient_send(avrSWAP,ierror) ! Send to server at each second (turbine-level) timestep
        ! IF(modulo(status,2)==0) call MPIClient_receive(avrSWAP,ierror) ! Receive command from server at each second (turbine-level) timestep, with one timestep delay 
        call MPIClient_send(avrSWAP,ierror) ! Send to server at each (turbine-level) timestep
        call MPIClient_receive(avrSWAP,ierror) ! Receive from server at each (turbine-level) timestep
    case(-1) !Send and stop if actual 
        if(verbose) write(output_unit_eff,*)'Client main: termination turbine nr: ', iT     
        ! call MPIClient_receive(avrSWAP,ierror)!Receive message from server as if not terminating
        avrSWAP(1)=-1
        call MPIClient_send(avrSWAP,ierror) !Inform server that this instance is actually terminating 
        if(verbose) write(output_unit_eff,*)'Client main: termination command sent turbine nr: ', iT     
        call MPIClient_receive(avrSWAP,ierror) !Upon termination message, the server immediately answers back
        if(verbose) write(output_unit_eff,*)'Client main: termination command confirmed received turbine nr: ', iT, 'status: ',avrSWAP(1)
        if(verbose .and. textout) call flush(output_unit_eff)
        if (avrSWAP(1)==-2) call MPIClient_stop(ierror) !After all instances have confirmed termination, the server answers with a global termination flag for the last instance to close the door.
        call cleanup()
    case default !Cycle
    endselect

#else

    ReadSC=.false.
    call SC_init(lfilename, SCinit_filename, avrSWAP(3), ierror)
        
#endif

    10 continue
    if (Mod_AmbWind==1) call VTKSync(status,avrSWAP(2))
    ! print*, iT
    if (ierror/=0) call cleanup()
    
end subroutine SC_MPI
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
! Read SC input parameters
subroutine SC_init(lfilename,SCinit_filename_C,DT,ierror)

    integer, intent(in)                     :: lfilename
    character(kind=c_char), intent(in)      :: SCinit_filename_C(lfilename)         ! The name of the parameter input .IN file
    real(C_FLOAT), intent(in)               :: DT                                   !Turbine-level timestep, uses if ReadSC=false
    integer, intent(inout)                  :: ierror
    character(lfilename-1)                  :: SCinit_filename
    character(255)                          :: temp
    logical                                 :: file_exists
    integer                                 :: ind
    integer                                 :: nunit
    integer                                 :: e
    character(len=3), dimension(2)          :: status
    integer                                 :: count

    SCinit_filename = TRANSFER(SCinit_filename_C,SCinit_filename)

    ! Read turbine number from parent folder name
    ! call getcwd(dir_ctrl)
    ! dir_ctrl=trim(adjustl(dir_ctrl))//'\'//SCinit_filename
    dir_ctrl=SCinit_filename ! It looks like in newer versions of OpenFAST the absolute path is given
    dir_fast=dir_ctrl
    count=1
    do while(.true.)
        ind=max(index(dir_fast,"/",back=.true.),index(dir_fast,"\",back=.true.))
        temp = dir_fast(ind+1:) ! stripped directory name
        dir_fast = dir_fast(1:ind-1) ! path to root directory
        if (count==1) dir_ctrl=dir_fast ! control directory is the root directory of controller input file (given by SCinit_filename)
        if (count==2 .and. .not. ReadSC) then ! one turbine. assume simulation directory is root directory of control directory
            it=1
            exit 
        endif

        ! print*, 'dir_fast=',trim(adjustl(dir_fast)), 'temp=', trim(adjustl(temp))
        if(temp(1:1)=="T") then
            read(temp(2:),*,IOSTAT=e) iT
            if(e==0) exit
        endif
        
        if(trim(adjustl(temp))=="") then
            print*, "DISCON: Did not manage to identify turbine number. The control file location must fulfill the architecture <Farm root folder>/<Farm root folder>_<Seed number>/OpenFAST/T<Turbine number>/ControlData "
            stop
        endif
        count=count+1
    enddo

    output_unit_eff = stdout
    if (verbose .and. textout) then ! Print to text file instead of stdout
        output_unit_eff = 100 
        ! inquire( file=trim(adjustl(dir_ctrl))//"/stdout_SCClientSubs.txt", exist=file_exists)
        open( unit=output_unit_eff, file=trim(adjustl(dir_ctrl))//"/stdout_SCClientSubs.txt")
    endif
    
    ! print*, 'dir_fast=',trim(adjustl(dir_fast)), ', temp=', trim(adjustl(temp))

    ! Read from SC_input.dat, located in the same root directory as the turbine-specific T*it* folders
    if (ReadSC) then
        nunit=40
        open(unit=nunit,file=trim(adjustl(dir_fast))//"/SC_input.dat",action='read')
            read(nunit,*) UseSC    ! 0: do not use SC, 1: use SC
            read(nunit,*) nT       ! Number of turbines
            read(nunit,*) SC_DT    ! Farm-level timestep
            read(nunit,*) MPIinit_filename ! path of input file for MPI connection
            read(nunit,*) Mod_AmbWind ! Ambient wind mode in FAST.farm
            read(nunit,*) Nseeds  ! Number of parallel simulations
        close(nunit)
    else
        UseSC=0
        nT=1
        SC_DT=DT
        MPIinit_filename="dummy"
        Mod_AmbWind=3
        Nseeds=1
    endif

    ! Identify which seed we are
    if (Nseeds==1) then
        Iseed=1
    else 
        if (.NOT. ReadSC) then
            dir_farm = dir_fast
        else
            ind=max(index(dir_fast,"/",back=.true.),index(dir_fast,"\",back=.true.))
            if (dir_fast(ind+1:)/="OpenFAST") then
                print*, "DISCON: Did not manage to identify OpenFAST folder. The control file location must fulfill the architecture <Farm root folder>/<Farm root folder>_<Seed number>/OpenFAST/T<Turbine number>/ControlData "
                stop
            endif
            dir_farm = dir_fast(1:ind-1) ! path to  farm root directory
            ind=index(dir_farm,"_",back=.true.)
            read(dir_farm(ind+1:),*,IOSTAT=e) Iseed
            if (e/=0) then
                print*, "DISCON: Did not manage to identify seed number. The control file location must fulfill the architecture <Farm root folder>/<Farm root folder>_<Seed number>/OpenFAST/T<Turbine number>/ControlData "
                stop
            endif         
        endif
    endif


    ! print*, iT, dir_fast
    ! print*, iT, dir_ctrl
    ! print*, iT, dir_farm


end subroutine SC_init

!-------------------------------------------------------------------------
! Get farm-level variables
subroutine get_TSC(SC_var)
    type(TSC), intent(inout) :: SC_var
    SC_var%useSC = useSC
    SC_var%SC_DT = SC_DT
    SC_var%iT = iT
    SC_var%nT = nT
    SC_var%dir_ctrl = dir_ctrl
    SC_var%dir_fast = dir_fast
endsubroutine get_TSC

!-------------------------------------------------------------------------
! Works in coordination with InflowWind_Driver to progressively provide vtk files to FAST.farm
subroutine VTKSync(status,time)
    integer, intent(in)                 :: status
    real, intent(in)                    :: time
    real, save                          :: VTK_DT ! timesteps for high-resolution domain vtk files
    character(3), save                  :: tnumstr
    character(255)                      :: path_src, path_dst
    character(511)                      :: command
    logical                             :: file_exists
    integer                             :: count_clock, count_rate
    real, save                          :: t0, count_rate_f 
    real                                :: tic, toc, t_write, t_check, t_delete,  t_checklo, t_writelo, t_deletelo
    logical                             :: warningflag
    integer                             :: i, inext, iprev
    integer, parameter                  :: Out = 2 ! Where does IfWDrv write its files? 0: local vtk folder; 1: farm-level WindData folder of current simulation; 2: other (typically remote drive)
    integer                             :: error
    logical, parameter                  :: OutputTiming = .false.
    character(255), save                :: winddir_src_t, winddir_src_f, winddir_dst ! path to high-resolution domain source files, path to low-resolution domain source files and wind data directory for the simulation (WindFilePath in .fstf file)
    character(2)                        :: seednr_str
    character(*), parameter             :: RootOut= &              ! Root path name for wind files if Out=2
   "C:\Users\valentinc\Workspace\OpenFAST\Simulations\Demo_TCWPP_SEWTC_HPC\NREL5MW\10mps_Debug\RemoteDrive" 

    if (status==0) then
        write(tnumstr,'(i3)') iT
        if (OutputTiming) then
            call system_clock(count_clock,count_rate)
            count_rate_f=real(count_rate)
            t0 = count_clock/count_rate_f
        endif
        winddir_src_t=adjustl(trim(adjustl(dir_ctrl))//"\..\WindData")
        winddir_src_f=adjustl(trim(adjustl(dir_ctrl))//"\..\..\..\WindData")
        if(Out/=2) then
            winddir_dst=winddir_src_f
        else
            winddir_dst=adjustl(trim(RootOut)//"\WindData")
        endif
        if (Nseeds/=1 .or. Out==2) then
            write(seednr_str,'(i2)') Iseed
            winddir_dst= trim(adjustl(winddir_dst))//"_"//trim(adjustl(seednr_str))
        endif
    endif

    ! Read output timestep in InflowWind_driver's input file at init
    if (status==0) then
        inquire(file=trim(winddir_src_t)//"\IfWDrv.inp",exist=file_exists)
        if (file_exists) then
            open (unit=30,file=trim(winddir_src_t)//"\IfWDrv.inp",status='old')
            do i=1,12
                read(30,*)
            enddo
            read(30,*) VTK_DT
            close(30)
        else
            print*, 'VTKSync: the InflowWind driver input file was not found'
        endif
    endif

    ! Update time for InflowWind_driver
    if (iT==1) call UpdateTime(floor(time/SC_DT),trim(winddir_src_f),1000,t_writelo) ! Low-resolution domain
    call UpdateTime(nint(time/VTK_DT),trim(winddir_src_t),1000+iT,t_write) 
    if (iT==1 .and. OutputTiming) t_write=t_write+t_writelo
    
    ! Make sure vtk files for the next farm-level timestep are in FAST.farm's WindData folder
    inext=ceiling((time+VTK_DT)/SC_DT) ! Take a margin of VTK_DT to make sure time=0s gives inext=1
    if (iT==1) call CheckInput(inext,  trim(winddir_dst)//"\Low",t_checklo)
    call CheckInput((inext+1)*ceiling(SC_DT/VTK_DT)+2,  trim(winddir_dst)//"\HighT"//trim(adjustl(tnumstr)),t_check) ! take a margin of 2 turbine-level timesteps as the controller is called after reading vtk files in FAST.farm, and there might be missing timesteps once in a while due to different sampling rates
    if (iT==1 .and. OutputTiming) t_check=t_check+t_checklo

    if (Out==0) then
        command = "@copy /y "//trim(adjustl(path_src))//" "//trim(adjustl(path_dst))//" >NUL"
        ! print*, iT, status, command
        if (OutputTiming) call cpu_time(tic)
        call execute_command_line(command)
        if (OutputTiming) call cpu_time(toc)
        ! print*, iT, status, 'create',toc-tic
    endif

    
    ! Remove old vtk files
    ! Delete functionality is more efficient (less time critical) in InflowWind_driver.exe or in GarbageCollector.bat
    ! iprev=floor(time/SC_DT)
    ! if (iT==1) call DeleteOld(iprev,  trim(winddir_dst)//"\Low",33,t_deletelo)
    ! call DeleteOld(iprev*ceiling(SC_DT/VTK_DT),  trim(winddir_dst)//"\HighT"//trim(adjustl(tnumstr)),33,t_delete)
    ! if (iT==1 .and. OutputTiming) t_delete=t_delete+t_deletelo
    t_delete = 0.0
    

    ! Monitor computing times
    if (OutputTiming) then
        if (status==0) open (unit=34,file=trim(winddir_src_t)//"\Timing_VTKSync.dat",action='write',status='replace')
        call system_clock(count_clock)
        write(34,'(f13.4, 2i5, 3f13.4)') (count_clock/count_rate_f)-t0, nint(time/VTK_DT), (inext+1)*ceiling(SC_DT/VTK_DT)+2, t_check, t_delete, t_write 
        if (status<0) close(34)
    endif

    contains

    subroutine UpdateTime(tstpnr,FolderName,fid,tcmp)
        integer, intent(in)                 :: tstpnr
        character(*), intent(in)            :: FolderName
        integer, intent(in)                 :: fid
        real, intent(out)                   :: tcmp
        character(255)                      :: Filename


        if (status==0) then

            Filename=FolderName//"\Timestep.dat"
            open (unit=fid,file=Filename,action='write',status='replace')
            ! print*, Filename
        endif
        rewind(fid)
        if (status<0) then
            write(fid,'(i6)') 0 ! reset to 0 at end of simulation
            close(unit=fid)
            return
        else 
            if (OutputTiming) call cpu_time(tic)    
            do while(.true.)
                write(fid,'(i6)',iostat=error) tstpnr
                if(error==0) exit
                print*, "Error writing to file ",trim(adjustl(Filename))
            enddo
            if (OutputTiming) then
                call cpu_time(toc)
                tcmp=toc-tic
            else
                tcmp=0.0
            endif
        endif
        ! print*, 'UpdateTime: ',iT, time, tstpnr
    endsubroutine UpdateTime

    subroutine CheckInput(tstpnr,FolderName,tcmp)
        integer, intent(in)                 :: tstpnr
        character(*), intent(in)            :: FolderName
        real, intent(out)                   :: tcmp
        character(6)                        :: tstpstr

        write(tstpstr,'(i6)') tstpnr
        if (status==0) then
            if (Out==0) then
                command = "@rd /s /q "//FolderName//" >NUL"
                call execute_command_line(command)
                command = "@mkdir "//FolderName//" >NUL"
                ! print*, iT, status, command
                call execute_command_line(command)
            endif
        endif
        path_dst = FolderName//"\Amb.t"//trim(adjustl(tstpstr))//".vtk"
        if (Out==0) then
            path_src = trim(winddir_src_t)//"\vtk\DiskYZ.t"//trim(adjustl(tstpstr))//".vtk"
        else
            path_src = path_dst
        endif
        ! print*, status, path_src
        call system_clock(count_clock)
        tic=count_clock/count_rate_f
        warningflag=.false.
        do while (.true.)
            inquire(file=path_src,exist=file_exists)
            call system_clock(count_clock)
            tcmp=(count_clock/count_rate_f)-tic
            if (file_exists) exit
            if (tcmp>2.0 .and. .not. warningflag) then
                print*, "VTKSync: waiting for ",trim(adjustl(path_src))," to be created"
                warningflag=.true.
            endif
        enddo
    endsubroutine CheckInput

    subroutine DeleteOld(tstpnr,FolderName,fid,tcmp)
        integer, intent(in)                 :: tstpnr
        character(*), intent(in)            :: FolderName
        integer, intent(in)                 :: fid
        real, intent(out)                   :: tcmp
        character(6)                        :: tstpstr
        if (OutputTiming) call cpu_time(tic)
        write(tstpstr,'(i6)') tstpnr
        
        if (Out==0) then ! delete file in source directory
            path_src = trim(winddir_src_t)//"\vtk\DiskYZ.t"//trim(adjustl(tstpstr))//".vtk"
            inquire(file=path_src,exist=file_exists)
            ! command="@del "//trim(adjustl(path_src))//"/q >NUL"
            ! if (file_exists) call execute_command_line(command)
            if (file_exists) then
                open(unit=fid,file=path_src)
                error=1
                do while(error/=0)
                    close(unit=fid,status="delete",iostat=error)
                enddo
            endif
        endif
        path_dst = FolderName//"\Amb.t"//trim(adjustl(tstpstr))//".vtk"
        inquire(file=path_dst,exist=file_exists)
        ! command="@del "//trim(adjustl(path_dst))//"/q >NUL"
        ! if (file_exists) call execute_command_line(command)
        if (file_exists) then
            open(unit=fid,file=path_dst)
            error=1
            do while(error/=0)
                close(unit=fid,status="delete",iostat=error)
            enddo
        endif
        if (OutputTiming) then
            call cpu_time(toc)
            tcmp=toc-tic
        else
            tcmp=0.0
        endif
    endsubroutine DeleteOld

end subroutine VTKSync

#if MPI

!-------------------------------------------------------------------------
! Sets up MPI communication
subroutine MPIClient_init(ierror)
    integer ierror
    character(MPI_MAX_PORT_NAME) port_name
    character(6) FMT     
    logical initflag, file_exists
    character(255) current_dir
    integer MPI_ProvidedSupport

    call MPI_INITIALIZED(initflag,ierror)
    ! if(.not. initflag) call MPI_INIT(ierror) ! Init mpi protocol
    if(.not. initflag) call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,MPI_ProvidedSupport,ierror)
    if(verbose) write(output_unit_eff,*) 'MPI provided support:',MPI_ProvidedSupport,MPI_THREAD_MULTIPLE,MPI_THREAD_SERIALIZED
    inquire( file=trim(adjustl(MPIinit_filename)), exist=file_exists) 
    if (.not. file_exists) then !The path is not absolute. Assume the path is relative.
        call getcwd(current_dir)
        MPIinit_filename=trim(adjustl(current_dir))//'\'//trim(adjustl(MPIinit_filename))
        inquire( file=trim(adjustl(MPIinit_filename)), exist=file_exists) 
        if (.not. file_exists) write(output_unit_eff,*) "The path for shared MPI setup file provided in the controller input file *.IN for turbine ",iT," does not correspond to an existing absolute or relative path"
    endif
    open(1, file=trim(adjustl(MPIinit_filename)), status='old')
    write(FMT,'(i3)') MPI_MAX_PORT_NAME
    FMT = '(A' // FMT(1:3) // ')'
    read(1,FMT) port_name
    !print *, "port name is ", port_name
    port_name=adjustl(port_name) ! Remove leading blank

    call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, server_comm, ierror) ! Request connection to server

    if(verbose) write(output_unit_eff,*) "connection established with communicators: ", server_comm, " (client/server), ", MPI_COMM_WORLD, "(global); error: ", ierror
    write(1,'(i11)') server_comm ! Broadcast communication ID to textfile for other instances to read at next call
    write(1,'(i8)') nT ! Broadcast number of turbines to textfile
    write(1,'(i8)') nc ! Broadcast message length to textfile
    write(1,'(e13.4)') SC_DT
    close(1)
    
    if (ierror/=0) call cleanup()
endsubroutine MPIClient_init
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
! Sends message to server
subroutine MPIClient_send(avrSWAP,ierror)
    integer ierror
    integer msg_status(MPI_STATUS_SIZE)
    real(C_FLOAT), intent(inout)                  :: avrSWAP(:)

    ! Retrieve communication ID from textfile at second call (except for turbine 1)
    if (server_comm==0) then
        open(1, file=trim(adjustl(MPIinit_filename)), status='old')
        read(1,*)
        read(1,'(i11)') server_comm
        close(1)
    endif

    !print*, 'Client send: turbine nr:', iT,'   timestep nr:', nint(avrSWAP(2)/SC_DT), '   communicators: ', server_comm, ' (client/server), ', MPI_COMM_WORLD, '(global), error: ', ierror
    if(verbose) write(output_unit_eff,*) 'Client send: turbine nr:', iT,'   timestep nr:', floor((avrSWAP(2)-avrSWAP(3))/SC_DT)+1, '   status:', nint(avrSWAP(1))
    if(verbose .and. textout) call flush(output_unit_eff)
    call MPI_SEND(avrSWAP, nc, MPI_FLOAT, 0, iT, server_comm, ierror) ! C_FLOAT and MPI_FLOAT are equivalent
    if(verbose) write(output_unit_eff,*) 'Client send: turbine nr:', iT, ' sending complete' ! Warning: completed MPI_SEND does not mean the message is received on the other side, it may just have been buffered https://docs.microsoft.com/en-us/message-passing-interface/mpi-send-function 
    if(verbose .and. textout) call flush(output_unit_eff)
    if (ierror/=0) call cleanup()
endsubroutine MPIClient_send
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
! Receives message from server
subroutine MPIClient_receive(avrSWAP,ierror)
    integer ierror
    integer msg_status(MPI_STATUS_SIZE)
    real(C_FLOAT), intent(inout)                  :: avrSWAP(:)
    logical, parameter                            :: blocking = .true.
    integer(C_INT)                                :: request
    logical                                       :: received
    real(C_DOUBLE)                                :: init_time, max_time, wait_time
 
    ! max_time=1.0
    ! init_time = MPI_WTIME()
    ! wait_time=0.0

    !Receive message
    if(verbose) write(output_unit_eff,*) 'MPIClient_receive: waiting ', iT, avrSWAP(1)
    if(verbose .and. textout) call flush(output_unit_eff)
    if (blocking) then
        call MPI_RECV(avrSWAP, nc, MPI_FLOAT, MPI_ANY_SOURCE, iT, server_comm, msg_status, ierror) !Accepts only messages with tag iT (i.e for this turbine).
    else ! This non-blocking alternative implements a timeout. There is however no straightforward way to force re-sending from the server after timeout, so the message is lost. As multi-threaded blocking communication seems to work, this is inactual and only kept here in case.
        do while (.true.)
            call MPI_IRECV(avrSWAP, nc, MPI_FLOAT, MPI_ANY_SOURCE, iT, server_comm, request, ierror)
            call MPI_TEST(request,received,msg_status,ierror)
            wait_time = MPI_WTIME()-init_time
            if (received) then
                exit
            elseif (wait_time>max_time) then
                call MPI_CANCEL(request,ierror)
                exit
            endif
        enddo
    endif
    if(verbose) write(output_unit_eff,*) 'Client receive: turbine nr:', iT,'   timestep nr:', floor((avrSWAP(2)-avrSWAP(3))/SC_DT)+1, '   status:', nint(avrSWAP(1)), 'iT check:',msg_status(4)
    if(verbose) write(output_unit_eff,*) 'Client receive: turbine nr:', iT,'   msg :', avrSWAP
    if(verbose .and. textout) call flush(output_unit_eff)

    if (ierror/=0) call cleanup()
endsubroutine MPIClient_receive
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
! Terminates MPI communication
subroutine MPIClient_stop(ierror)
    integer ierror
    logical term_flag
    if(verbose) write(output_unit_eff,*) 'Client stop: terminating'
    if(verbose .and. textout) call flush(output_unit_eff)
    call MPI_COMM_DISCONNECT(server_comm, ierror) ! Close connection (done from the client's side)
    if(verbose) write(output_unit_eff,*) 'Client stop: disconnected'
    if(verbose .and. textout) call flush(output_unit_eff)
    call MPI_FINALIZED(term_flag,ierror)
    if (.not. term_flag) call MPI_FINALIZE(ierror) ! Terminate mpi protocol. Note: this should be done by the thread used to call mpi_init_thread, but using any thread (most likely a different one) seems to complete without error...
    if(verbose) write(output_unit_eff,*) 'Client stop: finished with error ', ierror 
    if(verbose .and. textout) call flush(output_unit_eff)

endsubroutine MPIClient_stop
!-------------------------------------------------------------------------

#endif

!-------------------------------------------------------------------------
! Terminates module 
subroutine cleanup()
    if (verbose .and. textout) close(output_unit_eff)
endsubroutine cleanup
!-------------------------------------------------------------------------

end module SCClientSubs
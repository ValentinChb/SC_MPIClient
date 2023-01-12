module SCClient

implicit none
    
contains
!**************************************************************************************************
subroutine SCClient_DISCON (avrSWAP, aviFAIL, avcINFILE, avcOUTNAME, avcMSG) bind(c,name='SCClient_DISCON')
    use, intrinsic :: ISO_C_Binding
    use SCClientSubs
#if DTUWEC
    use dtu_we_controller_bladed 
#endif   
#if ROSCO
    use rosco 
#endif

    

    !! DLLEXPORT attribute seems to be needed here in order for SCClient_DISCON to be included in dll when DTUWEC is statically linked (C binding is not sufficient, why???) 
    !GCC$ ATTRIBUTES DLLEXPORT :: SCClient_DISCON

    ! New fields in avrSWAP (currently unused by ServoDyn):
    ! 85: turbine number
    ! 86: length of control directory path
    ! 87: estimated wind speed
 
    ! Passed in variables from simulation codes (OpenFAST or Bladed):
    real(c_float),   intent(inout)    :: avrSWAP    (*)   ! The swap array, used to send data to, and receive data from, the DLL controller.
    integer(c_int),  intent(inout)    :: aviFAIL          ! A flag used to indicate the success of this DLL call set as follows:
                                                            !        = 0 if the DLL call was successful.
                                                            !        > 0 if the DLL call was successful but cMessage should be issued as a warning messsage.
                                                            !        < 0 if the DLL call was unsuccessful or for any other reason the simulation
                                                            !            is to be stopped at this point with cMessage as the error message.

    character(Kind=c_char), intent(in   )  :: avcINFILE (nint(avrSWAP(50)))   ! An array of 1-byte CHARACTERs giving the name of the parameter input file, 'DISCON.IN'.

    character(Kind=c_char), intent(inout)  :: avcMSG    (nint(avrSWAP(49)))   ! An array of 1-byte CHARACTERS giving the message contained in cMessage,
                                                                                ! which will be displayed by the calling program if aviFAIL <> 0.

    character(Kind=c_char), intent(inout)  :: avcOUTNAME(nint(avrSWAP(51)))   ! An array of 1-byte CHARACTERS giving the simulation run name
                                                                                ! including path without extension.

#if DTUWEC 
    ! DTUWEC variables
    real(c_float), dimension(87)           :: avrSWAP_DTUWEC
#endif   

#if ROSCO
    ! ROSCO variables
    real(c_float), dimension(87)           :: avrSWAP_ROSCO
    integer(c_size_t), save                :: ROSCO_IN_LEN
    character(Kind=c_char), save           :: ROSCO_IN(255) ! ROSCO_DISCON.IN//C_NULL_CHAR should be shorter that the controller input file's name (typically controller_input.dat) 
#endif

    ! Local variables
    real(c_float), dimension(87)           :: avrSWAP_SC
    integer(c_int)                         :: iStatus
    integer, save                          :: fidout
    logical                                :: PrintFlag
    type(TSC), save                        :: SC_var
    integer                                :: i
    character(Kind=c_char), save           :: ctrl_dir_c(512)
    integer(c_size_t), save                :: ctrl_dir_c_len
    integer                                :: error
    integer, parameter                     :: NavrSWAP=87
    real(c_float), parameter               :: pi = 3.14159265358979
    character(len=size(avcINFILE))         :: cInFile    
    

    iStatus = nint(avrSWAP(1))

    PrintFlag = mod(dble(avrSWAP(2)),SC_var%SC_DT)==0.0 .or. iStatus==0 ! Outputs only at new farm-level timestep
    if(avrSWAP(87)==0.0) avrSWAP(87) = avrSWAP(27)

    ! Call SuperController
    avrSWAP_SC = avrSWAP(1:87)
    call SC_MPI(iStatus, avrSWAP_SC, nint(avrSWAP_SC(50)), avcINFILE , aviFAIL) ! If MPI is not defined as preprocessor flag, this will only read input file at init and run hard-coded supercontroller is applicable
    ! Extract desired outputs
    if(iStatus>=1) avrSWAP(1:87)=avrSWAP_SC

    if (iStatus==0) then
        call get_TSC(SC_var)
        ctrl_dir_c=c_null_char
        ctrl_dir_c_len=len(trim(adjustl(SC_var%dir_ctrl)))+1 ! Add 1 for c_null_char at the end
        ctrl_dir_c(1:ctrl_dir_c_len)=transfer(trim(adjustl(SC_var%dir_ctrl))//c_null_char,ctrl_dir_c(1:ctrl_dir_c_len))
    endif

    avrSWAP(85)=REAL(SC_Var%iT)
    avrSWAP(86)=REAL(ctrl_dir_c_len)

    if (iStatus==0) then
        cInFile = transfer(avcINFILE(1:len(cInFile)),cInFile)
        i = index(cInFile, C_NULL_CHAR) - 1       ! Find the c NULL character at the end of cInFile, if it has then remove it
        if (i>0) cInFile = cInFile(1:i)
    endif

#if DTUWEC
    ! Call DTUWEC
    avrSWAP_DTUWEC = avrSWAP(1:NavrSWAP)
    call DTUWEC_DISCON (avrSWAP_DTUWEC, aviFAIL, avcINFILE, avcOUTNAME, avcMSG, ctrl_dir_c(1:ctrl_dir_c_len))
    ! Extract desired outputs
    avrSWAP(1:NavrSWAP)=avrSWAP_DTUWEC
#endif

#if ROSCO
    ! Call ROSCO
    if (iStatus==0) then
        i = max(index(cInFile,"/",back=.true.),index(cInFile,"\",back=.true.)) ! Get working directory as the root folder of the specified controller input file
        ROSCO_IN=c_null_char
        ROSCO_IN_LEN=i+len("ROSCO.IN")+1 ! Append a c_null_char at the end
        ROSCO_IN(1:ROSCO_IN_LEN)=transfer(ROSCO_IN(1:i)//"ROSCO.IN"//c_null_char,ROSCO_IN(1:ROSCO_IN_LEN))
    endif
    avrSWAP_ROSCO = avrSWAP(1:NavrSWAP)
    avrSWAP_ROSCO(50)=real(ROSCO_IN_LEN,c_float)
    avrSWAP_ROSCO(65) = merge(1.0,0.0,PrintFlag)
    call ROSCO_DISCON(avrSWAP_ROSCO,aviFAIL,ROSCO_IN(1:ROSCO_IN_LEN),avcOUTNAME,avcMSG)
    ! Extract desired outputs
    ! Vobs=avrSWAP_ROSCO(27)
#endif

    ! VC edit: print out avrSWAP for debug
    if(iStatus==0) then
        fidout=101
        open(unit=fidout,file=trim(SC_var%dir_ctrl)//'\SCClient_out.dat',iostat=error)
    endif
    write(fidout,'(*(e13.5))') avrSWAP(2), avrSWAP(21), avrSWAP(45)*180.0/pi, avrSWAP(4)*180.0/pi, avrSWAP(47), avrSWAP(13), avrSWAP(15), avrSWAP(27), avrSWAP(87), avrSWAP(48), yawangle_cmd, yawangle_meas
    if(PrintFlag) call flush(fidout)
    if (iStatus < 0) close(fidout)
    ! if(SC_var%iT==1) print*, 'DISCON: ',iStatus, avrSWAP(13), avrSWAP(87)

endsubroutine SCClient_DISCON

endmodule SCClient
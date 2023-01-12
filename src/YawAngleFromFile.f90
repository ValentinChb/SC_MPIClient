module YawAngleFromFile

implicit none

contains

subroutine YawAngleFromFileSub(t,iT,dir_ctrl,y_out)

    real, intent(in)                :: t
    integer, intent(in)             :: iT
    character(255), intent(in)      :: dir_ctrl
    real, intent(out)               :: y_out

    character(*), parameter :: file="YawAngles_10mps90deg.dat"
    integer, parameter      :: fid=303
    
    real, dimension(32)     :: y,y_prev
    real                    :: t_loc, t_loc_prev
    integer                 :: iT_loc
    real                    :: alpha, dt
    integer                 :: error


    open(unit=fid,file=trim(adjustl(dir_ctrl))//'/'//file,iostat=error)
    read(fid,*)

    t_loc=0.0
    y=0.0
    do while (t_loc<t)
        t_loc_prev=t_loc
        do iT_loc=1,32
            if(y(iT_loc)/=0.0) y_prev(iT_loc)=y(iT_loc)
        enddo
        read(fid,'(33E8.2)') t_loc, y
    enddo

    dt=t_loc-t_loc_prev
    alpha=(t_loc-t)/dt
    y_out=y(iT)*(1-alpha)+y_prev(iT)*alpha

    close(fid)

endsubroutine YawAngleFromFileSub

endmodule YawAngleFromFile
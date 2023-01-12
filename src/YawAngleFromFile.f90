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
    
    real, dimension(32)     :: y, y_next
    real                    :: t_loc
    integer                 :: iT_loc
    real                    :: alpha, dt
    integer                 :: error


    open(unit=fid,file=trim(adjustl(dir_ctrl))//'/'//file,iostat=error)
    read(fid,*)

    t_loc=0.0
    y=0.0
    y_next=0.0
    do while (t_loc<t)
        do iT_loc=1,32
            if(y_next(iT_loc)/=0.0) y(iT_loc)=y_next(iT_loc)
        enddo
        read(fid,*) t_loc, y_next
    enddo

    y_out=y(iT)

    close(fid)

endsubroutine YawAngleFromFileSub

endmodule YawAngleFromFile
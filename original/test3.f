      
      ! z = (f1,f1'), y = (f2)
      subroutine fsub(x,z,y,f)
            double precision x, z(2), y(1), f(2)
            f(1) = y(1) + z(2) ! f_1 = f1''(x)
            f(2) = y(1) + z(1)*z(2) - x ! f_2 = 0 = algebraic
      end
      subroutine dfsub(x,z,y,df)
            double precision x, z(2), y(0), df(2,3)
            df(1,1)=0 ! df_1/df1
            df(1,2)=1 ! df_1/df1'
            df(1,3)=1 ! df_1/df2

            df(2,1)=z(2) ! df_2/df1
            df(2,2)=z(1) ! df_2/df1'
            df(2,3)=1    ! df_2/df2
      end

      subroutine gsub(i,z,g)
            integer i
            double precision z(2), g
            go to (10,20) i
   10       g = z(1)+z(2)+0.15   ! g_1 = 0 = z(0)-1   left
            return
   20       g = z(1)-6.5 ! g_2 = 0 = z(0)-1   left
            return  
      end
      subroutine dgsub(i,z,dg)
            integer i
            double precision z(2), dg(2)
            go to (10,20) i
   10       dg(1)=1 ! dg_1/df1  
            dg(2)=1 ! dg_1/df1' 
            return
   20       dg(1)=1 ! dg_2/df1  
            dg(2)=0 ! dg_2/df1' 
            return   
      end

      program test
      implicit none
            integer ncomp, ny, orders(1),i,j,subsample
            integer ipar(12), ispace(100000), ltol(2), iflag
            double precision left,right,fspace(1000000),x,z(2),y(1)
            double precision sides(2), tol(2), fixpnt(0)
            external fsub,dfsub,gsub,dgsub
            double precision start, finish
            character(len = 32) :: outFormat = '(32(e19.12,1x))'

            ncomp=1;
            ny=1;
            orders(1) = 2;
            left = 0.0;
            right = 5.0;
            sides(1) = 0.0; ! first bc is on the left
            sides(2) = 5.0; ! second bc is on the left
            ltol(1) = 1 ! 1st tolerance is for z_1
            ltol(2) = 2 ! 2nd tolerance is for z_2
            tol(1) = 0.0001; ! 1st tolerance
            tol(2) = 0.0001; ! 2nd tolerance
            
            ipar(1)=1 ! problem is linear
            ipar(2)=0 ! auto num of collocation points
            ipar(3)=0 ! auto num of subintervals
            ipar(4)=1 ! number of tolerances
            ipar(5)=1000000 ! dim of fspace
            ipar(6)=100000 ! dim of ispace
            ipar(7)=1 ! -1 for full printout 
            ipar(8)=0 ! generate new mesh
            ipar(9)=0 ! no guess provided
            ipar(10)=0 ! problem is regular
            ipar(11)=0 ! no fixed points
            ipar(12)=0 ! dae index

           
            ! call cpu_time(start)
                         
            call coldae(ncomp,ny,orders,left,right,sides,
     .                 ipar,ltol,tol,fixpnt,ispace,fspace,iflag,
     .                 fsub,dfsub,gsub,dgsub,0)
           
!             call cpu_time(finish)
!             print '("Time = ",f14.9," seconds for 50000.")',
!      .       (finish-start)

            

            write(*,*) "iflag = ", iflag

            if(iflag == 1) then
                  open(100,file=trim("result3.txt"))
                  do i = 1, ispace(1)+1
                        if(i.eq.ispace(1)+1) then
                              x = fspace(i)
                              call appsln(x,z,y,fspace,ispace)
                              write(100,outFormat) x,z(1),z(2),y(1)
                        else 
                              subsample=1
                              do j = 0, subsample-1
                              x = fspace(i)*(1.0-float(j)/subsample)+
     .                                 fspace(i+1)*float(j)/subsample
                              call appsln(x,z,y,fspace,ispace)
                              write(100,outFormat) x,z(1),z(2),y(1)
                              enddo
                        endif
                  enddo
                  close(100)
            endif                                     

      end

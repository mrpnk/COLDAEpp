      ! z''(x)=sin(0.6*z'(x)) + x
      ! z(0)=1, z'(1)=-0.1   
      ! z_i = (z,z') 
      subroutine fsub(x,z,y,f)
            double precision x, z(2), y(0), f(2)
            f(1)=z(2) ! z'(x)
            f(2)=sin(0.6*z(2)) + x ! z''(x)
      end
      subroutine dfsub(x,z,y,df)
            double precision x, z(2), y(0), df(2,2)
            df(1,1)=0 ! dz'/dz
            df(1,2)=1 ! dz'/dz'
            df(2,1)=0 ! dz''/dz
            df(2,2)=0.6*cos(0.6*z(2)) ! dz''/dz'
      end

      subroutine gsub(i,z,g)
            integer i
            double precision z(2), g
            go to (10,20) i
   10       g=z(1)-1 ! 0=z(0)-1 left
            return
   20       g=z(2)+0.1 ! 0=z(0)-1 right
            return  
      end
      subroutine dgsub(i,z,dg)
            integer i
            double precision z(2), dg(2)
            go to (10,20) i
   10       dg(1)=1 ! d/dz  g_left
            dg(2)=0 ! d/dz' g_left
            return
   20       dg(1)=0 ! d/dz  g_right
            dg(2)=1 ! d/dz' g_right
            return   
      end

      program test
      implicit none
            integer ncomp, ny, orders(2),i,j
            integer ipar(12), ispace(10000), ltol(1), iflag
            double precision left,right,fspace(10000),x,z(2),y(0)
            double precision sides(2), tol(1), fixpnt(0)
            external fsub,dfsub,gsub,dgsub
            double precision start, finish
            character(len = 32) :: outFormat = '(32(e19.12,1x))'

            ncomp=2;
            ny=0;
            orders(1) = 1;
            orders(2) = 1;
            left = 0.0;
            right = 1.0;
            sides(1) = 0.0; ! first bc is on the left
            sides(2) = 1.0; ! second bc is on the right
            ltol(1) = 1 ! 1st tolerance is for z_1
            tol(1) = 0.0001; ! 1st tolerance

            ipar(1)=1 ! problem is linear
            ipar(2)=0 ! auto num of collocation points
            ipar(3)=0 ! auto num of subintervals
            ipar(4)=1 ! number of tolerances
            ipar(5)=10000 ! dim of fspace
            ipar(6)=10000 ! dim of ispace
            ipar(7)=-1 ! full printout 
            ipar(8)=0 ! generate new mesh
            ipar(9)=0 ! no guess provided
            ipar(10)=0 ! problem is regular
            ipar(11)=0 ! no fixed points
            ipar(12)=0 ! dae index, ignored

           
            call cpu_time(start)
                         
            call coldae(ncomp,ny,orders,left,right,sides,
     .                 ipar,ltol,tol,fixpnt,ispace,fspace,iflag,
     .                 fsub,dfsub,gsub,dgsub,0)

            call cpu_time(finish)
            print '("Time = ",f6.3," seconds.")',finish-start

            

            write(*,*) "iflag = ", iflag

            if(iflag == 1) then
                  open(100,file=trim("result.txt"))
                  do i = 1, ispace(1)+1 
                        x = fspace(i)
                        call appsln(x,z,y,fspace,ispace)
                        write(100,outFormat) x,z(1)
                  enddo
                  close(100)
            endif                                     

      end

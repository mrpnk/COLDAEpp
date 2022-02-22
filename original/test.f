      ! z'(x)=x
      ! z(0)=1    
      subroutine fsub(x,z,y,f)
            double precision x, z(1), y(0), f(1)
            f(1)=x ! z'(x)
      end
      subroutine dfsub(x,z,y,df)
            double precision x, z(1), y(0), df(1,1)
            df(1,1)=0 ! dz'/dz
      end

      subroutine gsub(i,z,g)
            integer i
            double precision z(1), g
            g=z(1)-1 ! 0=z(0)-1
      end
      subroutine dgsub(i,z,dg)
            integer i
            double precision z(1), dg(1)
            dg(1)=1 ! d/dz (z(0)-1)
      end

      subroutine evil(mat)
            integer mat(7,1)

            do i=1,3
                  do j=1,3
                        mat(i,j)=999
                  enddo
            enddo
      end

      subroutine fox(i)
            integer i

            i=2345
      end

      program test
      implicit none
            integer ncomp, ny, orders(1),i,j
            integer ipar(12), ispace(10000), ltol(1), iflag
            double precision left,right,fspace(10000),x,z(1),y(0)
            double precision sides(1), tol(1), fixpnt(0)
            external fsub,dfsub,gsub,dgsub

            character(len = 32) :: outFormat = '(32(e19.12,1x))'


            ! i=5
            ! write(*,*) "i=",i
            ! call fox(i+1)
            ! write(*,*) "i=",i

            ! integer matrix(28,7)
            ! do i=1,28
            !       do j=1,7
            !             matrix(i,j)=0
            !       enddo
            ! enddo
            ! call evil(matrix(1,1))
            ! do, i=1,28
            !       write(*,"(100g15.5)") ( matrix(i,j), j=1,7 )
            ! enddo

            ! return


            ncomp=1;
            ny=0;
            orders(1) = 1;
            left = 0.0;
            right = 1.0;
            sides(1) = 0.0; ! first bc is on the left
            ltol(1) = 1 ! 1st tolerance is for z_1
            tol(1) = 0.1; ! 1st tolerance

            ipar(1)=0 ! problem is linear
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

            call coldae(ncomp,ny,orders,left,right,sides,
     .                 ipar,ltol,tol,fixpnt,ispace,fspace,iflag,
     .                 fsub,dfsub,gsub,dgsub,0)

            write(*,*) "iflag = ", iflag

            if(iflag == 1) then
                  open(100,file=trim("result.txt"))
                  do i = 1, ispace(1)+1 
                        x = fspace(i)
                        call appsln(x,z,y,fspace,ispace)
                        write(100,outFormat) x,z(1)
                  enddo
                  close(100)

                  x=0.14
                  ! write(*,*)"vorh x=",x
                  ! call appsln(x,z,y,fspace,ispace)
                  ! write(*,*)"nach x=",x

            endif                                     

      end

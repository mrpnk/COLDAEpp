program test

	implicit none

	integer i(2,2), j(1)

	i(1,1) = 1; i(1,2) = 3;
	i(2,1) = 2; i(2,2) = 4;

	j(1) = 77

	print*, i
	call func(j)
	print*, i

	end program test





	subroutine func(arg)
	integer arg(1,2)

	print*, arg
	arg(1,1)=235;



	if(3.eq.3)	go to 10
	print*, "nein"
10 print*, "ja?"
	print*, "ja"

end



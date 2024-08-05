module qsort1
use initialize
use particle
implicit none


contains

recursive subroutine QSort(A,nA)

    ! DUMMY ARGUMENTS
    integer, intent(in) :: nA
    type (group), dimension(nA), intent(in out) :: A

    ! LOCAL VARIABLES
    integer :: left, right
    real :: random
    integer :: pivot
    type (group) :: temp
    integer :: marker

    if (nA > 1) then

        call random_number(random)
        pivot = A(int(random*real(nA-1))+1)%value   ! random pivor (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1

        do while (left < right)
            right = right - 1
            do while (A(right)%value > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left)%value < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do

        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
        call QSort(A(:marker-1),marker-1)
        call QSort(A(marker:),nA-marker+1)

    end if
end subroutine QSort


subroutine sort(array1,array2,len) 
    use initialize
    use particle
    implicit none
    integer,intent(in) :: len
    real(dp),intent(inout) :: array1(1:len)
    integer,intent(inout) :: array2(1:len)
    
    integer :: is
    type(group) :: Alist(1:len)
    real(dp) :: temp(1:len)



    do is=1,len 
        Alist(is)%order=is
        Alist(is)%value=array2(is)
    end do

    call QSort(Alist,len)

    do is=1,len 
        array2(is)=Alist(is)%value
        temp(is)=array1(Alist(is)%order)
    end do

    array1(1:len)=temp(1:len)


    
end subroutine sort

end module qsort1
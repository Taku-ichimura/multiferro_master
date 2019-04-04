module iniset

  use global
  use mtmod

  implicit none
contains

  subroutine iniJ_ij(J_ij,seed)

    integer,intent(in) :: seed
    real(8),intent(inout) :: J_ij(:,:,:,:)
    integer i,j,k,m,count
    real(8) r1 ,r2

    call sgrnd(seed)  ! 乱数の種(seed)を指定して、ウォームアップ(warm)
    count = 1

    do i = 1,L
       do j = 1,L
          do k=1,L
             do m = 1,N_bond
                
                if(mod(count,2) /= 0)then
                   r1 = grnd()
                   r2 = grnd()
                   ! Box-Muller transformation
                   if(r2<=0.d0) then
                      
                      r2=grnd()
                   end if

                   r1=J_0*sqrt(-2.0d0*log(r2))*cos(twopi*r1)+ave
                   r2=J_0*sqrt(-2.0d0*log(r2))*sin(twopi*r1)+ave

                   J_ij(m,k,j,i)=r1
                   count = count +1
                  
                else
                   J_ij(m,k,j,i)=r2
                   count = count +1
                                       

                end if
             end do
          end do
       end do
    end do

  end subroutine iniJ_ij

  subroutine iniJ_ij_niti(J_ij,seed)

    integer,intent(in) :: seed
    real(8),intent(inout) :: J_ij(:,:,:,:)
    integer i,j,k,m
    real(8) r1 ,r2



    call sgrnd(seed)  ! 乱数の種(seed)を指定して、ウォームアップ(warm)
    

    do i = 1,L
       do j = 1,L
          do k=1,L
             do m = 1,N_bond
                
                   
                   r1 = grnd()
                   if(r1 <= Niratio)then
                      J_ij(m,k,j,i) = -1.0d0
                   else
                      J_ij(m,k,j,i) = 1.0d0
                   end if                          

                
             end do
          end do
       end do
    end do
  end subroutine iniJ_ij_niti
  

  

  subroutine initiarize_spin(spinA,spinB,seed)

    integer,intent(in) :: seed
    real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer i,j,k,T_i
    real(8) r1,r2,r3,r4,tmp1,tmp2

    call sgrnd(seed)  ! 乱数の種(seed)を指定して、ウォームアップ(warm)

    do T_i = 1,N_T

       do i = 1,L
          do j = 1,L
             do k = 1,L

                r1 = grnd()*2.0d0-1.0d0
                r2 = grnd()*twopi
                r3 = grnd()*2.0d0-1.0d0
                r4 = grnd()*twopi

                tmp1 = sqrt(1.0d0 - r1*r1)
                tmp2 = sqrt(1.0d0 - r3*r3)

                spinA(1,k,j,i,T_i) = tmp1*cos(r2)
                spinA(2,k,j,i,T_i) = tmp1*sin(r2)
                spinA(3,k,j,i,T_i) = r1
                spinB(1,k,j,i,T_i) = tmp2*cos(r4)
                spinB(2,k,j,i,T_i) = tmp2*sin(r4)
                spinB(3,k,j,i,T_i) = r3

             end do
          end do
       end do

    end do

  end subroutine initiarize_spin

   subroutine initiarize_spin_rep(spin_repA,spin_repB)

    real(8),intent(inout) :: spin_repA(:,:,:,:,:),spin_repB(:,:,:,:,:)
    integer i,j,k,T_i
    real(8) r1,r2

     do T_i = 1,N_T
        
       do i = 1,L
          do j = 1,L
             do k = 1,L

                
                r1 = grnd()*twopi
                r2 = grnd()*twopi

                spin_repA(1,k,j,i,T_i) = cos(r1)
                spin_repA(2,k,j,i,T_i) = sin(r1)
                spin_repB(1,k,j,i,T_i) = cos(r2)
                spin_repB(2,k,j,i,T_i) = sin(r2)

             end do
          end do
       end do

    end do

  end subroutine initiarize_spin_rep

  
  subroutine stock_random(randnum)

    real(8),intent(inout) :: randnum(:,:)
    integer T_i,i,j,k
    real(8) r1,r2,r3,r4

    do T_i = 1,N_T

       do i = 1,N_spin*2

          randnum(i,T_i) = grnd()

       end do
    end do

   


  end subroutine stock_random


  subroutine initiarize_temperature(T)

    real(8),intent(inout) :: T(:)
    integer T_i
    real(8) dT

    dT = (T_max - T_min)

    do T_i =1,N_T

       T(T_i) = dT*(T_i-1) + T_min

    end do

  end subroutine initiarize_temperature

  subroutine initiarize_T(T_ex,T_ex_rep)
    integer,intent(inout) :: T_ex(:),T_ex_rep(:)
    integer i

    do i = 1,N_T
       T_ex(i) = i
       T_ex_rep(i) = i
    end do

  end subroutine initiarize_T
  
  

end module iniset

    

                
                





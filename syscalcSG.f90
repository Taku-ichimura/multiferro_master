module syscalc

  use global


  implicit none
contains

  subroutine sysE(rg,spinA,spinB,T_i,T_ex,J_ij,E)
    real(8) ,intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:),J_ij(:,:,:,:)
    integer ,intent(in) ::T_i,T_ex(:)
    real(8) ,intent(in) :: rg(:)
    real(8) ,intent(out) :: E(:)
    real(8) local(3),kitaev(3),kitvec(3,3)
    integer i, j, k , size1 , size2 , size3
    integer a,b,c,d
    real(8)posi(3),posi2(3),base(4,3),local2(3)

    size1 = L
    size2 = L
    size3 = L
    E(T_i) = 0.0d0


    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    !**********bond vector*************
    kitvec(1,1) = 1.0d0
    kitvec(1,2) = 0.0d0
    kitvec(1,3) = 0.0d0
    kitvec(2,1) = 0.0d0
    kitvec(2,2) = 1.0d0
    kitvec(2,3) = 0.0d0
    kitvec(3,1) = 0.0d0
    kitvec(3,2) = 0.0d0
    kitvec(3,3) = 1.0d0
    !*********************************


    !*****************calc Energy*************************

    do i = 1, L

       !****************boundary condition*******************
       a = mod(i-1+size1-1,size1)+1
       !******************************************************

       do j = 1,L

          !*************BC*************
          d = mod(j-1+size2-1,size2)+1
          !*****************************

          do k =1,L

             !*************BC*************
             c = mod(k-1+size3-1,size3)+1
             !*****************************

             !***********position**************************************
             posi = (i-1)*base(1,:) +(j-1)*base(2,:) + (k-1)*base(4,:)
             posi = posi - rg
             posi = posi/dble(L)
             posi2 = (i-1)*base(1,:) +(j-1)*base(2,:) + base(3,:)+ (k-1)*base(4,:)
             posi2 = posi2 -rg
             posi2 = posi2/dble(L)
             !********************************************************

             kitaev = kitvec(1,:)*dot_product(spinB(:,k,j,i,T_ex(T_i)),kitvec(1,:)) &
                  + kitvec(2,:)*dot_product(spinB(:,k,d,i,T_ex(T_i)),kitvec(2,:))&
                  + kitvec(3,:)*dot_product(spinB(:,k,j,a,T_ex(T_i)),kitvec(3,:))

             local(:) = J_ij(3,k,j,i)*spinB(:,k,j,a,T_ex(T_i)) + J_ij(1,k,j,i)*spinB(:,k,j,i,T_ex(T_i))+&
                  J_ij(2,k,j,i)*spinB(:,k,d,i,T_ex(T_i)) + Jp*J_0*spinB(:,c,j,i,T_ex(T_i))

             local(1) = local(1)+XY_0*spinA(2,k,j,i,T_ex(T_i))
             local(2) = local(2)+XY_0*spinA(3,k,j,i,T_ex(T_i))
             local(3) = local(3)+XY_0*spinA(1,k,j,i,T_ex(T_i))

             local2(1) = local2(1)+XY_0*spinB(2,k,j,i,T_ex(T_i))
             local2(2) = local2(2)+XY_0*spinB(3,k,j,i,T_ex(T_i))
             local2(3) = local2(3)+XY_0*spinB(1,k,j,i,T_ex(T_i))

             local(1) = local(1)+(spinB(3,k,d,i,T_ex(T_i))+spinB(2,k,j,a,T_ex(T_i)))*gam
             local(2) = local(2)+(spinB(3,k,j,i,T_ex(T_i))+spinB(1,k,j,a,T_ex(T_i)))*gam
             local(3) = local(3)+(spinB(2,k,j,i,T_ex(T_i))+spinB(1,k,d,i,T_ex(T_i)))*gam


             local = local + kitaev*K_0 -h0*H

             local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
             local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

             local2(1)=local2(1)-(Ht(1)*H(2)*posi2(2)-Ht(2)*H(1)*posi2(2))
             local2(2)=local2(2)-(-Ht(1)*H(2)*posi2(1)+Ht(2)*H(1)*posi2(1))
             local2(3) = local2(3)+0.0d0

             E(T_i) = E(T_i) + dot_product(spinA(:,k,j,i,T_ex(T_i)),local(:)) -&
                  dot_product(spinB(:,k,j,i,T_ex(T_i)),h0*H)+dot_product(spinB(:,k,j,i,T_ex(T_i)),local2)

          end do

       end do

    end do
    !******************end calc  E*****************************

  end subroutine sysE


  subroutine local_EA(rg,spinA,spinB,J_ij,T_i,i,j,k,dE)!local subA E

    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:),J_ij(:,:,:,:)
    integer,intent(in) ::  i, j, k ,T_i
    real(8),intent(out) :: dE
    real(8), intent(in) :: rg(:)
    real(8) local(3),kitvec(3,3),kitaev(3)
    real(8)base(4,3),posi(3),tmp(2),tmx,tmy,tmz,P(3)
    integer size1 , size2 , size3 ,a,b,c,d

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    !**********bond vector*************
    kitvec(1,1) = 1.0d0
    kitvec(1,2) = 0.0d0
    kitvec(1,3) = 0.0d0
    kitvec(2,1) = 0.0d0
    kitvec(2,2) = 1.0d0
    kitvec(2,3) = 0.0d0
    kitvec(3,1) = 0.0d0
    kitvec(3,2) = 0.0d0
    kitvec(3,3) = 1.0d0
    !*********************************

    !*****************BC**********************
    a = mod(i-1+size1-1,size1)+1
    c = mod(k-1+size3-1,size3)+1
    d = mod(j-1+size2-1,size2)+1
    !******************************************


    !***********troidal**************************************
    posi = (i-1)*base(1,:) +(j-1)*base(2,:) + (k-1)*base(4,:)
    posi = posi -rg
    posi = posi/dble(L)
    !********************************************************

    kitaev = kitvec(1,:)*dot_product(spinB(:,k,j,i,T_i),kitvec(1,:)) &
         + kitvec(2,:)*dot_product(spinB(:,k,d,i,T_i),kitvec(2,:))&
         + kitvec(3,:)*dot_product(spinB(:,k,j,a,T_i),kitvec(3,:))

    local(:) = J_ij(3,k,j,i)*spinB(:,k,j,a,T_i) + J_ij(1,k,j,i)*spinB(:,k,j,i,T_i) +&
         J_ij(2,k,j,i)*spinB(:,k,d,i,T_i) + Jp*J_0*spinB(:,c,j,i,T_i)

    local(1) = local(1)+XY_0*spinA(2,k,j,i,T_i)
    local(2) = local(2)+XY_0*spinA(3,k,j,i,T_i)
    local(3) = local(3)+XY_0*spinA(1,k,j,i,T_i)

    local(1) = (spinB(3,k,d,i,T_i)+spinB(2,k,j,a,T_i))*gam
    local(2) = (spinB(3,k,j,i,T_i)+spinB(1,k,j,a,T_i))*gam
    local(3) = (spinB(2,k,j,i,T_i)+spinB(1,k,d,i,T_i))*gam

    local = local + kitaev*K_0 - h0*H

    local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
    local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

    dE = dot_product(spinA(:,k,j,i,T_i),local(:))

  end subroutine local_EA



  subroutine local_EB(rg,spinB,spinA,J_ij,T_i,i,j,k,dE)!local E

    real(8),intent(in) :: spinB(:,:,:,:,:),spinA(:,:,:,:,:),J_ij(:,:,:,:)
    integer,intent(in) ::  i, j, k ,T_i
    real(8),intent(out) :: dE
    real(8),intent(in) ::rg(:)
    real(8) local(3),kitvec(3,3),kitaev(3)
    real(8)base(4,3),posi(3),tmp(2),tmx,tmy,tmz,P(3)
    integer size1 , size2 , size3 ,a,b,c,d

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    !**********bond vector*************
    kitvec(1,1) = 1.0d0
    kitvec(1,2) = 0.0d0
    kitvec(1,3) = 0.0d0
    kitvec(2,1) = 0.0d0
    kitvec(2,2) = 1.0d0
    kitvec(2,3) = 0.0d0
    kitvec(3,1) = 0.0d0
    kitvec(3,2) = 0.0d0
    kitvec(3,3) = 1.0d0
    !*********************************

    !*****************BC**********************
    a = mod(i+1+size1-1,size1)+1
    c = mod(k+1+size3-1,size3)+1
    d = mod(j+1+size2-1,size2)+1
    !******************************************


    !***********troidal**************************************
    posi = (i-1)*base(1,:) +(j-1)*base(2,:) + base(3,:) + (k-1)*base(4,:)
    posi = posi -rg
    posi = posi/dble(L)
    !********************************************************

    kitaev = kitvec(1,:)*dot_product(spinA(:,k,j,i,T_i),kitvec(1,:)) &
         + kitvec(2,:)*dot_product(spinA(:,k,d,i,T_i),kitvec(2,:))&
         + kitvec(3,:)*dot_product(spinA(:,k,j,a,T_i),kitvec(3,:))

    local(:) = J_ij(3,k,j,a)*spinA(:,k,j,a,T_i) + J_ij(1,k,j,i)*spinA(:,k,j,i,T_i) +&
         J_ij(2,k,d,i)*spinA(:,k,d,i,T_i) + Jp*J_0*spinA(:,c,j,i,T_i)

    local(1) = local(1)+XY_0*spinB(2,k,j,i,T_i)
    local(2) = local(2)+XY_0*spinB(3,k,j,i,T_i)
    local(3) = local(3)+XY_0*spinB(1,k,j,i,T_i)

    local(1) = (spinA(3,k,d,i,T_i)+spinA(2,k,j,a,T_i))*gam
    local(2) = (spinA(3,k,j,i,T_i)+spinA(1,k,j,a,T_i))*gam
    local(3) = (spinA(2,k,j,i,T_i)+spinA(1,k,d,i,T_i))*gam

    local = local + kitaev*K_0 -h0*H

    local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
    local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

    dE = dot_product(spinB(:,k,j,i,T_i),local(:))

  end subroutine local_EB

  subroutine calc_vec(spinA,spinB,T_i,T_ex,vec,Svec)!sum spin-vec
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:),Svec(:,:)
    real(8) rsA(3),rsB(3)
    integer i, j, k,size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L

    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

            call realizespin(spinA(:,k,j,i,T_ex(T_i)),rsA)
            call realizespin(spinB(:,k,j,i,T_ex(T_i)),rsB)
            vec(T_i,:) = vec(T_i,:) + (rsA+rsB)
            Svec(T_i,:) = Svec(T_i,:) + (spinA(:,k,j,i,T_ex(T_i)) + spinB(:,k,j,i,T_ex(T_i)))

          end do
       end do
    end do

    vec(T_i,:) = vec(T_i,:) /dble(N_spin)
    Svec(T_i,:) = Svec(T_i,:) /dble(N_spin)

  end subroutine calc_vec


  subroutine plane_calc_vec(spinA,spinB,T_i,T_ex,vec,Svec)!sum spin-vec
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:),Svec(:,:)
    real(8) rsA(3),rsB(3)
    integer i, j, k,size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L

    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

            call realizespin(spinA(:,k,j,i,T_ex(T_i)),rsA)
            call realizespin(spinB(:,k,j,i,T_ex(T_i)),rsB)
            vec(T_i,:) = vec(T_i,:) + ((-1)**k)*(rsA+rsB)
            Svec(T_i,:) = Svec(T_i,:) + ((-1)**k)*(spinA(:,k,j,i,T_ex(T_i)) + spinB(:,k,j,i,T_ex(T_i)))

          end do
       end do
    end do

    vec(T_i,:) = vec(T_i,:) /dble(N_spin)
    Svec(T_i,:) = Svec(T_i,:) /dble(N_spin)

  end subroutine plane_calc_vec


  subroutine staggered(spinA,spinB,T_i,T_ex,vec,Svec)!スタッガード磁化を返す
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:),Svec(:,:)
    real(8) rsA(3),rsB(3)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

            call realizespin(spinA(:,k,j,i,T_ex(T_i)),rsA)
            call realizespin(spinB(:,k,j,i,T_ex(T_i)),rsB)
            vec(T_i,:) = vec(T_i,:) -  rsA + rsB
            Svec(T_i,:) = Svec(T_i,:) -  spinA(:,k,j,i,T_ex(T_i)) + spinB(:,k,j,i,T_ex(T_i))

          end do
       end do
    end do

    vec(T_i,:) = vec(T_i,:) /dble(N_spin)
    Svec(T_i,:) = Svec(T_i,:) /dble(N_spin)

  end subroutine staggered

  subroutine zigzag(spinA,spinB,T_i,T_ex,vec,Svec)!zigzag order parameter
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:,:),Svec(:,:,:)
    real(8) rsA(3),rsB(3)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

            call realizespin(spinA(:,k,j,i,T_ex(T_i)),rsA)
            call realizespin(spinB(:,k,j,i,T_ex(T_i)),rsB)

            vec(T_i,:,1) = vec(T_i,:,1) +(1.0d0)**k*((-1.0d0)**i*rsA+(-1.0d0)**i*rsB)
            vec(T_i,:,2) = vec(T_i,:,2) +(1.0d0)**k*((-1.0d0)**j*rsA+(-1.0d0)**j*rsB)
            vec(T_i,:,3) = vec(T_i,:,3) +(-1.0d0)**k*((-1.0d0)**j*(-1.0d0)**i*rsA-(-1.0d0)**j*(-1.0d0)**i*rsB)

             Svec(T_i,:,1) = Svec(T_i,:,1) +(1.0d0)**k*((-1.0d0)**i*spinA(:,k,j,i,T_ex(T_i)) &
                  +(-1.0d0)**i* spinB(:,k,j,i,T_ex(T_i)))
             Svec(T_i,:,2) = Svec(T_i,:,2) +(1.0d0)**k*((-1.0d0)**j*  spinA(:,k,j,i,T_ex(T_i)) &
                  +(-1.0d0)**j* spinB(:,k,j,i,T_ex(T_i)))
             Svec(T_i,:,3) = Svec(T_i,:,3) +(-1.0d0)**k*((-1.0d0)**j*(-1.0d0)**i*  spinA(:,k,j,i,T_ex(T_i))&
                  -(-1.0d0)**j*(-1.0d0)**i* spinB(:,k,j,i,T_ex(T_i)))

          end do
       end do
    end do

    vec(T_i,:,:) = vec(T_i,:,:) /dble(N_spin)
    Svec(T_i,:,:) = Svec(T_i,:,:) /dble(N_spin)

  end subroutine zigzag

  subroutine zig2(spinA,spinB,T_i,T_ex,vec)!zigzag order parameter
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

             vec(T_i,:) = vec(T_i,:) +(-1.0d0)**j*(-1.0d0)**i*  spinA(:,k,j,i,T_ex(T_i)) &
                  +(-1.0d0)**j*(-1.0d0)**i* spinB(:,k,j,i,T_ex(T_i))

          end do
       end do
    end do

    vec(T_i,:) = vec(T_i,:) /dble(N_spin)

  end subroutine zig2

  subroutine zig3(spinA,spinB,T_i,T_ex,vec)!zigzag order parameter
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

             vec(T_i,:) = vec(T_i,:) +(-1.0d0)**j*(-1.0d0)**i*  spinA(:,k,j,i,T_ex(T_i))&
                  -(-1.0d0)**j*(-1.0d0)**i* spinB(:,k,j,i,T_ex(T_i))

          end do
       end do
    end do

    vec(T_i,:) = vec(T_i,:) /dble(N_spin)

  end subroutine zig3

  subroutine stripy(spinA,spinB,T_i,T_ex,vec,Svec)!zigzag order parameter
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:,:),Svec(:,:,:)
    real(8) rsA(3),rsB(3)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

            call realizespin(spinA(:,k,j,i,T_ex(T_i)),rsA)
            call realizespin(spinB(:,k,j,i,T_ex(T_i)),rsB)

            vec(T_i,:,1) = vec(T_i,:,1) +(1.0d0)**k*((1.0d0)**i*rsA-(-1.0d0)**i*rsB)
            vec(T_i,:,2) = vec(T_i,:,2) +(-1.0d0)**k*((-1.0d0)**j*((-1.0d0)**i*rsA+(-1.0d0)**i*rsB))
            vec(T_i,:,3) = vec(T_i,:,3) +(1.0d0)**k*((1.0d0)**j*rsA-(-1.0d0)**j*rsB)

            Svec(T_i,:,1) = Svec(T_i,:,1) +(1.0d0)**k*((1.0d0)**i*spinA(:,k,j,i,T_ex(T_i))-(-1.0d0)**i*spinB(:,k,j,i,T_ex(T_i)))
            Svec(T_i,:,2) = Svec(T_i,:,2) +(-1.0d0)**k*((-1.0d0)**j*((-1.0d0)**i*spinA(:,k,j,i,T_ex(T_i))+(-1.0d0)**i*spinB(:,k,j,i,T_ex(T_i))))
            Svec(T_i,:,3) = Svec(T_i,:,3) +(1.0d0)**k*((1.0d0)**j*spinA(:,k,j,i,T_ex(T_i))-(-1.0d0)**j*spinB(:,k,j,i,T_ex(T_i)))

          end do
       end do
    end do

    vec(T_i,:,:) = vec(T_i,:,:) /dble(N_spin)
    Svec(T_i,:,:) = Svec(T_i,:,:) /dble(N_spin)

  end subroutine stripy

  subroutine stripe2(spinA,spinB,T_i,T_ex,vec)!zigzag order parameter
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

             vec(T_i,:) = vec(T_i,:) +(-1.0d0)**j*((-1.0d0)**i*spinA(:,k,j,i,T_ex(T_i)) +(-1.0d0)**i*spinB(:,k,j,i,T_ex(T_i)))

          end do
       end do
    end do

    vec(T_i,:) = vec(T_i,:) /dble(N_spin)

  end subroutine stripe2

  subroutine stripe3(spinA,spinB,T_i,T_ex,vec)!zigzag order parameter
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

             vec(T_i,:) = vec(T_i,:) +(-1.0d0)**j*spinA(:,k,j,i,T_ex(T_i)) -(-1.0d0)**j*spinB(:,k,j,i,T_ex(T_i))

          end do
       end do
    end do

    vec(T_i,:) = vec(T_i,:) /dble(N_spin)

  end subroutine stripe3

  subroutine realizespin(spin,spinr)
    real(8),intent(in) :: spin(3)
    real(8),intent(out) :: spinr(3)
    real(8) vec(3,3)

    !***********abc axis*************
    vec(1,:) = (/1.0d0,1.0d0,-2.0d0/)
    vec(2,:) = (/-1.0d0,1.0d0,0.0d0/)
    vec(3,:) = (/1.0d0,1.0d0,1.0d0/)
    vec(1,:) = vec(1,:)/sqrt(6.0d0)
    vec(2,:) = vec(2,:)/sqrt(2.0d0)
    vec(3,:) = vec(3,:)/sqrt(3.0d0)
    !********************************

    spinr(1) = dot_product(spin,vec(1,:))
    spinr(2) = dot_product(spin,vec(2,:))
    spinr(3) = dot_product(spin,vec(3,:))

  end subroutine realizespin

  subroutine calc_chila_pola(spinA,spinB,T_i,T_ex,Chi,P,SChi,SP)!カイラリティー,分極の計算,(崎山さんのものに従いsubA×subBの外積、論文とABが逆なので注意)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: Chi(N_T,4,3),P(N_T,3),SChi(N_T,4,3),SP(N_T,3)
    real(8) tmp1(3),tmp2(3),tmp3(3),tmp4(3)
    real(8)tmpr1(3),tmpr2(3),tmpr3(3),tmpr4(3),spinrA(3,4),spinrB(3)
    real(8) bondvec(4,3),rbv(4,3)
    integer i, j, k,size1 , size2 , size3,pm,a,b,c,d

    size1 = L
    size2 = L
    size3 = L

    !*************bond vec***********
    bondvec(1,1)=-sqrt(3.0d0)*0.5
    bondvec(1,2)=-0.5d0
    bondvec(1,3)=0
    bondvec(2,1)=sqrt(3.0d0)*0.5
    bondvec(2,2)=-0.5d0
    bondvec(2,3)=0.0d0
    bondvec(3,1)=0.0d0
    bondvec(3,2)=1.0d0
    bondvec(3,3)=0.0d0
    bondvec(4,1)=0.0d0
    bondvec(4,2)=0.0d0
    bondvec(4,3)=1.0d0
    !********************************
    !************realize bond vec**********
    call realizespin(bondvec(1,:),rbv(1,:))
    call realizespin(bondvec(2,:),rbv(2,:))
    call realizespin(bondvec(3,:),rbv(3,:))
    call realizespin(bondvec(4,:),rbv(4,:))
    !**************************************

    do i = 1,size1

       !******************BC**********************
       a = mod(i+1+size1-1,size1)+1
       !******************************************

       do j = 1 , size1

          !****************BC*********************
          c = mod(j+1+size2-1,size2)+1
          !***************************************

          do k = 1, size1

            b = mod(k+1+size3-1,size3)+1

             !***********realize**************
             call realizespin(spinA(:,k,j,a,T_ex(T_i)),spinrA(:,1))
             call realizespin(spinA(:,k,j,i,T_ex(T_i)),spinrA(:,2))
             call realizespin(spinA(:,k,c,i,T_ex(T_i)),spinrA(:,3))
             call realizespin(spinA(:,b,j,i,T_i),spinrA(:,4))
             call realizespin(spinB(:,k,j,i,T_ex(T_i)),spinrB)
            !********************************

             tmp1(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,j,a,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,j,a,T_ex(T_i))
             tmp1(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,j,a,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,j,a,T_ex(T_i))
             tmp1(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,a,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,a,T_ex(T_i))

             tmp2(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,j,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,j,i,T_ex(T_i))
             tmp2(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,j,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,j,i,T_ex(T_i))
             tmp2(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,i,T_ex(T_i))

             tmp3(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,c,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,c,i,T_ex(T_i))
             tmp3(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,c,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,c,i,T_ex(T_i))
             tmp3(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,c,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,c,i,T_ex(T_i))

             tmp4(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,b,j,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,b,j,i,T_ex(T_i))
             tmp4(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,b,j,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,b,j,i,T_ex(T_i))
             tmp4(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,b,j,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,b,j,i,T_ex(T_i))

             !****real
             tmpr1(1) = spinrB(2)*spinrA(3,1) - spinrB(3)*spinrA(2,1)
             tmpr1(2) = spinrB(3)*spinrA(1,1) - spinrB(1)*spinrA(3,1)
             tmpr1(3) = spinrB(1)*spinrA(2,1) - spinrB(2)*spinrA(1,1)

             tmpr2(1) = spinrB(2)*spinrA(3,2) - spinrB(3)*spinrA(2,2)
             tmpr2(2) = spinrB(3)*spinrA(1,2) - spinrB(1)*spinrA(3,2)
             tmpr2(3) = spinrB(1)*spinrA(2,2) - spinrB(2)*spinrA(1,2)

             tmpr3(1) = spinrB(2)*spinrA(3,3) - spinrB(3)*spinrA(2,3)
             tmpr3(2) = spinrB(3)*spinrA(1,3) - spinrB(1)*spinrA(3,3)
             tmpr3(3) = spinrB(1)*spinrA(2,3) - spinrB(2)*spinrA(1,3)

             tmpr4(1) = spinrB(2)*spinrA(3,4) - spinrB(3)*spinrA(2,4)
             tmpr4(2) = spinrB(3)*spinrA(1,4) - spinrB(1)*spinrA(3,4)
             tmpr4(3) = spinrB(1)*spinrA(2,4) - spinrB(2)*spinrA(1,4)

             !*****

             chi(T_i,3,:)  = chi(T_i,3,:) + tmpr1!bond3
             chi(T_i,1,:)  = chi(T_i,1,:) + tmpr2!bond1
             chi(T_i,2,:)  = chi(T_i,2,:) + tmpr3!bond2
             chi(T_i,4,:)  = chi(T_i,4,:) + tmpr4!inter

             Schi(T_i,3,:)  = Schi(T_i,3,:) + tmp1!bond3
             Schi(T_i,1,:)  = Schi(T_i,1,:) + tmp2!bond1
             Schi(T_i,2,:)  = Schi(T_i,2,:) + tmp3!bond2
             Schi(T_i,4,:)  = Schi(T_i,4,:) + tmp4!inter

             p(T_i,1) = tmpr1(3) - tmpr2(3)*0.5d0 - tmpr3(3)*0.5d0 - tmpr4(2) + p(T_i,1)
             p(T_i,2) = tmpr2(3)*(sqrt(3.0d0)*0.5d0) - tmpr3(3)*(sqrt(3.0d0)*0.5d0) + tmpr4(1) + p(T_i,2)
             p(T_i,3) = -tmpr1(1)-tmpr2(2)*(sqrt(3.0d0)*0.5d0)+ tmpr2(1)*(0.5d0)+&
                  tmpr3(2)*(sqrt(3.0d0)*0.5d0)+ tmpr3(1)*(0.5d0) + p(T_i,3)

             Sp(T_i,1) = tmp1(3)*rbv(3,2)-tmp1(2)*rbv(3,3)+tmp2(3)*rbv(1,2)-tmp2(2)*rbv(1,3)+&
             tmp3(3)*rbv(2,2)-tmp3(2)*rbv(2,3)+tmp4(3)*rbv(4,2)-tmp4(2)*rbv(4,3)  + Sp(T_i,1)

             Sp(T_i,2) = tmp1(1)*rbv(3,3)-tmp1(3)*rbv(3,1)+tmp2(1)*rbv(1,3)-tmp2(3)*rbv(1,1)+&
             tmp3(1)*rbv(2,3)-tmp3(3)*rbv(2,1)+tmp4(1)*rbv(4,3)-tmp4(3)*rbv(4,1) + Sp(T_i,2)

             Sp(T_i,3) = tmp1(2)*rbv(3,1)-tmp1(1)*rbv(3,2)+tmp2(2)*rbv(1,1)-tmp2(1)*rbv(1,2)+&
             tmp3(2)*rbv(2,1)-tmp3(1)*rbv(2,2)+tmp4(2)*rbv(4,1)-tmp4(1)*rbv(4,2)+ Sp(T_i,3)



          end do
       end do
    end do

  end subroutine calc_chila_pola

  subroutine scP_int_heisen(spinA,spinB,T_i,T_ex,P,SP)!カイラリティー,分極の計算,(崎山さんのものに従いsubA×subBの外積、論文とABが逆なので注意)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: P(N_T,3),SP(N_T,3)
    real(8) tmp1(3),tmp2(3),tmp3(3),tmp4(3),Chi(N_T,4,3),SChi(N_T,4,3)
    real(8)tmpr1(3),tmpr2(3),tmpr3(3),tmpr4(3),spinrA(3,4),spinrB(3)
    real(8) bondvec(4,3),rbv(4,3)
    integer i, j, k,size1 , size2 , size3,pm,a,b,c,d

    size1 = L
    size2 = L
    size3 = L

    !*************bond vec***********
    bondvec(1,1)=-sqrt(3.0d0)*0.5
    bondvec(1,2)=-0.5d0
    bondvec(1,3)=0
    bondvec(2,1)=sqrt(3.0d0)*0.5
    bondvec(2,2)=-0.5d0
    bondvec(2,3)=0.0d0
    bondvec(3,1)=0.0d0
    bondvec(3,2)=1.0d0
    bondvec(3,3)=0.0d0
    bondvec(4,1)=0.0d0
    bondvec(4,2)=0.0d0
    bondvec(4,3)=1.0d0
    !********************************
    !************realize bond vec**********
    call realizespin(bondvec(1,:),rbv(1,:))
    call realizespin(bondvec(2,:),rbv(2,:))
    call realizespin(bondvec(3,:),rbv(3,:))
    call realizespin(bondvec(4,:),rbv(4,:))
    !**************************************

    do i = 1,size1

       !******************BC**********************
       a = mod(i+1+size1-1,size1)+1
       !******************************************

       do j = 1 , size1

          !****************BC*********************
          c = mod(j+1+size2-1,size2)+1
          !***************************************

          do k = 1, size1

            b = mod(k+1+size3-1,size3)+1

             !***********realize**************
             call realizespin(spinA(:,k,j,a,T_ex(T_i)),spinrA(:,1))
             call realizespin(spinA(:,k,j,i,T_ex(T_i)),spinrA(:,2))
             call realizespin(spinA(:,k,c,i,T_ex(T_i)),spinrA(:,3))
             call realizespin(spinA(:,b,j,i,T_i),spinrA(:,4))
             call realizespin(spinB(:,k,j,i,T_ex(T_i)),spinrB)
            !********************************

             tmp1(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,j,a,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,j,a,T_ex(T_i))
             tmp1(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,j,a,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,j,a,T_ex(T_i))
             tmp1(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,a,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,a,T_ex(T_i))

             tmp2(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,j,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,j,i,T_ex(T_i))
             tmp2(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,j,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,j,i,T_ex(T_i))
             tmp2(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,i,T_ex(T_i))

             tmp3(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,c,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,c,i,T_ex(T_i))
             tmp3(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,c,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,c,i,T_ex(T_i))
             tmp3(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,c,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,c,i,T_ex(T_i))

             tmp4(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,b,j,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,b,j,i,T_ex(T_i))
             tmp4(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,b,j,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,b,j,i,T_ex(T_i))
             tmp4(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,b,j,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,b,j,i,T_ex(T_i))

             !****real
             tmpr1(1) = spinrB(2)*spinrA(3,1) - spinrB(3)*spinrA(2,1)
             tmpr1(2) = spinrB(3)*spinrA(1,1) - spinrB(1)*spinrA(3,1)
             tmpr1(3) = spinrB(1)*spinrA(2,1) - spinrB(2)*spinrA(1,1)

             tmpr2(1) = spinrB(2)*spinrA(3,2) - spinrB(3)*spinrA(2,2)
             tmpr2(2) = spinrB(3)*spinrA(1,2) - spinrB(1)*spinrA(3,2)
             tmpr2(3) = spinrB(1)*spinrA(2,2) - spinrB(2)*spinrA(1,2)

             tmpr3(1) = spinrB(2)*spinrA(3,3) - spinrB(3)*spinrA(2,3)
             tmpr3(2) = spinrB(3)*spinrA(1,3) - spinrB(1)*spinrA(3,3)
             tmpr3(3) = spinrB(1)*spinrA(2,3) - spinrB(2)*spinrA(1,3)

             tmpr4(1) = spinrB(2)*spinrA(3,4) - spinrB(3)*spinrA(2,4)
             tmpr4(2) = spinrB(3)*spinrA(1,4) - spinrB(1)*spinrA(3,4)
             tmpr4(3) = spinrB(1)*spinrA(2,4) - spinrB(2)*spinrA(1,4)

             !*****

             chi(T_i,3,:)  = chi(T_i,3,:) + tmp1!bond3
             chi(T_i,1,:)  = chi(T_i,1,:) + tmp2!bond1
             chi(T_i,2,:)  = chi(T_i,2,:) + tmp3!bond2
             chi(T_i,4,:)  = chi(T_i,4,:) + tmp4!inter

             Schi(T_i,3,:)  = Schi(T_i,3,:) + tmp1!bond3
             Schi(T_i,1,:)  = Schi(T_i,1,:) + tmp2!bond1
             Schi(T_i,2,:)  = Schi(T_i,2,:) + tmp3!bond2
             Schi(T_i,4,:)  = Schi(T_i,4,:) + tmp4!inter

             p(T_i,1) = - tmpr4(2) + p(T_i,1)
             p(T_i,2) =  tmpr4(1) + p(T_i,2)
             p(T_i,3) =  + p(T_i,3)


             Sp(T_i,1) = tmp4(3)*rbv(4,2)-tmp4(2)*rbv(4,3)  + Sp(T_i,1)

             Sp(T_i,2) = tmp4(1)*rbv(4,3)-tmp4(3)*rbv(4,1) + Sp(T_i,2)

             Sp(T_i,3) = tmp4(2)*rbv(4,1)-tmp4(1)*rbv(4,2)+ Sp(T_i,3)



          end do
       end do
    end do

  end subroutine scP_int_heisen


  subroutine calc_chila_pola_bond(spinA,spinB,T_i,T_ex,P)!カイラリティー,分極の計算,(崎山さんのものに従いsubA×subBの外積、論文とABが逆なので注意)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: P(N_T,3,3)
    real(8) tmp1(3),tmp2(3),tmp3(3),tmp4(3)
    real(8)tmpr1(3),tmpr2(3),tmpr3(3),tmpr4(3),spinrA(3,4),spinrB(3)
    real(8) bondvec(4,3),rbv(4,3)
    integer i, j, k,size1 , size2 , size3,pm,a,b,c,d

    size1 = L
    size2 = L
    size3 = L

    !*************bond vec***********
    bondvec(1,1)=-sqrt(3.0d0)*0.5
    bondvec(1,2)=-0.5d0
    bondvec(1,3)=0
    bondvec(2,1)=sqrt(3.0d0)*0.5
    bondvec(2,2)=-0.5d0
    bondvec(2,3)=0.0d0
    bondvec(3,1)=0.0d0
    bondvec(3,2)=1.0d0
    bondvec(3,3)=0.0d0
    bondvec(4,1)=0.0d0
    bondvec(4,2)=0.0d0
    bondvec(4,3)=1.0d0
    !********************************
    !************realize bond vec**********
    call realizespin(bondvec(1,:),rbv(1,:))
    call realizespin(bondvec(2,:),rbv(2,:))
    call realizespin(bondvec(3,:),rbv(3,:))
    call realizespin(bondvec(4,:),rbv(4,:))
    !**************************************

    do i = 1,size1

       !******************BC**********************
       a = mod(i+1+size1-1,size1)+1
       !******************************************

       do j = 1 , size1

          !****************BC*********************
          c = mod(j+1+size2-1,size2)+1
          !***************************************

          do k = 1, size1

            b = mod(k+1+size3-1,size3)+1

             !***********realize**************
             call realizespin(spinA(:,k,j,a,T_ex(T_i)),spinrA(:,1))
             call realizespin(spinA(:,k,j,i,T_ex(T_i)),spinrA(:,2))
             call realizespin(spinA(:,k,c,i,T_ex(T_i)),spinrA(:,3))
             call realizespin(spinA(:,b,j,i,T_i),spinrA(:,4))
             call realizespin(spinB(:,k,j,i,T_ex(T_i)),spinrB)
            !********************************

             tmp1(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,j,a,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,j,a,T_ex(T_i))
             tmp1(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,j,a,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,j,a,T_ex(T_i))
             tmp1(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,a,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,a,T_ex(T_i))

             tmp2(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,j,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,j,i,T_ex(T_i))
             tmp2(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,j,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,j,i,T_ex(T_i))
             tmp2(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,i,T_ex(T_i))

             tmp3(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,k,c,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,k,c,i,T_ex(T_i))
             tmp3(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,k,c,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,k,c,i,T_ex(T_i))
             tmp3(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,c,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,c,i,T_ex(T_i))

             tmp4(1) = spinB(2,k,j,i,T_ex(T_i))*spinA(3,b,j,i,T_ex(T_i)) - spinB(3,k,j,i,T_ex(T_i))*spinA(2,b,j,i,T_ex(T_i))
             tmp4(2) = spinB(3,k,j,i,T_ex(T_i))*spinA(1,b,j,i,T_ex(T_i)) - spinB(1,k,j,i,T_ex(T_i))*spinA(3,b,j,i,T_ex(T_i))
             tmp4(3) = spinB(1,k,j,i,T_ex(T_i))*spinA(2,b,j,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,b,j,i,T_ex(T_i))

             !****real
             tmpr1(1) = spinrB(2)*spinrA(3,1) - spinrB(3)*spinrA(2,1)
             tmpr1(2) = spinrB(3)*spinrA(1,1) - spinrB(1)*spinrA(3,1)
             tmpr1(3) = spinrB(1)*spinrA(2,1) - spinrB(2)*spinrA(1,1)

             tmpr2(1) = spinrB(2)*spinrA(3,2) - spinrB(3)*spinrA(2,2)
             tmpr2(2) = spinrB(3)*spinrA(1,2) - spinrB(1)*spinrA(3,2)
             tmpr2(3) = spinrB(1)*spinrA(2,2) - spinrB(2)*spinrA(1,2)

             tmpr3(1) = spinrB(2)*spinrA(3,3) - spinrB(3)*spinrA(2,3)
             tmpr3(2) = spinrB(3)*spinrA(1,3) - spinrB(1)*spinrA(3,3)
             tmpr3(3) = spinrB(1)*spinrA(2,3) - spinrB(2)*spinrA(1,3)

             tmpr4(1) = spinrB(2)*spinrA(3,4) - spinrB(3)*spinrA(2,4)
             tmpr4(2) = spinrB(3)*spinrA(1,4) - spinrB(1)*spinrA(3,4)
             tmpr4(3) = spinrB(1)*spinrA(2,4) - spinrB(2)*spinrA(1,4)

             !*****

             p(T_i,1,1) = tmpr2(3)*0.5d0 + p(T_i,1,1)
             p(T_i,2,1) = tmpr2(3)*(sqrt(3.0d0)*0.5d0)+ p(T_i,2,1)
             p(T_i,3,1) = -tmpr2(2)*(sqrt(3.0d0)*0.5d0)+ tmpr2(1)*(0.5d0) + p(T_i,3,1)

             p(T_i,1,2) =  - tmpr3(3)*0.5d0 + p(T_i,1,2)
             p(T_i,2,2) =  - tmpr3(3)*(sqrt(3.0d0)*0.5d0) + p(T_i,2,2)
             p(T_i,3,2) = tmpr3(2)*(sqrt(3.0d0)*0.5d0)+ tmpr3(1)*(0.5d0) + p(T_i,3,2)

             p(T_i,1,3) = tmpr1(3)+ p(T_i,1,3)
             p(T_i,2,3) =  p(T_i,2,3)
             p(T_i,3,3) = -tmpr1(1) + p(T_i,3,3)


          end do
       end do
    end do

  end subroutine calc_chila_pola_bond



  subroutine troidal(rg,spinA,spinB,T_i,T_ex,m,tm,P,Stm,SP)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:),m(:,:),rg(:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: tm(N_T,3),P(N_T,3),Stm(N_T,3),SP(N_T,3)
    real(8) base(4,3),posi(3),tmp(3),tmx,tmy,tmz
    real(8)rposi(3),tmpr(3),rtmx,rtmy,rtmz,rH(3)
    integer i, j, k, size1 , size2 , size3 , a, b

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    !*******************subA toidal**************************************
    do i = 1,size1

       a = i-1
       b = 0

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             posi = posi - rg
             tmp = spinA(:,k,j,i,T_ex(T_i)) !- m(T_i,:)

             tmx = posi(2)*tmp(3) - posi(3)*tmp(2)
             tmy = posi(3)*tmp(1) - posi(1)*tmp(3)
             tmz = posi(1)*tmp(2) - posi(2)*tmp(1)

             Stm(T_i,1) = Stm(T_i,1) + tmx
             Stm(T_i,2) = Stm(T_i,2) + tmy
             Stm(T_i,3) = Stm(T_i,3) + tmz

             SP(T_i,1) = SP(T_i,1) -tmz*H(2) + tmy*H(3)
             SP(T_i,2) = SP(T_i,2) + tmz*H(1) -tmx*H(3)
             SP(T_i,3) = SP(T_i,3) + tmx*H(2) - tmy*H(1)

             call realizespin(posi,rposi)
             call realizespin(tmp,tmpr)
             call realizespin(H,rH)

             rtmx = rposi(2)*tmpr(3) - rposi(3)*tmpr(2)
             rtmy = rposi(3)*tmpr(1) - rposi(1)*tmpr(3)
             rtmz = rposi(1)*tmpr(2) - rposi(2)*tmpr(1)

             tm(T_i,1) = tm(T_i,1) + rtmx
             tm(T_i,2) = tm(T_i,2) + rtmy
             tm(T_i,3) = tm(T_i,3) + rtmz

             P(T_i,1) = P(T_i,1) -rtmz*rH(2) + rtmy*rH(3)
             P(T_i,2) = P(T_i,2) + rtmz*rH(1) -rtmx*rH(3)
             P(T_i,3) = P(T_i,3) + rtmx*rH(2) - rtmy*rH(1)



          end do
       end do
    end do
    !************************************************************************


    !*********************sunB troidal***************************************
    do i = 1,size1

       a = i -1 !base1
       b = 1 !base3

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             posi = posi -rg
             tmp = spinB(:,k,j,i,T_ex(T_i)) !- m(T_i,:)

             tmx = posi(2)*tmp(3) - posi(3)*tmp(2)
             tmy = posi(3)*tmp(1) - posi(1)*tmp(3)
             tmz = posi(1)*tmp(2) - posi(2)*tmp(1)

             Stm(T_i,1) = Stm(T_i,1) + tmx
             Stm(T_i,2) = Stm(T_i,2) + tmy
             Stm(T_i,3) = Stm(T_i,3) + tmz

             SP(T_i,1) = SP(T_i,1) -tmz*H(2) + tmy*H(3)
             SP(T_i,2) = SP(T_i,2) + tmz*H(1) -tmx*H(3)
             SP(T_i,3) = SP(T_i,3) + tmx*H(2) - tmy*H(1)

             call realizespin(posi,rposi)
             call realizespin(tmp,tmpr)
             call realizespin(H,rH)

             rtmx = rposi(2)*tmpr(3) - rposi(3)*tmpr(2)
             rtmy = rposi(3)*tmpr(1) - rposi(1)*tmpr(3)
             rtmz = rposi(1)*tmpr(2) - rposi(2)*tmpr(1)

             tm(T_i,1) = tm(T_i,1) + rtmx
             tm(T_i,2) = tm(T_i,2) + rtmy
             tm(T_i,3) = tm(T_i,3) + rtmz

             P(T_i,1) = P(T_i,1) -rtmz*rH(2) + rtmy*rH(3)
             P(T_i,2) = P(T_i,2) + rtmz*rH(1) -rtmx*rH(3)
             P(T_i,3) = P(T_i,3) + rtmx*rH(2) - rtmy*rH(1)


          end do
       end do
    end do
    !***************************************************************************

  end subroutine troidal


  subroutine troidal2(spinA,spinB,T_i,T_ex,m,tm,P)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:),m(:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: tm(N_T,3),P(N_T,3)
    real(8) base(4,3),posi(3),tmp(2),tmx,tmy,tmz
    integer i, j, k, size1 , size2 , size3 , a, b

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    !*******************subA toidal**************************************
    do i = 1,size1

       a = i-1
       b = 0

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             tmp = spinA(:,k,j,i,T_ex(T_i))

             tmx = -posi(3)*tmp(2)
             tmy = posi(3)*tmp(1)
             tmz = posi(1)*tmp(2) - posi(2)*tmp(1)

             tm(T_i,1) = tm(T_i,1) + tmx
             tm(T_i,2) = tm(T_i,2) + tmy
             tm(T_i,3) = tm(T_i,3) + tmz

             P(T_i,1) = P(T_i,1) -tmz*H(2)
             P(T_i,2) = P(T_i,2) + tmz*H(1)
             P(T_i,3) = P(T_i,3) + tmx*H(2) - tmy*H(1)

          end do
       end do
    end do
    !************************************************************************


    !*********************sunB troidal***************************************
    do i = 1,size1

       a = i -1 !base1
       b = 1 !base3

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             tmp = spinB(:,k,j,i,T_ex(T_i))

             tmx = -posi(3)*tmp(2)
             tmy = posi(3)*tmp(1)
             tmz = posi(1)*tmp(2) - posi(2)*tmp(1)

             tm(T_i,1) = tm(T_i,1) + tmx
             tm(T_i,2) = tm(T_i,2) + tmy
             tm(T_i,3) = tm(T_i,3) + tmz

             P(T_i,1) = P(T_i,1) -tmz*H(2)
             P(T_i,2) = P(T_i,2) + tmz*H(1)
             P(T_i,3) = P(T_i,3) + tmx*H(2) - tmy*H(1)

          end do
       end do
    end do
    !***************************************************************************

  end subroutine troidal2



  subroutine spin_overlap(spinA,spin_repA,spinB,spin_repB,T_i,T_ex,T_ex_rep,q_uv,q_uv_kmin) !spin overlap
    real(8),intent(in) :: spinA(:,:,:,:,:),spin_repA(:,:,:,:,:), spinB(:,:,:,:,:),spin_repB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:),T_ex_rep(:)
    real(8),intent(out) :: q_uv(:,:),q_uv_kmin(:,:,:) !(N_T,4) tensor
    real(8) posi(3),base(4,3),r
    integer i,j,k,a,b,size1,size2,size3

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    size1 = L
    size2 = L
    size3 = L

    !****************************************spin overlap subA*********************************************

    do i = 1,size1

       a = i-1
       b = 0

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             r = twopi/dble(L)*posi(1)

             q_uv(T_i,1) = spinA(1,k,j,i,T_ex(T_i))*spin_repA(1,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,1)!xx
             q_uv(T_i,2) = spinA(1,k,j,i,T_ex(T_i))*spin_repA(2,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,2)!xy
             q_uv(T_i,3) = spinA(2,k,j,i,T_ex(T_i))*spin_repA(1,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,3)!yx
             q_uv(T_i,4) = spinA(2,k,j,i,T_ex(T_i))*spin_repA(2,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,4)!yy

             q_uv_kmin(T_i,1,1) = spinA(1,k,j,i,T_ex(T_i))*spin_repA(1,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,1,1)!xx
             q_uv_kmin(T_i,1,2) = spinA(1,k,j,i,T_ex(T_i))*spin_repA(1,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,1,2)!xx
             q_uv_kmin(T_i,2,1) = spinA(1,k,j,i,T_ex(T_i))*spin_repA(2,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,2,1)!xy
             q_uv_kmin(T_i,2,2) = spinA(1,k,j,i,T_ex(T_i))*spin_repA(2,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,2,2)!xy
             q_uv_kmin(T_i,3,1) = spinA(2,k,j,i,T_ex(T_i))*spin_repA(1,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,3,1)!yx
             q_uv_kmin(T_i,3,2) = spinA(2,k,j,i,T_ex(T_i))*spin_repA(1,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,3,2)!yx
             q_uv_kmin(T_i,4,1) = spinA(2,k,j,i,T_ex(T_i))*spin_repA(2,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,4,1)!yy
             q_uv_kmin(T_i,4,2) = spinA(2,k,j,i,T_ex(T_i))*spin_repA(2,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,4,2)!yy

          end do

       end do

    end do
    !********************************************************************************************************

    !****************************************spin overlap subB*********************************************

    do i = 1,size1

       a = i-1
       b = 1

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             r = twopi/dble(L)*posi(1)

             q_uv(T_i,1) = spinB(1,k,j,i,T_ex(T_i))*spin_repB(1,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,1)!xx
             q_uv(T_i,2) = spinB(1,k,j,i,T_ex(T_i))*spin_repB(2,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,2)!xy
             q_uv(T_i,3) = spinB(2,k,j,i,T_ex(T_i))*spin_repB(1,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,3)!yx
             q_uv(T_i,4) = spinB(2,k,j,i,T_ex(T_i))*spin_repB(2,k,j,i,T_ex_rep(T_i))+ q_uv(T_i,4)!yy

             q_uv_kmin(T_i,1,1) = spinB(1,k,j,i,T_ex(T_i))*spin_repB(1,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,1,1)!xx
             q_uv_kmin(T_i,1,2) = spinB(1,k,j,i,T_ex(T_i))*spin_repB(1,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,1,2)!xx
             q_uv_kmin(T_i,2,1) = spinB(1,k,j,i,T_ex(T_i))*spin_repB(2,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,2,1)!xy
             q_uv_kmin(T_i,2,2) = spinB(1,k,j,i,T_ex(T_i))*spin_repB(2,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,2,2)!xy
             q_uv_kmin(T_i,3,1) = spinB(2,k,j,i,T_ex(T_i))*spin_repB(1,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,3,1)!yx
             q_uv_kmin(T_i,3,2) = spinB(2,k,j,i,T_ex(T_i))*spin_repB(1,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,3,2)!yx
             q_uv_kmin(T_i,4,1) = spinB(2,k,j,i,T_ex(T_i))*spin_repB(2,k,j,i,T_ex_rep(T_i))*cos(r)+q_uv_kmin(T_i,4,1)!yy
             q_uv_kmin(T_i,4,2) = spinB(2,k,j,i,T_ex(T_i))*spin_repB(2,k,j,i,T_ex_rep(T_i))*sin(r)+q_uv_kmin(T_i,4,2)!yy

          end do

       end do

    end do
    !********************************************************************************************************
    q_uv(T_i,:) = q_uv(T_i,:)/dble(N_spin)
    q_uv_kmin(T_i,:,:) = q_uv_kmin(T_i,:,:)/dble(N_spin)

  end subroutine spin_overlap






  subroutine chiral_overlap(spinA,spin_repA,spinB,spin_repB,T_i,T_ex,T_ex_rep,q_uv,q_uv_kmin)
    real(8),intent(in) :: spinA(:,:,:,:,:),spin_repA(:,:,:,:,:), spinB(:,:,:,:,:),spin_repB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:),T_ex_rep(:)
    integer size1,size2,size3,i,j,k,a,c
    real(8),intent(out) :: q_uv(:,:),q_uv_kmin(:,:,:) !(N_T,4) tensor
    real(8) posi(3),base(4,3),r,r1,r2,r3
    real(8)chi,chi1,chi2,chi3,chi_rep,chi1_rep,chi2_rep,chi3_rep

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    !****************************************chiral overlap subB*********************************************

    do i = 1,size1

        !******************BC**********************
       a = mod(i+1+size1-1,size1)+1
       !******************************************


       do j = 1 , size1

          !****************BC*********************
          c = mod(j+1+size2-1,size2)+1
          !***************************************


          do k = 1,size3

             posi = (i-1)*base(1,:) +(j-1)*base(2,:) + base(3,:) + (k-1)*base(4,:)
             r = twopi/dble(L)*posi(1)!site posi_x
             r1 = twopi/dble(L)*(posi(1)-sqrt(3.0d0)/2)!bond1 posi_x
             r2 = twopi/dble(L)*(posi(1)+sqrt(3.0d0)/2)!bond2 posi_x
             r3 = twopi/dble(L)*posi(1)!bond3 posi_x

             chi1 = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,a,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,a,T_ex(T_i))
             chi2 = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,j,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,j,i,T_ex(T_i))
             chi3 = spinB(1,k,j,i,T_ex(T_i))*spinA(2,k,c,i,T_ex(T_i)) - spinB(2,k,j,i,T_ex(T_i))*spinA(1,k,c,i,T_ex(T_i))
             chi = chi1 + chi2 + chi3

             chi1_rep = spin_repB(1,k,j,i,T_ex_rep(T_i))*spin_repA(2,k,j,a,T_ex_rep(T_i))&
                  - spin_repB(2,k,j,i,T_ex_rep(T_i))*spin_repA(1,k,j,a,T_ex_rep(T_i))
             chi2_rep = spin_repB(1,k,j,i,T_ex_rep(T_i))*spin_repA(2,k,j,i,T_ex_rep(T_i))&
                  - spin_repB(2,k,j,i,T_ex_rep(T_i))*spin_repA(1,k,j,i,T_ex_rep(T_i))
             chi3_rep = spin_repB(1,k,j,i,T_ex_rep(T_i))*spin_repA(2,k,c,i,T_ex_rep(T_i))&
                  - spin_repB(2,k,j,i,T_ex_rep(T_i))*spin_repA(1,k,c,i,T_ex_rep(T_i))
             chi_rep = chi1_rep+chi2_rep+chi3_rep

             q_uv(T_i,1)= q_uv(T_i,1)+ (chi1*chi1_rep) + (chi2*chi2_rep) +(chi3*chi3_rep)  !bond overlap
             q_uv(T_i,2)= q_uv(T_i,2)+ (chi*chi_rep) !site overlap
             q_uv_kmin(T_i,1,1) =q_uv_kmin(T_i,1,1)+ chi1*chi1_rep*cos(r1) + chi2*chi2_rep*cos(r2) +chi3*chi3_rep*cos(r3)
             q_uv_kmin(T_i,1,2) =q_uv_kmin(T_i,1,2)+ chi1*chi1_rep*sin(r1) + chi2*chi2_rep*sin(r2) +chi3*chi3_rep*sin(r3)
             q_uv_kmin(T_i,2,1) =q_uv_kmin(T_i,2,1)+ chi*chi_rep*cos(r)
             q_uv_kmin(T_i,2,2) =q_uv_kmin(T_i,2,2)+ chi*chi_rep*sin(r)



          end do
       end do
    end do

    q_uv(T_i,1) = q_uv(T_i,1) / dble(L*L*L*3)
    q_uv(T_i,2) = q_uv(T_i,2) / dble(L*L*L*2)
    q_uv_kmin(T_i,1,:) = q_uv_kmin(T_i,1,:)/dble(L*L*L*3)
    q_uv_kmin(T_i,2,:) = q_uv_kmin(T_i,2,:)/dble(L*L*L*2)


  end subroutine chiral_overlap


  subroutine monpole(rg,spinA,spinB,T_i,T_ex,m,mono)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:),m(:,:),rg(:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: mono(:)
    real(8) base(4,3),posi(3),tmp(3),tmx,tmy,tmz
    integer i, j, k, size1 , size2 , size3 , a, b

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base/sqrt(3.0d0)
    tmp(3) =0.0d0
    mono(T_i) =0.0d0

    !*******************subA monopole**************************************
    do i = 1,size1

       a = i-1
       b = 0

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + k*base(4,:)
             posi = posi -rg
             tmp(1) = spinA(1,k,j,i,T_ex(T_i)) - m(T_i,1)
             tmp(2) = spinA(2,k,j,i,T_ex(T_i)) - m(T_i,2)
             tmp(3) = spinA(3,k,j,i,T_ex(T_i)) - m(T_i,3)
             mono(T_i) = dot_product(posi,tmp) + mono(T_i)

          end do
       end do
    end do
    !************************************************************************


    !*********************sunB monopole***************************************
    do i = 1,size1

       a = i -1 !base1
       b = 1 !base3

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + k*base(4,:)
             posi = posi -rg
             tmp(1) = spinB(1,k,j,i,T_ex(T_i)) - m(T_i,1)
             tmp(2) = spinB(2,k,j,i,T_ex(T_i)) - m(T_i,2)
             tmp(3) = spinB(3,k,j,i,T_ex(T_i)) - m(T_i,3)

             mono(T_i) = dot_product(posi,tmp) + mono(T_i)

          end do
       end do
    end do
    !***************************************************************************

  end subroutine monpole



  subroutine quadpole(rg,spinA,spinB,T_i,T_ex,m,quad)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:),m(:,:),rg(:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: quad(:,:)
    real(8) base(4,3),posi(3),tmp(3),tmx,tmy,tmz
    integer i, j, k, size1 , size2 , size3 , a, b

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base/sqrt(3.0d0)
    tmp(3) =0.0d0
    quad(T_i,:) =0.0d0

    !*******************subA quadpole**************************************
    do i = 1,size1

       a = i-1
       b = 0

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + k*base(4,:)
             posi = posi -rg
             tmp(1) = spinA(1,k,j,i,T_ex(T_i)) - m(T_i,1)
             tmp(2) = spinA(2,k,j,i,T_ex(T_i)) - m(T_i,2)
             tmp(3) = spinA(3,k,j,i,T_ex(T_i)) - m(T_i,3)

             quad(T_i,1) = quad(T_i,1)+ 2*posi(1)*tmp(1) -2.0d0/3.0d0*dot_product(posi,tmp)
             quad(T_i,2) = quad(T_i,2)+ posi(1)*tmp(2) + posi(2)*tmp(1)
             quad(T_i,3) = quad(T_i,3)+ posi(1)*tmp(3) + posi(3)*tmp(1)
             quad(T_i,4) = quad(T_i,4)+ posi(2)*tmp(1) + posi(1)*tmp(2)
             quad(T_i,5) = quad(T_i,5)+ 2*posi(2)*tmp(2) -2.0d0/3.0d0*dot_product(posi,tmp)
             quad(T_i,6) = quad(T_i,6)+ posi(2)*tmp(3) + posi(3)*tmp(2)
             quad(T_i,7) = quad(T_i,7)+ posi(3)*tmp(1) + posi(1)*tmp(3)
             quad(T_i,8) = quad(T_i,8)+ posi(3)*tmp(2) + posi(2)*tmp(3)
             quad(T_i,9) = quad(T_i,9)+ 2*posi(3)*tmp(3) -2.0d0/3.0d0*dot_product(posi,tmp)


          end do
       end do
    end do
    !************************************************************************


    !*********************sunB quadpaole***************************************
    do i = 1,size1

       a = i -1 !base1
       b = 1 !base3

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + k*base(4,:)
             posi = posi -rg
             tmp(1) = spinB(1,k,j,i,T_ex(T_i)) - m(T_i,1)
             tmp(2) = spinB(2,k,j,i,T_ex(T_i)) - m(T_i,2)
             tmp(3) = spinB(3,k,j,i,T_ex(T_i)) - m(T_i,3)

             quad(T_i,1) = quad(T_i,1)+ 2*posi(1)*tmp(1) -2.0d0/3.0d0*dot_product(posi,tmp)
             quad(T_i,2) = quad(T_i,2)+ posi(1)*tmp(2) + posi(2)*tmp(1)
             quad(T_i,3) = quad(T_i,3)+ posi(1)*tmp(3) + posi(3)*tmp(1)
             quad(T_i,4) = quad(T_i,4)+ posi(2)*tmp(1) + posi(1)*tmp(2)
             quad(T_i,5) = quad(T_i,5)+ 2*posi(2)*tmp(2) -2.0d0/3.0d0*dot_product(posi,tmp)
             quad(T_i,6) = quad(T_i,6)+ posi(2)*tmp(3) + posi(3)*tmp(2)
             quad(T_i,7) = quad(T_i,7)+ posi(3)*tmp(1) + posi(1)*tmp(3)
             quad(T_i,8) = quad(T_i,8)+ posi(3)*tmp(2) + posi(2)*tmp(3)
             quad(T_i,9) = quad(T_i,9)+ 2*posi(3)*tmp(3) -2.0d0/3.0d0*dot_product(posi,tmp)


          end do
       end do
    end do
    !***************************************************************************

  end subroutine quadpole

   subroutine r_gcalc(r_g)
    real(8),intent(inout) :: r_g(:)
    real(8) base(4,3),posi(3)
    integer i, j, k, size1 , size2 , size3 , a, b

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)
    posi = 0.0d0

    !*******************subA toidal**************************************
    do i = 1,size1

       a = i-1
       b = 0

       do j = 1,size2
          do k = 1,size3

             posi = posi +  a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)

          end do
       end do
    end do
    !************************************************************************


    !*********************sunB troidal***************************************
    do i = 1,size1

       a = i -1 !base1
       b = 1 !base3

       do j = 1,size2
          do k = 1,size3

             posi = posi + a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)

          end do
       end do
    end do
    !***************************************************************************

    r_g = posi/N_spin

  end subroutine r_gcalc


  subroutine staggeredvec(spinA,spinB,T_i,T_ex,vec)!スタッガード磁化を返す
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: vec(:,:)
    integer i, j, k, size1 , size2 , size3

    size1 = L
    size2 = L
    size3 = L


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

             vec(T_i,:) = vec(T_i,:) -  spinA(:,k,j,i,T_ex(T_i)) + spinB(:,k,j,i,T_ex(T_i))

          end do
       end do
    end do

  end subroutine staggeredvec

  subroutine ac(rg,spinA,spinB,spinA0,spinB0,m0,T_i,T_ex,m,spin_ac,toro_ac)
    real(8),intent(in) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:),m(:,:),rg(:)
    real(8),intent(in) :: spinA0(:,:,:,:),spinB0(:,:,:,:),m0(:)
    integer,intent(in) :: T_i,T_ex(:)
    real(8),intent(out) :: spin_ac,toro_ac(:)
    real(8) base(4,3),posi(3),tmp(2),tmx,tmy,tmz
    real(8) tmx2,tmy2,tmz2,tmp2(3)
    integer i, j, k, size1 , size2 , size3 , a, b

    size1 = L
    size2 = L
    size3 = L

    !********並進ベクトル******
    base(1,2) = 3.0d0/2
    base(1,1) = sqrt(3.0d0)/2
    base(1,3) = 0.0d0
    base(2,1) = sqrt(3.0d0)
    base(2,2) = 0.0d0
    base(2,3) = 0.0d0
    base(3,1) = sqrt(3.0d0)/2
    base(3,2) = 1.0d0/2
    base(3,3) = 0.0d0
    base(4,1) = sqrt(3.0d0)/2
    base(4,2) = 0.5d0
    base(4,3) = 1.0d0
    !*************************

    base = base /sqrt(3.0d0)

    spin_ac = 0.0d0
    toro_ac = 0.0d0

    !*******************subA toidal**************************************
    do i = 1,size1

       a = i-1
       b = 0

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             posi = posi - rg
             posi = posi/dble(L)
             tmp = spinA(:,k,j,i,T_ex(T_i)) - m(T_i,:)

             tmx = -posi(3)*tmp(2)
             tmy = posi(3)*tmp(1)
             tmz = posi(1)*tmp(2) - posi(2)*tmp(1)

             tmp2 = spinA0(:,k,j,i) - m0

             tmx2 = -posi(3)*tmp2(2)
             tmy2 = posi(3)*tmp2(1)
             tmz2 = posi(1)*tmp2(2) - posi(2)*tmp2(1)

             toro_ac(1) = toro_ac(1) + tmx*tmx2 + tmy*tmy2 + tmz*tmz2
             toro_ac(2) = toro_ac(2) + tmx*tmx2
             toro_ac(3) = toro_ac(3) + tmy*tmy2
             toro_ac(4) = toro_ac(4) + tmz*tmz2

             spin_ac = spin_ac + dot_product(spinA(:,k,j,i,T_ex(T_i)),spinA0(:,k,j,i))


          end do
       end do
    end do
    !************************************************************************


    !*********************sunB troidal***************************************
    do i = 1,size1

       a = i -1 !base1
       b = 1 !base3

       do j = 1,size2
          do k = 1,size3

             posi = a*base(1,:) +(j-1)*base(2,:) + b*base(3,:) + (k-1)*base(4,:)
             posi = posi -rg
             posi = posi/dble(L)
             tmp = spinB(:,k,j,i,T_ex(T_i)) - m(T_i,:)

             tmx = -posi(3)*tmp(2)
             tmy = posi(3)*tmp(1)
             tmz = posi(1)*tmp(2) - posi(2)*tmp(1)



             tmp2 = spinB0(:,k,j,i) - m0

             tmx2 = -posi(3)*tmp2(2)
             tmy2 = posi(3)*tmp2(1)
             tmz2 = posi(1)*tmp2(2) - posi(2)*tmp2(1)

             toro_ac(1) = toro_ac(1) + tmx*tmx2 + tmy*tmy2 + tmz*tmz2
             toro_ac(2) = toro_ac(2) + tmx*tmx2
             toro_ac(3) = toro_ac(3) + tmy*tmy2
             toro_ac(4) = toro_ac(4) + tmz*tmz2

             spin_ac = spin_ac + dot_product(spinB(:,k,j,i,T_ex(T_i)),spinB0(:,k,j,i))

          end do
       end do
    end do
    !***************************************************************************

    spin_ac = spin_ac/dble(N_spin)
    toro_ac = toro_ac/dble(N_spin)



  end subroutine ac


end module syscalc

module metropolis_method
  use mtmod
  use syscalc
  !$ use omp_lib
  implicit none
contains

  subroutine metropolis(rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
    integer,intent(in) :: T_i
    real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    real(8),intent(in) ::J_ij(:,:,:,:),alpha,randnum(:,:),T(:),rg(:)
    real(8) Eold,Enew,dE,temp(3),rpi,r1,r2,r3,tmp
    real(8) tmp1,tmp2
    integer i, j, k, size1 , size2 , size3 ,count

    size1 = L
    size2 = L
    size3 = L
    count = 1

    !tmp = 1.0d0/T_min*(alpha**(T(T_i)-1))
    tmp = T(N_T)

    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

             !***********spinAの更新**************
             r1 = randnum(count,T_i)*2.0d0-1.0d0
             r2 = randnum(count+1,T_i)*twopi
             r3 = randnum(count+2,T_i)
             tmp1 = sqrt(1.0d0 - r1*r1)
             count = count + 3        !乱数で方向を改めてみる
             call local_EA(rg,spinA,spinB,J_ij,T_i,i,j,k,Eold)
             temp(1) = spinA(1,k,j,i,T_i)!前の情報を一時保持
             temp(2) = spinA(2,k,j,i,T_i)
             temp(3) = spinA(3,k,j,i,T_i)
             spinA(1,k,j,i,T_i) = tmp1*cos(r2)
             spinA(2,k,j,i,T_i) = tmp1*sin(r2)
             spinA(3,k,j,i,T_i) = r1
             call local_EA(rg,spinA,spinB,J_ij,T_i,i,j,k,Enew)
             dE = Enew - Eold
             !反転の判定
             if(dE <= 0.0d0) then

             else if (r3 < exp(-dE/tmp)) then

             else !rの確率で以前の状態にする
                spinA(1,k,j,i,T_i) = temp(1)
                spinA(2,k,j,i,T_i) = temp(2)
                spinA(3,k,j,i,T_i) = temp(3)
             end if
             !**********spinAの更新終わり***********

          end do
       end do
    end do


    do i = 1,size1
       do j = 1,size2
          do k = 1,size3

             !***********spinBの更新**************
             r1 = randnum(count,T_i)*2.0d0-1.0d0
             r2 = randnum(count+1,T_i)*twopi
             r3 = randnum(count+2,T_i)
             tmp1 = sqrt(1.0d0 - r1*r1)
             count = count + 3 !乱数で方向を改めてみる
             call local_EB(rg,spinB,spinA,J_ij,T_i,i,j,k,Eold)
             temp(1) = spinB(1,k,j,i,T_i)!前の情報を一時保持
             temp(2) = spinB(2,k,j,i,T_i)
             temp(3) = spinB(3,k,j,i,T_i)
             spinB(1,k,j,i,T_i) = tmp1*cos(r2)
             spinB(2,k,j,i,T_i) = tmp1*sin(r2)
             spinB(3,k,j,i,T_i) = r1
             call local_EB(rg,spinB,spinA,J_ij,T_i,i,j,k,Enew)
             dE = Enew - Eold
             !反転の判定
             if(dE <= 0.0d0) then

             else if (r3 < exp(-dE/tmp)) then

             else !rの確率で以前の状態にする
                spinB(1,k,j,i,T_i) = temp(1)
                spinB(2,k,j,i,T_i) = temp(2)
                spinB(3,k,j,i,T_i) = temp(3)
             end if
             !**********spinBの更新終わり***********

          end do
       end do
    end do

  end subroutine metropolis


  subroutine over_relaxationA(i,j,k,rg,spinA,spinB,T_i,J_ij)
    real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    real(8),intent(in) :: J_ij(:,:,:,:),rg(:)
    integer ,intent(in) :: T_i,i,j,k
    real(8) local(3),local2,dot,kitvec(3,3),kitaev(3)
    real(8)base(4,3),posi(3),tmp(2),tmx,tmy,tmz,P(3)
    integer  size1 , size2 , size3,R_i
    integer a,c,d

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

       !*********************OR subA***********************************************


          !**********周期境界条件**********************(詳しくは別プログラム)
          a = mod(i-1+L-1,L)+1
          !******************************************



             !**********周期境界条件**********************(詳しくは別プログラム)
             d = mod(j-1+L-1,L)+1
             !******************************************



                !********周期境界条件*********************

                c = mod(k-1+L-1,L)+1
                !***************************************


                !***********position**************************************
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

                local(1) = local(1)+(spinB(3,k,d,i,T_i)+spinB(2,k,j,a,T_i))*gam
                local(2) = local(2)+(spinB(3,k,j,i,T_i)+spinB(1,k,j,a,T_i))*gam
                local(3) = local(3)+(spinB(2,k,j,i,T_i)+spinB(1,k,d,i,T_i))*gam


                local = local + kitaev*K_0-h0*H

                local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
                local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

                local2 = dot_product(local,local)
                dot = dot_product(spinA(:,k,j,i,T_i),local(:))/local2



                spinA(:,k,j,i,T_i) = 2.0d0*dot*local(:) - spinA(:,k,j,i,T_i)!スピンの向きをエネルギーを変えずに変更

              end subroutine over_relaxationA


       !***************************************************************************

       !*********************OR subB***********************************************

       subroutine over_relaxationB(i,j,k,rg,spinA,spinB,T_i,J_ij)
         real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
         real(8),intent(in) :: J_ij(:,:,:,:),rg(:)
         integer ,intent(in) :: T_i,i,j,k
         real(8) local(3),local2,dot,kitvec(3,3),kitaev(3)
         real(8)base(4,3),posi(3),tmp(2),tmx,tmy,tmz,P(3)
         integer  size1 , size2 , size3,R_i
         integer a,c,d

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

          !**********周期境界条件**********************(詳しくは別プログラム)
          a = mod(i+1+L-1,L)+1
          !******************************************



             !**********周期境界条件**********************(詳しくは別プログラム)
             d = mod(j+1+size2-1,size2)+1
             !******************************************



                !********周期境界条件*********************

                c = mod(k+1+L-1,L)+1
                !***************************************

                !***********position**************************************
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

                local(1) = local(1)+(spinA(3,k,d,i,T_i)+spinA(2,k,j,a,T_i))*gam
                local(2) = local(2)+(spinA(3,k,j,i,T_i)+spinA(1,k,j,a,T_i))*gam
                local(3) = local(3)+(spinA(2,k,j,i,T_i)+spinA(1,k,d,i,T_i))*gam

                local = local + kitaev*K_0 -h0*H

                local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
                local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))



                local2 = dot_product(local,local)
                dot = dot_product(spinB(:,k,j,i,T_i),local(:))/local2



                spinB(:,k,j,i,T_i) = 2.0d0*dot*local(:) - spinB(:,k,j,i,T_i)!スピンの向きをエネルギーを変えずに変更


       !***************************************************************************

  end subroutine over_relaxationB



  subroutine tempereture_exchange(T,T_ex,alpha,E,MC_i)

    integer,intent(inout) :: T_ex(:),T(:)
    integer,intent(in) :: MC_i
    real(8),intent(in) :: E(:),alpha
    real(8) lowtmp,hightmp
    real(8) beta,dE
    integer T_i,j,k,tmp,tmp2,tmp3

    do T_i = 1+mod(MC_i,2),N_T-1,2

       lowtmp = (alpha**(T_i-1))
       hightmp =(alpha**(T_i))
       beta = 1.0d0/T_min*(lowtmp - hightmp)
       dE = E(T_i +1) - E(T_i)

       !****************tmperature exchange******************
       if(grnd() <= exp(-beta*dE)) then !if exchange accepted
          !******only exchange temperature******
          T(T_ex(T_i)) = T_i+1
          T(T_ex(T_i+1)) = T_i
          tmp = T_ex(T_i)
          T_ex(T_i) = T_ex(T_i+1)
          T_ex(T_i+1) = tmp
          !***************************************
       end if

    end do

  end subroutine tempereture_exchange

  subroutine rep_counter(T,T_i,tmp_max,tmp_min,counter)
    integer,intent(in) ::T(:),T_i
    integer,intent(inout) ::counter(:),tmp_max(:),tmp_min(:)

    if(T(T_i) > tmp_max(T_i))tmp_max(T_i)=T(T_i)
    if(T(T_i) < tmp_min(T_i))tmp_min(T_i)=T(T_i)

    if( tmp_max(T_i)-tmp_min(T_i)==N_T-1 ) then
       counter(T_i) = counter(T_i) + 1
       if(T(T_i) == N_T)then
          tmp_max(T_i)=N_T
          tmp_min(T_i)=N_T
       else if(T(T_i) == 1) then
          tmp_max(T_i)=1
          tmp_min(T_i)=1
       end if

    end if

  end subroutine rep_counter


  subroutine HeatBathA(i,j,k,rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
    integer,intent(in) :: T_i
    real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    real(8),intent(in) ::J_ij(:,:,:,:),alpha,randnum(:,:),T(:),rg(:)
    integer,intent(in) :: i,j,k
    integer count
    integer m,n,o,p, n_o
    real(8) nn(3), nndm(3), r, S_tmp(3), H_L, exp2H_L, spin_xy, norm2
    real(8) local(3),kitvec(3,3),kitaev(3)
    real(8) base(4,3),posi(3),beta
    integer a,c,d

    beta = -1.0d0/T(N_T)
    count = 2*(k-1)*L*L+2*(j-1)*L+2*i-1
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
    !*********************Heatbath subA***********************************************


       !**********周期境界条件**********************(詳しくは別プログラム)
       a = mod(i-1+L-1,L)+1
       !******************************************



          !**********周期境界条件**********************(詳しくは別プログラム)
          d = mod(j-1+L-1,L)+1
          !******************************************



             !********周期境界条件*********************

             c = mod(k-1+L-1,L)+1
             !***************************************


             !***********position**************************************
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


             local(1) = local(1)+(spinB(3,k,d,i,T_i)+spinB(2,k,j,a,T_i))*gam
             local(2) = local(2)+(spinB(3,k,j,i,T_i)+spinB(1,k,j,a,T_i))*gam
             local(3) = local(3)+(spinB(2,k,j,i,T_i)+spinB(1,k,d,i,T_i))*gam

             local = local + kitaev*K_0-h0*H

             local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
             local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

             nn(:) = local(:)
             !**********************************************
             norm2 = nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3)
             if(norm2>HB_limit)then
                H_L = beta*sqrt(norm2)              !H_L = beta*J*|H_local|、いまsqrt(nn**2)=<J_ij>*|H_local|
                exp2H_L = exp(-H_L-H_L)                                  !高速化のための置き換え
                spinA(3,k,j,i,T_i) = 1d0 + log(exp2H_L+randnum(count,T_i)*(1d0-exp2H_L))/H_L
                !spinA(3,k,j,i,T_i) = 1d0 + log(exp2H_L+grnd()*(1d0-exp2H_L))/H_L
                spin_xy = sqrt(1d0-spinA(3,k,j,i,T_i)*spinA(3,k,j,i,T_i))

                r = randnum(count+1,T_i)*twopi
                !r = grnd()*twopi
                spinA(1,k,j,i,T_i) = spin_xy*cos(r)
                spinA(2,k,j,i,T_i) = spin_xy*sin(r)
                !***H_local方向をZ軸に取っていたのを系の量子化軸に直す
                S_tmp(:) = true_spin(nn(:),spinA(1,k,j,i,T_i),spinA(2,k,j,i,T_i),spinA(3,k,j,i,T_i))

                spinA(:,k,j,i,T_i) = S_tmp(:)
             else    !局所磁場が小さすぎるときはランダムに更新

                spinA(3,k,j,i,T_i) = grnd()*2 - 1d0
                spin_xy = sqrt(1d0-spinA(3,k,j,i,T_i)*spinA(3,k,j,i,T_i))

                r = grnd()*twopi
                spinA(1,k,j,i,T_i) = spin_xy*cos(r)
                spinA(2,k,j,i,T_i) = spin_xy*sin(r)
             endif


  end subroutine HeatBathA

    subroutine HeatBathB(i,j,k,rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
      integer,intent(in) :: T_i
      real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
      real(8),intent(in) ::J_ij(:,:,:,:),alpha,randnum(:,:),T(:),rg(:)
      integer,intent(in) :: i,j,k
      integer count
      integer m,n,o,p, n_o
      real(8) nn(3), nndm(3), r, S_tmp(3), H_L, exp2H_L, spin_xy, norm2
      real(8) local(3),kitvec(3,3),kitaev(3)
      real(8) base(4,3),posi(3),beta
      integer a,c,d

      beta = -1.0d0/T(N_T)
      count = 2*(k-1)*L*L+2*(j-1)*L+2*i-1+2*L*L*L
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

    !*********************Heatbath subB***********************************************


          !**********周期境界条件**********************(詳しくは別プログラム)
          a = mod(i+1+L-1,L)+1
          !******************************************


             !**********周期境界条件**********************(詳しくは別プログラム)
             d = mod(j+1+L-1,L)+1
             !******************************************


                !********周期境界条件*********************

                c = mod(k+1+L-1,L)+1
                !***************************************

                !***********position**************************************
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

                local(1) = local(1)+(spinA(3,k,d,i,T_i)+spinA(2,k,j,a,T_i))*gam
                local(2) = local(2)+(spinA(3,k,j,i,T_i)+spinA(1,k,j,a,T_i))*gam
                local(3) = local(3)+(spinA(2,k,j,i,T_i)+spinA(1,k,d,i,T_i))*gam

                local = local + kitaev*K_0 -h0*H

                local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
                local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

             nn(:) = local(:)
             !**********************************************
             norm2 = nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3)
             if(norm2>HB_limit)then
                H_L = beta*sqrt(norm2)              !H_L = beta*J*|H_local|、いまsqrt(nn**2)=<J_ij>*|H_local|
                exp2H_L = exp(-H_L-H_L)                                  !高速化のための置き換え
                spinB(3,k,j,i,T_i) = 1d0 + log(exp2H_L+randnum(count,T_i)*(1d0-exp2H_L))/H_L
                !spinB(3,k,j,i,T_i) = 1d0 + log(exp2H_L+grnd()*(1d0-exp2H_L))/H_L
                spin_xy = sqrt(1d0-spinB(3,k,j,i,T_i)*spinB(3,k,j,i,T_i))

                r = randnum(count+1,T_i)*twopi
                !r = grnd()*twopi
                spinB(1,k,j,i,T_i) = spin_xy*cos(r)
                spinB(2,k,j,i,T_i) = spin_xy*sin(r)
                !***H_local方向をZ軸に取っていたのを系の量子化軸に直す
                S_tmp(:) = true_spin(nn(:),spinB(1,k,j,i,T_i),spinB(2,k,j,i,T_i),spinB(3,k,j,i,T_i))

                spinB(:,k,j,i,T_i) = S_tmp(:)
             else    !局所磁場が小さすぎるときはランダムに更新
                spinB(3,k,j,i,T_i) = grnd()*2 - 1d0
                spin_xy = sqrt(1d0-spinB(3,k,j,i,T_i)*spinB(3,k,j,i,T_i))

                r = grnd()*twopi
                spinB(1,k,j,i,T_i) = spin_xy*cos(r)
                spinB(2,k,j,i,T_i) = spin_xy*sin(r)
             endif


  end subroutine HeatBathB

   function true_spin(unnrmlzd_nn,S_x,S_y,S_z) result(S_t)
        real(8), intent(in) :: unnrmlzd_nn(3), S_x, S_y, S_z
        real(8) cos_theta, sin_theta, cos_phi, sin_phi
        real(8) nn(3), S_t(3), TS_tmp

        nn(:) = unnrmlzd_nn(:)/sqrt(unnrmlzd_nn(1)*unnrmlzd_nn(1)&
                                  &+unnrmlzd_nn(2)*unnrmlzd_nn(2)+unnrmlzd_nn(3)*unnrmlzd_nn(3))

        cos_theta = nn(3) ! cos(z軸となす角度) = z成分
        sin_theta = sqrt(1d0 - cos_theta*cos_theta)

        if(cos_theta>-1d0.and.cos_theta<1d0)then
            cos_phi = nn(1)/sin_theta
            sin_phi = nn(2)/sin_theta
        elseif(cos_theta==1d0)then    ! xy平面上でx軸となす角度
            S_t(1) = S_x                        !局所磁場がz軸と平行ならスピンを更新しない
            S_t(2) = S_y
            S_t(3) = S_z
            return
        elseif(cos_theta==-1d0)then             !局所磁場がz軸と反平行ならスピンを反転する
            S_t(1) = -S_x
            S_t(2) = -S_y
            S_t(3) = -S_z
            return
        else
!            write(*,*) unnrmlzd_nn(:), S_x, S_y, S_z
!            write(*,*) nn(:), cos_theta, sin_theta, cos_phi, sin_phi
            stop 'HeatBath:H_localの角度を求める関数のエラー'
        endif

        TS_tmp = S_x*cos_theta+S_z*sin_theta

        S_t(1) = -S_y*sin_phi + TS_tmp*cos_phi
        S_t(2) =  S_y*cos_phi + TS_tmp*sin_phi
        S_t(3) = -S_x*sin_theta + S_z*cos_theta
    end function

subroutine HeatBath(rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
    integer,intent(in) :: T_i
    real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    real(8),intent(in) ::J_ij(:,:,:,:),alpha,randnum(:,:),T(:),rg(:)
    integer m,n,o,p,i, n_o,count
    real(8) nn(3), nndm(3), r, S_tmp(3), H_L, exp2H_L, spin_xy, norm2
    real(8) local(3),kitvec(3,3),kitaev(3)
    real(8) base(4,3),posi(3),beta
    integer j, k
    integer a,c,d

    beta = -1.0d0/T(N_T)
    count = 1
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
    !*********************Heatbath subA***********************************************
    do i = 1,L

       !**********周期境界条件**********************(詳しくは別プログラム)
       a = mod(i-1+L-1,L)+1
       !******************************************

       do j = 1,L

          !**********周期境界条件**********************(詳しくは別プログラム)
          d = mod(j-1+L-1,L)+1
          !******************************************

          do k = 1,L

             !********周期境界条件*********************

             c = mod(k-1+L-1,L)+1
             !***************************************


             !***********position**************************************
             posi = (i-1)*base(1,:) +(j-1)*base(2,:) + (k-1)*base(4,:)
             posi = posi -rg
             posi = posi/dble(L)
             !********************************************************

             kitaev = kitvec(1,:)*dot_product(spinB(:,k,j,i,T_i),kitvec(1,:)) &
                  + kitvec(2,:)*dot_product(spinB(:,k,d,i,T_i),kitvec(2,:))&
                  + kitvec(3,:)*dot_product(spinB(:,k,j,a,T_i),kitvec(3,:))

             local(:) = J_ij(3,k,j,i)*spinB(:,k,j,a,T_i) + J_ij(1,k,j,i)*spinB(:,k,j,i,T_i) +&
                  J_ij(2,k,j,i)*spinB(:,k,d,i,T_i) + Jp*J_0*spinB(:,c,j,i,T_i)

             local(1) = local(1)*XY_0
             local(2) = local(2)*XY_0
             local(3) = local(3)*Z_0

             local = local + kitaev-h0*H

             local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
             local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

             nn(:) = local(:)
             !**********************************************
             norm2 = nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3)
             if(norm2>HB_limit)then
                H_L = beta*sqrt(norm2)              !H_L = beta*J*|H_local|、いまsqrt(nn**2)=<J_ij>*|H_local|
                exp2H_L = exp(-H_L-H_L)                                  !高速化のための置き換え
                spinA(3,k,j,i,T_i) = 1d0 + log(exp2H_L+randnum(count,T_i)*(1d0-exp2H_L))/H_L
                spin_xy = sqrt(1d0-spinA(3,k,j,i,T_i)*spinA(3,k,j,i,T_i))

                r = randnum(count+1,T_i)*twopi
                spinA(1,k,j,i,T_i) = spin_xy*cos(r)
                spinA(2,k,j,i,T_i) = spin_xy*sin(r)
                !***H_local方向をZ軸に取っていたのを系の量子化軸に直す
                S_tmp(:) = true_spin(nn(:),spinA(1,k,j,i,T_i),spinA(2,k,j,i,T_i),spinA(3,k,j,i,T_i))

                spinA(:,k,j,i,T_i) = S_tmp(:)
             else    !局所磁場が小さすぎるときはランダムに更新
                
                spinA(3,k,j,i,T_i) = grnd()*2 - 1d0
                spin_xy = sqrt(1d0-spinA(3,k,j,i,T_i)*spinA(3,k,j,i,T_i))

                r = grnd()*twopi
                spinA(1,k,j,i,T_i) = spin_xy*cos(r)
                spinA(2,k,j,i,T_i) = spin_xy*sin(r)
             endif

             count = count +2
             
          enddo
       enddo
    enddo

    !*********************Heatbath subB***********************************************
   do i = 1,L

          !**********周期境界条件**********************(詳しくは別プログラム)
          a = mod(i+1+L-1,L)+1
          !******************************************

          do j = 1,L

             !**********周期境界条件**********************(詳しくは別プログラム)
             d = mod(j+1+L-1,L)+1
             !******************************************

             do k = 1,L

                !********周期境界条件*********************

                c = mod(k+1+L-1,L)+1
                !***************************************

                !***********position**************************************
                posi = (i-1)*base(1,:) +(j-1)*base(2,:) + base(3,:) + (k-1)*base(4,:)
                posi = posi -rg
                posi = posi/dble(L)
                !********************************************************

                kitaev = kitvec(1,:)*dot_product(spinA(:,k,j,i,T_i),kitvec(1,:)) &
                     + kitvec(2,:)*dot_product(spinA(:,k,d,i,T_i),kitvec(2,:))&
                     + kitvec(3,:)*dot_product(spinA(:,k,j,a,T_i),kitvec(3,:))

                local(:) = J_ij(3,k,j,a)*spinA(:,k,j,a,T_i) + J_ij(1,k,j,i)*spinA(:,k,j,i,T_i) +&
                     J_ij(2,k,d,i)*spinA(:,k,d,i,T_i) + Jp*J_0*spinA(:,c,j,i,T_i)

                local(1) = local(1)*XY_0
                local(2) = local(2)*XY_0
                local(3) = local(3)*Z_0

                local = local + kitaev -h0*H

                local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
                local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

             nn(:) = local(:)
             !**********************************************
             norm2 = nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3)
             if(norm2>HB_limit)then
                H_L = beta*sqrt(norm2)              !H_L = beta*J*|H_local|、いまsqrt(nn**2)=<J_ij>*|H_local|
                exp2H_L = exp(-H_L-H_L)                                  !高速化のための置き換え
                spinB(3,k,j,i,T_i) = 1d0 + log(exp2H_L+randnum(count,T_i)*(1d0-exp2H_L))/H_L
                spin_xy = sqrt(1d0-spinB(3,k,j,i,T_i)*spinB(3,k,j,i,T_i))

                r = randnum(count+1,T_i)*twopi
                spinB(1,k,j,i,T_i) = spin_xy*cos(r)
                spinB(2,k,j,i,T_i) = spin_xy*sin(r)
                !***H_local方向をZ軸に取っていたのを系の量子化軸に直す
                S_tmp(:) = true_spin(nn(:),spinB(1,k,j,i,T_i),spinB(2,k,j,i,T_i),spinB(3,k,j,i,T_i))

                spinB(:,k,j,i,T_i) = S_tmp(:)
             else    !局所磁場が小さすぎるときはランダムに更新
                spinB(3,k,j,i,T_i) = grnd()*2 - 1d0
                spin_xy = sqrt(1d0-spinB(3,k,j,i,T_i)*spinB(3,k,j,i,T_i))

                r = grnd()*twopi
                spinB(1,k,j,i,T_i) = spin_xy*cos(r)
                spinB(2,k,j,i,T_i) = spin_xy*sin(r)
             endif

             count = count + 2
          enddo
       enddo
    enddo
  end subroutine HeatBath

subroutine over_relaxation(rg,spinA,spinB,T_i,J_ij)
    real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
    real(8),intent(in) :: J_ij(:,:,:,:),rg(:)
    integer ,intent(in) :: T_i
    real(8) local(3),local2,dot,kitvec(3,3),kitaev(3)
    real(8)base(4,3),posi(3),tmp(2),tmx,tmy,tmz,P(3)
    integer i, j, k, size1 , size2 , size3,R_i
    integer a,c,d

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

    do R_i = 1,N_relax


       !*********************OR subA***********************************************
       do i = 1,L

          !**********周期境界条件**********************(詳しくは別プログラム)
          a = mod(i-1+L-1,L)+1
          !******************************************

          do j = 1,L

             !**********周期境界条件**********************(詳しくは別プログラム)
             d = mod(j-1+L-1,L)+1
             !******************************************

             do k = 1,L

                !********周期境界条件*********************

                c = mod(k-1+L-1,L)+1
                !***************************************


                !***********position**************************************
                posi = (i-1)*base(1,:) +(j-1)*base(2,:) + (k-1)*base(4,:)
                posi = posi -rg
                posi = posi/dble(L)
                !********************************************************

                kitaev = kitvec(1,:)*dot_product(spinB(:,k,j,i,T_i),kitvec(1,:)) &
                     + kitvec(2,:)*dot_product(spinB(:,k,d,i,T_i),kitvec(2,:))&
                     + kitvec(3,:)*dot_product(spinB(:,k,j,a,T_i),kitvec(3,:))

                local(:) = J_ij(3,k,j,i)*spinB(:,k,j,a,T_i) + J_ij(1,k,j,i)*spinB(:,k,j,i,T_i) +&
                     J_ij(2,k,j,i)*spinB(:,k,d,i,T_i) + Jp*J_0*spinB(:,c,j,i,T_i)

                local(1) = local(1)*XY_0
                local(2) = local(2)*XY_0
                local(3) = local(3)*Z_0


                local = local + kitaev-h0*H

                local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
                local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))

                local2 = dot_product(local,local)
                dot = dot_product(spinA(:,k,j,i,T_i),local(:))/local2



                spinA(:,k,j,i,T_i) = 2.0d0*dot*local(:) - spinA(:,k,j,i,T_i)!スピンの向きをエネルギーを変えずに変更


             end do
          end do
       end do
       !***************************************************************************

       !*********************OR subB***********************************************
       do i = 1,L

          !**********周期境界条件**********************(詳しくは別プログラム)
          a = mod(i+1+L-1,L)+1
          !******************************************

          do j = 1,L

             !**********周期境界条件**********************(詳しくは別プログラム)
             d = mod(j+1+size2-1,size2)+1
             !******************************************

             do k = 1,L

                !********周期境界条件*********************

                c = mod(k+1+L-1,L)+1
                !***************************************

                !***********position**************************************
                posi = (i-1)*base(1,:) +(j-1)*base(2,:) + base(3,:) + (k-1)*base(4,:)
                posi = posi -rg
                posi = posi/dble(L)
                !********************************************************

                kitaev = kitvec(1,:)*dot_product(spinA(:,k,j,i,T_i),kitvec(1,:)) &
                     + kitvec(2,:)*dot_product(spinA(:,k,d,i,T_i),kitvec(2,:))&
                     + kitvec(3,:)*dot_product(spinA(:,k,j,a,T_i),kitvec(3,:))

                local(:) = J_ij(3,k,j,a)*spinA(:,k,j,a,T_i) + J_ij(1,k,j,i)*spinA(:,k,j,i,T_i) +&
                     J_ij(2,k,d,i)*spinA(:,k,d,i,T_i) + Jp*J_0*spinA(:,c,j,i,T_i)

                local(1) = local(1)*XY_0
                local(2) = local(2)*XY_0
                local(3) = local(3)*Z_0

                local = local + kitaev -h0*H

                local(1)=local(1)-(Ht(1)*H(2)*posi(2)-Ht(2)*H(1)*posi(2))
                local(2)=local(2)-(-Ht(1)*H(2)*posi(1)+Ht(2)*H(1)*posi(1))



                local2 = dot_product(local,local)
                dot = dot_product(spinB(:,k,j,i,T_i),local(:))/local2



                spinB(:,k,j,i,T_i) = 2.0d0*dot*local(:) - spinB(:,k,j,i,T_i)!スピンの向きをエネルギーを変えずに変更

             end do
          end do
       end do
       !***************************************************************************

    end do


  end subroutine over_relaxation

  subroutine metropolisA(i,j,k,rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
      integer,intent(in) :: T_i,i,j,k
      real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
      real(8),intent(in) ::J_ij(:,:,:,:),alpha,randnum(:,:),T(:),rg(:)
      real(8) Eold,Enew,dE,temp(3),rpi,r1,r2,r3,tmp
      real(8) tmp1,tmp2
      integer size1 , size2 , size3 ,count

      count = 3*(k-1)*L*L+3*(j-1)*L+3*i-2

      !tmp = 1.0d0/T_min*(alpha**(T(T_i)-1))
      tmp = T(N_T)

               !***********spinAの更新**************
               r1 = randnum(count,T_i)*2.0d0-1.0d0
               r2 = randnum(count+1,T_i)*twopi
               r3 = randnum(count+2,T_i)
               tmp1 = sqrt(1.0d0 - r1*r1)       !乱数で方向を改めてみる
               call local_EA(rg,spinA,spinB,J_ij,T_i,i,j,k,Eold)
               temp(1) = spinA(1,k,j,i,T_i)!前の情報を一時保持
               temp(2) = spinA(2,k,j,i,T_i)
               temp(3) = spinA(3,k,j,i,T_i)
               spinA(1,k,j,i,T_i) = tmp1*cos(r2)
               spinA(2,k,j,i,T_i) = tmp1*sin(r2)
               spinA(3,k,j,i,T_i) = r1
               call local_EA(rg,spinA,spinB,J_ij,T_i,i,j,k,Enew)
               dE = Enew - Eold
               !反転の判定
               if(dE <= 0.0d0) then

               else if (r3 < exp(-dE/tmp)) then

               else !rの確率で以前の状態にする
                  spinA(1,k,j,i,T_i) = temp(1)
                  spinA(2,k,j,i,T_i) = temp(2)
                  spinA(3,k,j,i,T_i) = temp(3)
               end if
               !**********spinAの更新終わり***********

    end subroutine metropolisA

      subroutine metropolisB(i,j,k,rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
        integer,intent(in) :: T_i,i,j,k
        real(8),intent(inout) :: spinA(:,:,:,:,:),spinB(:,:,:,:,:)
        real(8),intent(in) ::J_ij(:,:,:,:),alpha,randnum(:,:),T(:),rg(:)
        real(8) Eold,Enew,dE,temp(3),rpi,r1,r2,r3,tmp
        real(8) tmp1,tmp2
        integer size1 , size2 , size3 ,count

        count = 3*(k-1)*L*L+3*(j-1)*L+3*i-2+3*L*L*L

        !tmp = 1.0d0/T_min*(alpha**(T(T_i)-1))
        tmp = T(N_T)

               !***********spinBの更新**************
               r1 = randnum(count,T_i)*2.0d0-1.0d0
               r2 = randnum(count+1,T_i)*twopi
               r3 = randnum(count+2,T_i)
               tmp1 = sqrt(1.0d0 - r1*r1)
               call local_EB(rg,spinB,spinA,J_ij,T_i,i,j,k,Eold)
               temp(1) = spinB(1,k,j,i,T_i)!前の情報を一時保持
               temp(2) = spinB(2,k,j,i,T_i)
               temp(3) = spinB(3,k,j,i,T_i)
               spinB(1,k,j,i,T_i) = tmp1*cos(r2)
               spinB(2,k,j,i,T_i) = tmp1*sin(r2)
               spinB(3,k,j,i,T_i) = r1
               call local_EB(rg,spinB,spinA,J_ij,T_i,i,j,k,Enew)
               dE = Enew - Eold
               !反転の判定
               if(dE <= 0.0d0) then

               else if (r3 < exp(-dE/tmp)) then

               else !rの確率で以前の状態にする
                  spinB(1,k,j,i,T_i) = temp(1)
                  spinB(2,k,j,i,T_i) = temp(2)
                  spinB(3,k,j,i,T_i) = temp(3)
               end if
               !**********spinBの更新終わり***********


    end subroutine metropolisB

  end module metropolis_method

module boot
  use mtmod
  implicit none
contains

  subroutine boot_strap(dataname,dataname2,n,fn,beg,L,alpha,T_min)!n:N_T,fn:N_sample

    character(80) ,intent(in)::dataname(:),dataname2(:)
    integer,intent(in)::n,fn,beg,L
    real(8),intent(in):: alpha,T_min
    integer i,j,k,seed,fnum,fnum2,T_i,rand
    real(8),allocatable:: temp(:,:),er(:,:),temp2(:,:,:)
    real(8),allocatable :: bsample(:,:,:),sample(:,:,:)
    real(8) ,allocatable :: sgL(:),sgLxx(:),sgLyy(:) ,cgLb(:),cgLs(:),gL(:,:,:)
    real(8) tmpg,T
    character(80)syssize ,fname3
    character(80),allocatable::fname2(:) 


    fnum = 41!ファイル一つの長さ(行の数)
    fnum2 = size(dataname2)
    allocate(temp(n,fnum),er(n,fnum2),fname2(fnum))!データの数(何列かとデータの長さ（何行か)で割付け
    allocate(sample(fn,(fnum-beg+1),n))
    allocate(bsample(2*fn,(fnum-beg+1),n))
    allocate(sgL(n),sgLxx(n),sgLyy(n) ,cgLb(n),cgLs(n),gL(fnum2,2*fn,n))

    sgL(:) = 0.0d0
    sgLxx(:) = 0.0d0
    sgLyy(:) = 0.0d0
    cgLb(:) = 0.0d0
    cgLs(:) = 0.0d0

    do i = 1,fnum2
       fname3 = trim(adjustl(dataname2(i)))//'.d'
       open(30+i,file = fname3,action='write')
    end do

    do i= 1,fn
       open(21, file = dataname(i) ,action ='read')
       do T_i =1,n
          read(21,*)temp(T_i,:)
          sample(i,:,T_i) = temp(T_i,beg:fnum)
       end do
    end do

    close(21)


    bsample = 0.0d0

    do T_i = 1,n
       do i = 1,fn*2

          do j = 1,fn

             

             rand = grnd() * (fn) + 1 

             bsample(i,:,T_i) = bsample(i,:,T_i) + sample(rand,:,T_i)
          end do
           bsample(i,:,T_i) = bsample(i,:,T_i)/dble(fn)
       end do
    end do


    do T_i=1,n
       do i = 1,fn*2
          gL(1,i,T_i) = 1.0d0/(2.0d0*sin(acos(-1.0d0)/L)) * sqrt(abs(bsample(i,5,T_i)/bsample(i,6,T_i)-1.0d0))/dble(L)
          gL(2,i,T_i) = 1.0d0/(2.0d0*sin(acos(-1.0d0)/L)) * sqrt(abs(bsample(i,1,T_i)/bsample(i,7,T_i)-1.0d0))/dble(L)
          gL(3,i,T_i) = 1.0d0/(2.0d0*sin(acos(-1.0d0)/L)) * sqrt(abs(bsample(i,4,T_i)/bsample(i,8,T_i)-1.0d0))/dble(L)
          gL(4,i,T_i) = 1.0d0/(2.0d0*sin(acos(-1.0d0)/L)) * sqrt(abs(bsample(i,9,T_i)/bsample(i,11,T_i)-1.0d0))/dble(L)
          gL(5,i,T_i) = 1.0d0/(2.0d0*sin(acos(-1.0d0)/L)) * sqrt(abs(bsample(i,10,T_i)/bsample(i,12,T_i)-1.0d0))/dble(L)
          
       end do
    end do


    do T_i = 1,n
       do i = 1,fn*2
          sgL(T_i) = sgL(T_i) +  gL(1,i,T_i)
          sgLxx(T_i) = sgLxx(T_i) +  gL(2,i,T_i)
          sgLyy(T_i)= sgLyy(T_i) +  gL(3,i,T_i)
          cgLb(T_i)= cgLb(T_i) +  gL(4,i,T_i)
          cgLs(T_i)= cgLs(T_i) +  gL(5,i,T_i)
       end do
    end do

    sgL = sgL/dble(2*fn)
    sgLxx = sgLxx/dble(2*fn)
    sgLyy = sgLyy/dble(2*fn)
    cgLb = cgLb/dble(2*fn)
    cgLs = cgLs/dble(2*fn)

    er = 0.0d0

    do T_i = 1,n
       do i = 1,fn
          er(T_i,1) = er(T_i,1) + (gL(1,i,T_i)-sgL(T_i))**2
          er(T_i,2) = er(T_i,2) + (gL(2,i,T_i)-sgLxx(T_i))**2
          er(T_i,3) = er(T_i,3) + (gL(3,i,T_i)-sgLyy(T_i))**2
          er(T_i,4) = er(T_i,4) + (gL(4,i,T_i)-cgLb(T_i))**2
          er(T_i,5) = er(T_i,5) + (gL(5,i,T_i)-cgLs(T_i))**2 
       end do
       er(T_i,:) = er(T_i,:)/dble(fn*2-1)
       er(T_i,:) = sqrt(er(T_i,:))
    end do
    

  
       do j=1, n
           T = T_min*(alpha**dble(-j+1))
           write(1+30,*)T,sgL(j),er(j,1)
           write(2+30,*)T,sgLxx(j),er(j,2)
           write(3+30,*)T,sgLyy(j),er(j,3)
           write(4+30,*)T,cgLb(j),er(j,4)
           write(5+30,*)T,cgLs(j),er(j,5)!一列目は温度や位置など基準のものを書く
       end do

    do i = 1,fnum
       close(i+30)
    end do

  end subroutine boot_strap

end module boot



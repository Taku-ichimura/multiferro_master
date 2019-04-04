module toukei
  implicit none
contains
  subroutine statistics(fname,dataname,n,fn)!統計を取る物理量の名前、読み込むファイルの名前、データの量(何行か)、ファイルの個数(seedの数など)
    character(80) ,intent(in)::fname(:),dataname(:)
    integer,intent(in)::n,fn
    integer i,j,k,seed,fnum
    real(8),allocatable:: data(:,:),data2(:,:),temp(:,:),er(:,:)
    real(8) dT,T
    character(80)syssize ,fname3
    character(80),allocatable::fname2(:) 
    

    fnum = size(fname)!ファイル一つの長さ(行の数)
    allocate(data(n,fnum),data2(n,fnum),temp(n,fnum),er(n,fnum),fname2(fnum))!データの数(何列かとデータの長さ（何行か)で割付け

    do i = 1,fnum
      fname3 = trim(adjustl(fname(i)))//'.d'

      ! write(fname2(i),'(a,a,a)')
       open(30+i,file = fname3,action='write')
    end do
    


    data(:,:) =0.0d0
    data2(:,:) = 0.0d0

    do i= 1,fn
       open(21, file = dataname(i) ,action ='read')

       do j = 1,n
          read(21,*)temp(j,:)
       end do

       data = data + temp
       data2 = data2 + temp**2
    end do

    close(21)



    data(:,:) = data(:,:)/dble(fn)
    data2 = data2 /dble(fn)
    er = fn/dble(fn-1)*abs(data2 - data**2)
    er = er/dble(fn)
    er = sqrt(er)

    
    do i = 1, fnum
       do j=1, n
          write(i+30,*)temp(j,1),data(j,i),er(j,i) !一列目は温度や位置など基準のものを書く
       end do
    end do

    do i = 1,fnum
       close(i+30)
    end do

    

    

  end subroutine statistics
end module toukei



     


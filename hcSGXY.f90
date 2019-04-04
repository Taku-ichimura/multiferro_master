!analyzing XY honeycomb lattice of SG phase by tempereture exchange method
program main

  !$ use omp_lib
  use global
  use toukei
  use iniset
  use mtmod
  use boot
  use syscalc
  use metropolis_method
  !  use MPI

  implicit none

  !***********************increment*********************
  integer T_i , S_i !for tempereture and sample
  integer i , j ,kk,ii,r_i,jj
  integer sweep2
  integer count
  integer , parameter :: fo2= 3000,fo3=4000
  !*****************************************************

  !***********************system************************
  real(8),allocatable:: spinA(:,:,:,:,:) ,spinB(:,:,:,:,:) ,spin_repB(:,:,:,:,:), spin_repA(:,:,:,:,:)
  real(8) ,allocatable::J_ij(:,:,:,:)
  real(8) dT ,alpha!delta tmp
  real(8) tmp_i
  real(8) ,allocatable ::randnum(:,:),randnum_rep(:,:)
  real(8) t1,t2,rg(3)
  real(8) cyc,T(N_T)
  integer T_rep(N_T)!temperature point
  integer T_ex(N_T),T_ex_rep(N_T)
  integer counter(N_T),tmp_max(N_T),tmp_min(N_T)
  !******************************************************

  !*******************physical quantiteis****************
  real(8) E(N_T),E_rep(N_T),sum_e(N_T,2)
  real(8) C(N_T),test,dE1,dE2
  real(8) m(N_T,3), sum_rmsm(N_T),rmsm(N_T), sum_rmsm2(N_T),m_rep(N_T,2)
  real(8) AF(N_T,3),xAF(N_T,3),sum_AF(N_T,6),sum_rmsAF(N_T,2)
  real(8) X_SG(N_T,5),sumX_SG(N_T,2),X_SG_kmin(N_T,3)
  real(8) X_CG(N_T,4)
  real(8) k(N_T,4,3),sum_k(N_T,4,3),sum_k_norm(N_T,4) !chirality
  real(8) tm(N_T,3),sum_tm(N_T,3) !troidal moment
  real(8) scP(N_T,3),absscP(N_T,3),sum_scP(N_T,6),sum_scP2(N_T,3) !spin current polarization
  real(8) tmP(N_T,3),abstmP(N_T,3),sum_tmP(N_T,6),sum_tmP2(N_T,3), tmP_rep(N_T,3) !troidal moment polarization
  real(8) scPerm(N_T,3),scPermx,scPermy
  real(8) tmPerm(N_T,3),tmPermx,tmPermy,tmPermz
  real(8) q_uv(N_T,4),q_uv_kmin(N_T,4,2),sum_q_uv_kmin(N_T,4),q2_uv_kmin(N_T)
  real(8) sum_gL(N_T,3),g_L(N_T,3) !corelation-ratio
  real(8) q_k(N_T,2),q_k_kmin(N_T,2,2),sum_q_k(N_T,2),sum_q_k_kmin(N_T,2)
  real(8) mono(N_T),sum_mono(N_T), quad(N_T,9),sum_quad(N_T,9)
  real(8) sum_vmAF(N_T,3),sum_vm(N_T,3),sus
  real(8),allocatable ::spinA0(:,:,:,:) ,spinB0(:,:,:,:)
  real(8) m0(2),spin_ac,toro_ac(4)
  real(8) zig(N_T,3,3),sum_zig1(N_T,4),sum_zig2(N_T,4),sum_zig3(N_T,4),sum_zig(N_T,4)
  real(8) stripe(N_T,3,3),sum_stripe1(N_T,4),sum_stripe2(N_T,4),sum_stripe3(N_T,4),sum_stripe(N_T,4)
  real(8) pm(N_T,3),sum_pm(N_T,4),sum_m(N_T,3)
  real(8) scP_intr(N_T,3),absscP_intr(N_T,3),sum_scP_intr(N_T,6),sum_scP2_intr(N_T,3) !spin current polarization
  real(8) scP_para(N_T,3),absscP_para(N_T,3),sum_scP_para(N_T,6),sum_scP2_para(N_T,3) !spin current polarization
  real(8) scPerm_para(N_T,3),scPerm_intr(N_T,3)
  real(8) scP_bond(N_T,3,3),absscP_bond(N_T,3,3),sum_scP_bond(N_T,3,3),sum_scP2_bond(N_T,3)
  !*******************************************************

  !***********************spin space***************************
  real(8) Sm(N_T,3), sum_Srmsm(N_T),Srmsm(N_T), sum_Srmsm2(N_T)
  real(8) Spm(N_T,3),sum_Spm(N_T,4),sum_Sm(N_T,3)
  real(8) SAF(N_T,3),SxAF(N_T,3),sum_SAF(N_T,6),sum_SrmsAF(N_T,2)
  real(8) Sk(N_T,4,3),sum_Sk(N_T,4,3),sum_Sk_norm(N_T,4) !chirality
  real(8) Stm(N_T,3),sum_Stm(N_T,3) !troidal moment
  real(8) SscP(N_T,3),SabsscP(N_T,3),sum_SscP(N_T,6),sum_SscP2(N_T,3) !spin current polarization
  real(8) StmP(N_T,3),SabstmP(N_T,3),sum_StmP(N_T,6),sum_StmP2(N_T,3) !troidal moment polarization
  real(8) SscPerm(N_T,3),SscPermx,SscPermy
  real(8) StmPerm(N_T,3),StmPermx,StmPermy,StmPermz
  real(8) Szig(N_T,3,3),sum_Szig1(N_T,4),sum_Szig2(N_T,4),sum_Szig3(N_T,4),sum_Szig(N_T,4)
  real(8) Sstripe(N_T,3,3),sum_Sstripe1(N_T,4),sum_Sstripe2(N_T,4),sum_Sstripe3(N_T,4),sum_Sstripe(N_T,4)
  real(8) SscP_intr(N_T,3),SabsscP_intr(N_T,3),sum_SscP_intr(N_T,6),sum_SscP2_intr(N_T,3) !spin current polarization
  real(8) SscP_para(N_T,3),SabsscP_para(N_T,3),sum_SscP_para(N_T,6),sum_SscP2_para(N_T,3) !spin current polarization
  real(8) SscPerm_para(N_T,3),SscPerm_intr(N_T,3)
  !*************************file**************************
  character(80) fname(N_sample), dataname(41),dataname2(5)
  character(80) syssize,dataname3(32),fname2(N_sample)
  character(80) fname3(N_sample),fname4(N_sample),dataname4(65)
  character(80) dataname5(82),fname5(N_sample)
  !*******************************************************



  !************************data name**************************************
  write(syssize,'(a,i2.2,a,i2.2,a,i2.2,a)')'L',L,'x',int(H(1)),'y',int(H(2)),'K'
  dataname(1) = 'T'//syssize
  dataname(2) = 'e'//syssize
  dataname(3) = 'C'//syssize
  dataname(4) = 'rmsm'//syssize
  dataname(5) = 'k1_norm'//syssize
  dataname(6) = 'k2_norm'//syssize
  dataname(7) = 'k3_norm'//syssize
  dataname(8) = 'scPa'//syssize
  dataname(9) = 'scPb'//syssize
  dataname(10) = 'absscPa'//syssize
  dataname(11) = 'absscPb'//syssize
  dataname(12) = 'xAFa'//syssize
  dataname(13) = 'xAFb'//syssize
  dataname(14) = 'mAFa'//syssize
  dataname(15) = 'mAFb'//syssize
  dataname(16) = 'ta'//syssize
  dataname(17) = 'tb'//syssize
  dataname(18) = 'tc'//syssize
  dataname(19) = 'tmPa'//syssize
  dataname(20) = 'tmPb'//syssize
  dataname(21) = 'tmPc'//syssize
  dataname(22) = 'abstmPa'//syssize
  dataname(23) = 'abstmPb'//syssize
  dataname(24) = 'abstmPc'//syssize
  dataname(25) = 'scperma'//syssize
  dataname(26) = 'scpermb'//syssize
  dataname(27) = 'tmperma'//syssize
  dataname(28) = 'tmpermb'//syssize
  dataname(29) = 'tmpermc'//syssize
  dataname(30) = 'q_0_xx'//syssize
  dataname(31) = 'q_0_xy'//syssize
  dataname(32) = 'q_0_yx'//syssize
  dataname(33) = 'q_0_yy'//syssize
  dataname(34) = 'q_0'//syssize
  dataname(35) = 'q_kmin'//syssize
  dataname(36) = 'q_kmin_xx'//syssize
  dataname(37) = 'q_kmin_yy'//syssize
  dataname(38) = 'q_k_b'//syssize !b:bond
  dataname(39) = 'q_k_s'//syssize !s:site
  dataname(40) = 'q_k_kmin_b'//syssize
  dataname(41) = 'q_k_kmin_s'//syssize

  dataname2(1) = 'gL'//syssize
  dataname2(2) = 'gLxx'//syssize
  dataname2(3) = 'gLyy'//syssize
  dataname2(4) = 'cgLb'//syssize
  dataname2(5) = 'cgLs'//syssize

  dataname3(1) = 'T'//syssize
  dataname3(2) = 'mono'//syssize
  dataname3(3) = 'quad11'//syssize
  dataname3(4) = 'quad12'//syssize
  dataname3(5) = 'quad13'//syssize
  dataname3(6) = 'quad21'//syssize
  dataname3(7) = 'quad22'//syssize
  dataname3(8) = 'quad23'//syssize
  dataname3(9) = 'quad31'//syssize
  dataname3(10) = 'quad32'//syssize
  dataname3(11) = 'quad33'//syssize
  dataname3(12) = 'vAFx'//syssize
  dataname3(13) = 'vAFy'//syssize
  dataname3(14) = 'vmx'//syssize
  dataname3(15) = 'vmy'//syssize
  dataname3(16) = 'sus'//syssize
  dataname3(17) = 'k1_a'//syssize
  dataname3(18) = 'k1_b'//syssize
  dataname3(19) = 'k1_c'//syssize
  dataname3(20) = 'k2_a'//syssize
  dataname3(21) = 'k2_b'//syssize
  dataname3(22) = 'k2_c'//syssize
  dataname3(23) = 'k3_a'//syssize
  dataname3(24) = 'k3_b'//syssize
  dataname3(25) = 'k3_c'//syssize
  dataname3(26) = 'scPc'//syssize
  dataname3(27) = 'absscPc'//syssize
  dataname3(28) = 'scpermc'//syssize
  dataname3(29) = 'mAFc'//syssize
  dataname3(30) = 'xAFc'//syssize
  dataname3(31) = 'mAF'//syssize
  dataname3(32) = 'susAF'//syssize

  dataname4(1) = 'T'//syssize
  dataname4(2) = 'ziga'//syssize
  dataname4(3) = 'zigb'//syssize
  dataname4(4) = 'zigc'//syssize
  dataname4(5) = 'zig'//syssize
  dataname4(6) = 'zig1a'//syssize
  dataname4(7) = 'zig1b'//syssize
  dataname4(8) = 'zig1c'//syssize
  dataname4(9) = 'zig1'//syssize
  dataname4(10) = 'zig2a'//syssize
  dataname4(11) = 'zig2b'//syssize
  dataname4(12) = 'zig2c'//syssize
  dataname4(13) = 'zig2'//syssize
  dataname4(14) = 'zig3a'//syssize
  dataname4(15) = 'zig3b'//syssize
  dataname4(16) = 'zig3c'//syssize
  dataname4(17) = 'zig3'//syssize
  dataname4(18) = 'stripea'//syssize
  dataname4(19) = 'stripeb'//syssize
  dataname4(20) = 'stripec'//syssize
  dataname4(21) = 'stripe'//syssize
  dataname4(22) = 'stripe1a'//syssize
  dataname4(23) = 'stripe1b'//syssize
  dataname4(24) = 'stripe1c'//syssize
  dataname4(25) = 'stripe1'//syssize
  dataname4(26) = 'stripe2a'//syssize
  dataname4(27) = 'stripe2b'//syssize
  dataname4(28) = 'stripe2c'//syssize
  dataname4(29) = 'stripe2'//syssize
  dataname4(30) = 'stripe3a'//syssize
  dataname4(31) = 'stripe3b'//syssize
  dataname4(32) = 'stripe3c'//syssize
  dataname4(33) = 'stripe3'//syssize
  dataname4(34) = 'pma'//syssize
  dataname4(35) = 'pmb'//syssize
  dataname4(36) = 'pmc'//syssize
  dataname4(37) = 'pm'//syssize
  dataname4(38) = 'ma'//syssize
  dataname4(39) = 'mb'//syssize
  dataname4(40) = 'mc'//syssize
  dataname4(41) = 'k4_a'//syssize
  dataname4(42) = 'k4_b'//syssize
  dataname4(43) = 'k4_c'//syssize
  dataname4(44) = 'k4_norm'//syssize
  dataname4(45) = 'scPa_intr'//syssize
  dataname4(46) = 'scPb_intr'//syssize
  dataname4(47) = 'scPc_intr'//syssize
  dataname4(48) = 'scPa_para'//syssize
  dataname4(49) = 'scPb_para'//syssize
  dataname4(50) = 'scPc_para'//syssize
  dataname4(51) = 'scintrperma'//syssize
  dataname4(52) = 'scintrpermb'//syssize
  dataname4(53) = 'scintrpermc'//syssize
  dataname4(54) = 'scparaperma'//syssize
  dataname4(55) = 'scparapermb'//syssize
  dataname4(56) = 'scparapermc'//syssize
  dataname4(57) = 'scPa_bond1'//syssize
  dataname4(58) = 'scPb_bond1'//syssize
  dataname4(59) = 'scPc_bond1'//syssize
  dataname4(60) = 'scPa_bond2'//syssize
  dataname4(61) = 'scPb_bond2'//syssize
  dataname4(62) = 'scPc_bond2'//syssize
  dataname4(63) = 'scPa_bond3'//syssize
  dataname4(64) = 'scPb_bond3'//syssize
  dataname4(65) = 'scPc_bond3'//syssize

  dataname5(1) = 'spinspace/T'//syssize
  dataname5(2) = 'spinspace/zigx'//syssize
  dataname5(3) = 'spinspace/zigy'//syssize
  dataname5(4) = 'spinspace/zigz'//syssize
  dataname5(5) = 'spinspace/zig1x'//syssize
  dataname5(6) = 'spinspace/zig1y'//syssize
  dataname5(7) = 'spinspace/zig1z'//syssize
  dataname5(8) = 'spinspace/zig2x'//syssize
  dataname5(9) = 'spinspace/zig2y'//syssize
  dataname5(10) = 'spinspace/zig2z'//syssize
  dataname5(11) = 'spinspace/zig3x'//syssize
  dataname5(12) = 'spinspace/zig3y'//syssize
  dataname5(13) = 'spinspace/zig3z'//syssize
  dataname5(14) = 'spinspace/stripex'//syssize
  dataname5(15) = 'spinspace/stripey'//syssize
  dataname5(16) = 'spinspace/stripez'//syssize
  dataname5(17) = 'spinspace/stripe1x'//syssize
  dataname5(18) = 'spinspace/stripe1y'//syssize
  dataname5(19) = 'spinspace/stripe1z'//syssize
  dataname5(20) = 'spinspace/stripe2x'//syssize
  dataname5(21) = 'spinspace/stripe2y'//syssize
  dataname5(22) = 'spinspace/stripe2z'//syssize
  dataname5(23) = 'spinspace/stripe3x'//syssize
  dataname5(24) = 'spinspace/stripe3y'//syssize
  dataname5(25) = 'spinspace/stripe3z'//syssize
  dataname5(26) = 'spinspace/k1_x'//syssize
  dataname5(27) = 'spinspace/k1_y'//syssize
  dataname5(28) = 'spinspace/k1_z'//syssize
  dataname5(29) = 'spinspace/k2_x'//syssize
  dataname5(30) = 'spinspace/k2_y'//syssize
  dataname5(31) = 'spinspace/k2_z'//syssize
  dataname5(32) = 'spinspace/k3_x'//syssize
  dataname5(33) = 'spinspace/k3_y'//syssize
  dataname5(34) = 'spinspace/k3_z'//syssize
  dataname5(35) = 'spinspace/scPx'//syssize
  dataname5(36) = 'spinspace/scPy'//syssize
  dataname5(37) = 'spinspace/scPz'//syssize
  dataname5(38) = 'spinspace/absscPx'//syssize
  dataname5(39) = 'spinspace/absscPy'//syssize
  dataname5(40) = 'spinspace/absscPz'//syssize
  dataname5(41) = 'spinspace/xAFx'//syssize
  dataname5(42) = 'spinspace/xAFy'//syssize
  dataname5(43) = 'spinspace/xAFz'//syssize
  dataname5(44) = 'spinspace/mAFx'//syssize
  dataname5(45) = 'spinspace/mAFy'//syssize
  dataname5(46) = 'spinspace/mAFz'//syssize
  dataname5(47) = 'spinspace/tx'//syssize
  dataname5(48) = 'spinspace/ty'//syssize
  dataname5(49) = 'spinspace/tz'//syssize
  dataname5(50) = 'spinspace/tmPx'//syssize
  dataname5(51) = 'spinspace/tmPy'//syssize
  dataname5(52) = 'spinspace/tmPz'//syssize
  dataname5(53) = 'spinspace/abstmPx'//syssize
  dataname5(54) = 'spinspace/abstmPy'//syssize
  dataname5(55) = 'spinspace/abstmPz'//syssize
  dataname5(56) = 'spinspace/scpermx'//syssize
  dataname5(57) = 'spinspace/scpermy'//syssize
  dataname5(58) = 'spinspace/scpermz'//syssize
  dataname5(59) = 'spinspace/tmpermx'//syssize
  dataname5(60) = 'spinspace/tmpermy'//syssize
  dataname5(61) = 'spinspace/tmpermz'//syssize
  dataname5(62) = 'spinspace/pmx'//syssize
  dataname5(63) = 'spinspace/pmy'//syssize
  dataname5(64) = 'spinspace/pmz'//syssize
  dataname5(65) = 'spinspace/mx'//syssize
  dataname5(66) = 'spinspace/my'//syssize
  dataname5(67) = 'spinspace/mz'//syssize
  dataname5(68) = 'spinspace/k4_x'//syssize
  dataname5(69) = 'spinspace/k4_y'//syssize
  dataname5(70) = 'spinspace/k4_z'//syssize
  dataname5(71) = 'spinspace/scPa_intr'//syssize
  dataname5(72) = 'spinspace/scPb_intr'//syssize
  dataname5(73) = 'spinspace/scPc_intr'//syssize
  dataname5(74) = 'spinspace/scPa_para'//syssize
  dataname5(75) = 'spinspace/scPb_para'//syssize
  dataname5(76) = 'spinspace/scPc_para'//syssize
  dataname5(77) = 'spinspace/scintrperma'//syssize
  dataname5(78) = 'spinspace/scintrpermb'//syssize
  dataname5(79) = 'spinspace/scintrpermc'//syssize
  dataname5(80) = 'spinspace/scparaperma'//syssize
  dataname5(81) = 'spinspace/scparapermb'//syssize
  dataname5(82) = 'spinspace/scparapermc'//syssize

  !****************************************************************************

  ! if(MPIswitch) then
  !   !******************MPI*************************
  !  integer::isize,irank
  ! integer::ierr,mstaus(MPI_STATUS_SIZE)

  !call mpi_init(ierr)
  !call mpi_comm_size(MPI_COMM_WORLD,isize,ierr)
  !call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr)
  !!**********************************************
  !end if

  !***********************sample loop**********************
  do seed = 1,N_sample


     allocate(spinA(3,L,L,L,N_T))
     allocate(spinB(3,L,L,L,N_T))
     !allocate(spin_repA(2,L,L,L,N_T))
     !allocate(spin_repB(2,L,L,L,N_T))
     allocate(J_ij(3,L,L,L))
     allocate(randnum(N_spin*2,N_T))
     !allocate(randnum_rep(N_spin*2,N_T))

     fo = seed
     H = H_ini
     Ht = Ht_ini
     h0 =1.0d0

     call r_gcalc(rg)

     !************initialize**************
     !call iniJ_ij(J_ij,seed)
     !if(nitiswitch)call iniJ_ij_niti(J_ij,seed)
     call initiarize_spin(spinA,spinB,seed)
     !call initiarize_spin_rep(spin_repA,spin_repB)
     call initiarize_T(T_ex,T_ex_rep)
     !     alpha = (T_min/T_max)**(1.0d0/dble(N_T-1))
     dT = (T_max - T_min)/dble(mesh-1)
     !************************************
     J_ij = regJ


     write(fname(seed),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a)') 'hcxySGL',L,'x',int(H(1)),'y',int(H(2)),'seed',seed,'K.d' !output format(system size,magnetic field,seed)
     fname(seed) =  adjustl(fname(seed))
     open(seed, file = fname(seed) ,action ='write')

     write(fname2(seed),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a)') 'hcxySG2L',L,'x',int(H(1)),'y',int(H(2)),'seed',seed,'K.d' !output format(system size,magnetic field,seed)
     fname2(seed) =  adjustl(fname2(seed))
     open(fo2, file = fname2(seed) ,action ='write')

     write(fname3(seed),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a)') 'hcxySG3L',L,'x',int(H(1)),'y',int(H(2)),'seed',seed,'K.d' !output format(system size,magnetic field,seed)
     fname3(seed) =  adjustl(fname3(seed))
     open(fo3, file = fname3(seed) ,action ='write')

     write(fname5(seed),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a)') 'hcxySG4L',L,'x',int(H(1)),'y',int(H(2)),'seed',seed,'K.d' !output format(system size,magnetic field,seed)
     fname5(seed) =  adjustl(fname5(seed))
     open(fo3+1000, file = fname5(seed) ,action ='write')

     do S_i=1,mesh

        T(N_T) = -(S_i-1)*dT + T_max

        !*******************zero shift****************
        sweep2 = 0
        tmp_max = 1
        tmp_min = N_T
        counter = 0
        sum_e =0.0d0
        sum_rmsm = 0.0d0
        sum_rmsAF = 0.0d0
        sum_AF = 0.0d0
        sum_k = 0.0d0
        sum_tm = 0.0d0
        sum_scP = 0.0d0
        sum_scP2 = 0.0d0
        sum_scP_intr = 0.0d0
        sum_scP2_intr = 0.0d0
        sum_scP_para = 0.0d0
        sum_scP2_para = 0.0d0
        sum_tmP = 0.0d0
        sum_tmP2 = 0.0d0
        sum_q_uv_kmin = 0.0d0
        sumX_SG = 0.0d0
        X_SG = 0.0d0
        X_CG = 0.0d0
        sum_q_k = 0.0d0
        sum_q_k_kmin = 0.0d0
        q_uv = 0.0d0
        q_uv_kmin = 0.0d0
        sum_rmsm2 = 0.0d0
        sum_mono=0.0d0
        sum_quad=0.0d0
        sum_vmAF=0.0d0
        sum_vm=0.0d0
        sum_k_norm=0.0d0
        sum_zig1=0.0d0
        sum_zig2=0.0d0
        sum_zig3=0.0d0
        sum_zig=0.0d0
        sum_stripe1=0.0d0
        sum_stripe2=0.0d0
        sum_stripe3=0.0d0
        sum_stripe=0.0d0
        sum_pm=0.0d0
        sum_m=0.0d0
        sum_scP_bond=0.0d0

        sum_Srmsm = 0.0d0
        sum_SrmsAF = 0.0d0
        sum_SAF = 0.0d0
        sum_Sk = 0.0d0
        sum_Stm = 0.0d0
        sum_SscP = 0.0d0
        sum_SscP2 = 0.0d0
        sum_SscP_intr = 0.0d0
        sum_SscP2_intr = 0.0d0
        sum_SscP_para = 0.0d0
        sum_SscP2_para = 0.0d0
        sum_StmP = 0.0d0
        sum_StmP2 = 0.0d0
        sum_Srmsm2 = 0.0d0
        sum_Szig1=0.0d0
        sum_Szig2=0.0d0
        sum_Szig3=0.0d0
        sum_Szig=0.0d0
        sum_Sstripe1=0.0d0
        sum_Sstripe2=0.0d0
        sum_Sstripe3=0.0d0
        sum_Sstripe=0.0d0
        sum_Spm=0.0d0
        sum_Sm=0.0d0
        !*********************************************

        H = H_ini
        Ht = Ht_ini


        !************************MC sweep********************************
        do i = 1,N_sweep

           if(i > N_sweep/2)then !dvide by N_sweep/2

              !**zero shift**   to inistiarize observals
              m = 0.0d0
              AF = 0.0d0
              k = 0.0d0
              scP = 0.0d0
              scP_intr = 0.0d0
              scP_para = 0.0d0
              tm = 0.0d0
              tmP = 0.0d0
              q_uv = 0.0d0
              q_uv_kmin = 0.0d0
              q_k = 0.0d0
              q_k_kmin = 0.0d0
              mono=0.0d0
              quad=0.0d0
              zig=0.0d0
              stripe=0.0d0
              pm=0.0d0
              scP_bond=0.0d0

              Sm = 0.0d0
              SAF = 0.0d0
              Sk = 0.0d0
              SscP = 0.0d0
              SscP_intr = 0.0d0
              SscP_para = 0.0d0
              Stm = 0.0d0
              StmP = 0.0d0
              Szig=0.0d0
              Sstripe=0.0d0
              Spm=0.0d0
              !*************

              sweep2 = sweep2 + 1

              !$omp parallel do num_threads(N_core)
              do T_i = 1,N_T
                 call calc_vec(spinA,spinB,T_i,T_ex,m,Sm)
                 call sysE(rg,spinA,spinB,T_i,T_ex,J_ij,E)
                 call staggered(spinA,spinB,T_i,T_ex,AF,SAF)
                 call calc_chila_pola(spinA,spinB,T_i,T_ex,k,scP,Sk,SscP)
                 call scP_int_heisen(spinA,spinB,T_i,T_ex,scP_intr,SscP_intr)
                 call calc_chila_pola_bond(spinA,spinB,T_i,T_ex,scP_bond)
                 ! call spin_overlap(spinA,spin_repA,spinB,spin_repB,T_i,T_ex,T_ex_rep,q_uv,q_uv_kmin)
                 !call chiral_overlap(spinA,spin_repA,spinB,spin_repB,T_i,T_ex,T_ex_rep,q_k,q_k_kmin)
                 call troidal(rg,spinA,spinB,T_i,T_ex,m,tm,tmP,Stm,StmP)
                 call monpole(rg,spinA,spinB,T_i,T_ex,m,mono)
                 call quadpole(rg,spinA,spinB,T_i,T_ex,m,quad)
                 call zigzag(spinA,spinB,T_i,T_ex,zig,Szig)
                 call stripy(spinA,spinB,T_i,T_ex,stripe,Sstripe)
                 call plane_calc_vec(spinA,spinB,T_i,T_ex,pm,Spm)
              end do
              !$omp end parallel do


          !write(*,*)E/dble(N_spin)

              sum_e(1:N_T,1) = sum_e(1:N_T,1) + E(1:N_T)/dble(N_spin)
              sum_e(1:N_T,2) = sum_e(1:N_T,2) + E(1:N_T)/dble(N_spin)*E/dble(N_spin)

              do T_i = 1,N_T
                 sum_rmsm(T_i) = dot_product(m(T_i,:),m(T_i,:)) + sum_rmsm(T_i)
                 sum_rmsm2(T_i) =  sum_rmsm2(T_i) + sqrt(dot_product(m(T_i,:),m(T_i,:)))

                 sum_rmsAF(T_i,1) = dot_product(AF(T_i,:),AF(T_i,:)) + sum_rmsAF(T_i,1)
                 sum_rmsAF(T_i,2) =  sum_rmsAF(T_i,2) + sqrt(dot_product(AF(T_i,:),AF(T_i,:)))

                 sum_pm(T_i,4) = dot_product(pm(T_i,:),pm(T_i,:)) + sum_pm(T_i,4)

                 sum_Srmsm(T_i) = dot_product(Sm(T_i,:),Sm(T_i,:)) + sum_Srmsm(T_i)
                 sum_Srmsm2(T_i) =  sum_Srmsm2(T_i) + sqrt(dot_product(Sm(T_i,:),Sm(T_i,:)))

                 sum_SrmsAF(T_i,1) = dot_product(SAF(T_i,:),SAF(T_i,:)) + sum_SrmsAF(T_i,1)
                 sum_SrmsAF(T_i,2) =  sum_SrmsAF(T_i,2) + sqrt(dot_product(AF(T_i,:),AF(T_i,:)))

                 sum_Spm(T_i,4) = dot_product(Spm(T_i,:),Spm(T_i,:)) + sum_Spm(T_i,4)

              end do

              sum_m(:,1) =  sum_m(:,1) + (m(:,1))**2!staggered_x^2
              sum_m(:,2) =  sum_m(:,2) + (m(:,2))**2!staggered_y^2
              sum_m(:,3) =  sum_m(:,3) + (m(:,3))**2!staggered_z^2
              sum_pm(:,1) =  sum_pm(:,1) + (pm(:,1))**2!staggered_x^2
              sum_pm(:,2) =  sum_pm(:,2) + (pm(:,2))**2!staggered_y^2
              sum_pm(:,3) =  sum_pm(:,3) + (pm(:,3))**2!staggered_z^2
              sum_k(:,1,:) = sum_k(:,1,:) + (k(:,1,:)/dble(N_spin))**2!chiral-bond1
              sum_k(:,2,:) = sum_k(:,2,:) + (k(:,2,:)/dble(N_spin))**2!chiral-bond2
              sum_k(:,3,:) = sum_k(:,3,:) +(k(:,3,:)/dble(N_spin))**2!chiral-bond3
              sum_k(:,4,:) = sum_k(:,3,:) + (k(:,4,:)/dble(N_spin))**2!chiral-inter
              sum_AF(:,1) =  sum_AF(:,1) + abs(AF(:,1))!staggered_x
              sum_AF(:,2) =  sum_AF(:,2) + abs(AF(:,2))!staggered_y
              sum_AF(:,3) =  sum_AF(:,3) + abs(AF(:,3))!staggered_z
              sum_AF(:,4) =  sum_AF(:,4) + (AF(:,1))**2!staggered_x^2
              sum_AF(:,5) =  sum_AF(:,5) + (AF(:,2))**2!staggered_y^2
              sum_AF(:,6) =  sum_AF(:,6) + (AF(:,3))**2!staggered_z^2
              sum_tm(:,1) = sum_tm(:,1) + (tm(:,1)/dble(N_spin))**2!troidal-moment_x
              sum_tm(:,2) = sum_tm(:,2) + (tm(:,2)/dble(N_spin))**2!troidal-moment_y
              sum_tm(:,3) = sum_tm(:,3) + (tm(:,3)/dble(N_spin))**2!troidal-moment_z

              sum_Sm(:,1) =  sum_Sm(:,1) + (Sm(:,1))**2!staggered_x^2
              sum_Sm(:,2) =  sum_Sm(:,2) + (Sm(:,2))**2!staggered_y^2
              sum_Sm(:,3) =  sum_Sm(:,3) + (Sm(:,3))**2!staggered_z^2
              sum_Spm(:,1) =  sum_Spm(:,1) + (Spm(:,1))**2!staggered_x^2
              sum_Spm(:,2) =  sum_Spm(:,2) + (Spm(:,2))**2!staggered_y^2
              sum_Spm(:,3) =  sum_Spm(:,3) + (Spm(:,3))**2!staggered_z^2
              sum_Sk(:,1,:) = sum_Sk(:,1,:) + Sk(:,1,:)/dble(N_spin)!chiral-bond1
              sum_Sk(:,2,:) = sum_Sk(:,2,:) + Sk(:,2,:)/dble(N_spin)!chiral-bond2
              sum_Sk(:,3,:) = sum_Sk(:,3,:) + Sk(:,3,:)/dble(N_spin)!chiral-bond3
              sum_Sk(:,4,:) = sum_Sk(:,3,:) + Sk(:,4,:)/dble(N_spin)!chiral-inter
              sum_SAF(:,1) =  sum_SAF(:,1) + abs(SAF(:,1))!staggered_x
              sum_SAF(:,2) =  sum_SAF(:,2) + abs(SAF(:,2))!staggered_y
              sum_SAF(:,3) =  sum_SAF(:,3) + abs(SAF(:,3))!staggered_z
              sum_SAF(:,4) =  sum_SAF(:,4) + (SAF(:,1))**2!staggered_x^2
              sum_SAF(:,5) =  sum_SAF(:,5) + (SAF(:,2))**2!staggered_y^2
              sum_SAF(:,6) =  sum_SAF(:,6) + (SAF(:,3))**2!staggered_z^2
              sum_Stm(:,1) = sum_Stm(:,1) + (Stm(:,1)/dble(N_spin))**2!troidal-moment_x
              sum_Stm(:,2) = sum_Stm(:,2) + (Stm(:,2)/dble(N_spin))**2!troidal-moment_y
              sum_Stm(:,3) = sum_Stm(:,3) + (Stm(:,3)/dble(N_spin))**2!troidal-moment_z

              !*************zigzag and stripe********************
              sum_zig(:,1) = (zig(:,1,1)**2 + zig(:,1,2)**2 + zig(:,1,3)**2) + sum_zig(:,1)
              sum_zig(:,2) = (zig(:,2,1)**2 + zig(:,2,2)**2 + zig(:,2,3)**2) + sum_zig(:,2)
              sum_zig(:,3) = (zig(:,3,1)**2 + zig(:,3,2)**2 + zig(:,3,3)**2) + sum_zig(:,3)

              sum_zig1(:,1) = zig(:,1,1)**2 + sum_zig1(:,1)
              sum_zig1(:,2) = zig(:,2,1)**2 + sum_zig1(:,2)
              sum_zig1(:,3) = zig(:,3,1)**2 + sum_zig1(:,3)

              sum_zig2(:,1) = zig(:,1,2)**2 + sum_zig2(:,1)
              sum_zig2(:,2) = zig(:,2,2)**2 + sum_zig2(:,2)
              sum_zig2(:,3) = zig(:,3,2)**2 + sum_zig2(:,3)

              sum_zig3(:,1) = zig(:,1,3)**2 + sum_zig3(:,1)
              sum_zig3(:,2) = zig(:,2,3)**2 + sum_zig3(:,2)
              sum_zig3(:,3) = zig(:,3,3)**2 + sum_zig3(:,3)

              sum_stripe(:,1) = stripe(:,1,1)**2 + stripe(:,1,2)**2 + stripe(:,1,3)**2 + sum_stripe(:,1)
              sum_stripe(:,2) = stripe(:,2,1)**2 + stripe(:,2,2)**2 + stripe(:,2,3)**2 + sum_stripe(:,2)
              sum_stripe(:,3) = stripe(:,3,1)**2 + stripe(:,3,2)**2 + stripe(:,3,3)**2 + sum_stripe(:,3)

              sum_stripe1(:,1) = stripe(:,1,1)**2 + sum_stripe1(:,1)
              sum_stripe1(:,2) = stripe(:,2,1)**2 + sum_stripe1(:,2)
              sum_stripe1(:,3) = stripe(:,3,1)**2 + sum_stripe1(:,3)

              sum_stripe2(:,1) = stripe(:,1,2)**2 + sum_stripe2(:,1)
              sum_stripe2(:,2) = stripe(:,2,2)**2 + sum_stripe2(:,2)
              sum_stripe2(:,3) = stripe(:,3,2)**2 + sum_stripe2(:,3)

              sum_stripe3(:,1) = stripe(:,1,3)**2 + sum_stripe3(:,1)
              sum_stripe3(:,2) = stripe(:,2,3)**2 + sum_stripe3(:,2)
              sum_stripe3(:,3) = stripe(:,3,3)**2 + sum_stripe3(:,3)

              !**************************************************

              !*************spin space zigzag and stripe********************
              sum_Szig(:,1) = (Szig(:,1,1)**2 + Szig(:,1,2)**2 + Szig(:,1,3)**2) + sum_Szig(:,1)
              sum_Szig(:,2) = (Szig(:,2,1)**2 + Szig(:,2,2)**2 + Szig(:,2,3)**2) + sum_Szig(:,2)
              sum_Szig(:,3) = (Szig(:,3,1)**2 + Szig(:,3,2)**2 + Szig(:,3,3)**2) + sum_Szig(:,3)

              sum_Szig1(:,1) = Szig(:,1,1)**2 + sum_Szig1(:,1)
              sum_Szig1(:,2) = Szig(:,2,1)**2 + sum_Szig1(:,2)
              sum_Szig1(:,3) = Szig(:,3,1)**2 + sum_Szig1(:,3)

              sum_Szig2(:,1) = Szig(:,1,2)**2 + sum_Szig2(:,1)
              sum_Szig2(:,2) = Szig(:,2,2)**2 + sum_Szig2(:,2)
              sum_Szig2(:,3) = Szig(:,3,2)**2 + sum_Szig2(:,3)

              sum_Szig3(:,1) = Szig(:,1,3)**2 + sum_Szig3(:,1)
              sum_Szig3(:,2) = Szig(:,2,3)**2 + sum_Szig3(:,2)
              sum_Szig3(:,3) = Szig(:,3,3)**2 + sum_Szig3(:,3)

              sum_Sstripe(:,1) = Sstripe(:,1,1)**2 + Sstripe(:,1,2)**2 + Sstripe(:,1,3)**2 + sum_Sstripe(:,1)
              sum_Sstripe(:,2) = Sstripe(:,2,1)**2 + Sstripe(:,2,2)**2 + Sstripe(:,2,3)**2 + sum_Sstripe(:,2)
              sum_Sstripe(:,3) = Sstripe(:,3,1)**2 + Sstripe(:,3,2)**2 + Sstripe(:,3,3)**2 + sum_Sstripe(:,3)

              sum_Sstripe1(:,1) = Sstripe(:,1,1)**2 + sum_Sstripe1(:,1)
              sum_Sstripe1(:,2) = Sstripe(:,2,1)**2 + sum_Sstripe1(:,2)
              sum_Sstripe1(:,3) = Sstripe(:,3,1)**2 + sum_Sstripe1(:,3)

              sum_Sstripe2(:,1) = Sstripe(:,1,2)**2 + sum_Sstripe2(:,1)
              sum_Sstripe2(:,2) = Sstripe(:,2,2)**2 + sum_Sstripe2(:,2)
              sum_Sstripe2(:,3) = Sstripe(:,3,2)**2 + sum_Sstripe2(:,3)

              sum_Sstripe3(:,1) = Sstripe(:,1,3)**2 + sum_Sstripe3(:,1)
              sum_Sstripe3(:,2) = Sstripe(:,2,3)**2 + sum_Sstripe3(:,2)
              sum_Sstripe3(:,3) = Sstripe(:,3,3)**2 + sum_Sstripe3(:,3)

              !**************************************************

              !**********************zig str all*****************
              do T_i = 1,N_T
                sum_zig(T_i,4) = dot_product(zig(T_i,:,1),zig(T_i,:,1)) + dot_product(zig(T_i,:,2),zig(T_i,:,2))&
                 + dot_product(zig(T_i,:,3),zig(T_i,:,3)) + sum_zig(T_i,4)

                sum_zig1(T_i,4) = dot_product(zig(T_i,:,1),zig(T_i,:,1)) + sum_zig1(T_i,4)
                sum_zig2(T_i,4) = dot_product(zig(T_i,:,2),zig(T_i,:,2)) + sum_zig2(T_i,4)
                sum_zig3(T_i,4) = dot_product(zig(T_i,:,3),zig(T_i,:,3)) + sum_zig3(T_i,4)

                sum_stripe(T_i,4) = dot_product(stripe(T_i,:,1),stripe(T_i,:,1)) + dot_product(stripe(T_i,:,2),stripe(T_i,:,2))&
                 + dot_product(stripe(T_i,:,3),stripe(T_i,:,3)) + sum_stripe(T_i,4)

                sum_stripe1(T_i,4) = dot_product(stripe(T_i,:,1),stripe(T_i,:,1)) + sum_stripe1(T_i,4)
                sum_stripe2(T_i,4) = dot_product(stripe(T_i,:,2),stripe(T_i,:,2)) + sum_stripe2(T_i,4)
                sum_stripe3(T_i,4) = dot_product(stripe(T_i,:,3),stripe(T_i,:,3)) + sum_stripe3(T_i,4)
              end do
              !***************************************************

              !**********************spin space zig str all*****************
              do T_i = 1,N_T
                sum_Szig(T_i,4) = dot_product(Szig(T_i,:,1),Szig(T_i,:,1)) + dot_product(Szig(T_i,:,2),Szig(T_i,:,2))&
                 + dot_product(Szig(T_i,:,3),Szig(T_i,:,3)) + sum_Szig(T_i,4)

                sum_Szig1(T_i,4) = dot_product(Szig(T_i,:,1),Szig(T_i,:,1)) + sum_Szig1(T_i,4)
                sum_Szig2(T_i,4) = dot_product(Szig(T_i,:,2),Szig(T_i,:,2)) + sum_Szig2(T_i,4)
                sum_Szig3(T_i,4) = dot_product(Szig(T_i,:,3),Szig(T_i,:,3)) + sum_Szig3(T_i,4)

                sum_Sstripe(T_i,4) = dot_product(Sstripe(T_i,:,1),Sstripe(T_i,:,1)) + &
                dot_product(Sstripe(T_i,:,2),Sstripe(T_i,:,2))&
                 + dot_product(Sstripe(T_i,:,3),Sstripe(T_i,:,3)) + sum_Sstripe(T_i,4)

                sum_Sstripe1(T_i,4) = dot_product(Sstripe(T_i,:,1),Sstripe(T_i,:,1)) + sum_Sstripe1(T_i,4)
                sum_Sstripe2(T_i,4) = dot_product(Sstripe(T_i,:,2),Sstripe(T_i,:,2)) + sum_Sstripe2(T_i,4)
                sum_Sstripe3(T_i,4) = dot_product(Sstripe(T_i,:,3),Sstripe(T_i,:,3)) + sum_Sstripe3(T_i,4)
              end do
              !***************************************************


              sum_scP(:,1) = sum_scP(:,1) + scP(:,1)/dble(N_spin)!sc-polarization_x
              sum_scP(:,2) = sum_scP(:,2) + scP(:,2)/dble(N_spin)!sc-polarization_y
              sum_scP(:,3) = sum_scP(:,3) + scP(:,3)/dble(N_spin)!sc-polarization_z
              sum_scP(:,4) = sum_scP(:,4) + (scP(:,1)/dble(N_spin))**2!sc-polarization_x^2
              sum_scP(:,5) = sum_scP(:,5) + (scP(:,2)/dble(N_spin))**2!sc-polarization_y^2
              sum_scP(:,6) = sum_scP(:,6) + (scP(:,3)/dble(N_spin))**2!sc-polarization_z^2
              sum_scP2(:,1) = sum_scP2(:,1) + abs(scP(:,1)/dble(N_spin))!sc-polarization_x-abs
              sum_scP2(:,2) = sum_scP2(:,2) + abs(scP(:,2)/dble(N_spin))!sc-polarization_y-abs
              sum_scP2(:,3) = sum_scP2(:,3) + abs(scP(:,3)/dble(N_spin))!sc-polarization_y-abs

              sum_scP_bond(:,1,1) = sum_scP_bond(:,1,1) + (scP_bond(:,1,1)/dble(N_spin))**2!sc-polarization_x
              sum_scP_bond(:,2,1) = sum_scP_bond(:,2,1) + (scP_bond(:,2,1)/dble(N_spin))**2!sc-polarization_y
              sum_scP_bond(:,3,1) = sum_scP_bond(:,3,1) + (scP_bond(:,3,1)/dble(N_spin))**2!sc-polarization_z
              sum_scP_bond(:,1,2) = sum_scP_bond(:,1,2) + (scP_bond(:,1,2)/dble(N_spin))**2!sc-polarization_x^2
              sum_scP_bond(:,2,2) = sum_scP_bond(:,2,2) + (scP_bond(:,2,2)/dble(N_spin))**2!sc-polarization_y^2
              sum_scP_bond(:,3,2) = sum_scP_bond(:,3,2) + (scP_bond(:,3,2)/dble(N_spin))**2!sc-polarization_z^2
              sum_scP_bond(:,1,3) = sum_scP_bond(:,1,3) + (scP_bond(:,1,3)/dble(N_spin))**2!sc-polarization_x^2
              sum_scP_bond(:,2,3) = sum_scP_bond(:,2,3) + (scP_bond(:,2,3)/dble(N_spin))**2!sc-polarization_y^2
              sum_scP_bond(:,3,3) = sum_scP_bond(:,3,3) + (scP_bond(:,3,3)/dble(N_spin))**2!sc-polarization_z^2

              !****intr****
              sum_scP_intr(:,1) = sum_scP_intr(:,1) + scP_intr(:,1)/dble(N_spin)!sc-polarization_x
              sum_scP_intr(:,2) = sum_scP_intr(:,2) + scP_intr(:,2)/dble(N_spin)!sc-polarization_y
              sum_scP_intr(:,3) = sum_scP_intr(:,3) + scP_intr(:,3)/dble(N_spin)!sc-polarization_z
              sum_scP_intr(:,4) = sum_scP_intr(:,4) + (scP_intr(:,1)/dble(N_spin))**2!sc-polarization_x^2
              sum_scP_intr(:,5) = sum_scP_intr(:,5) + (scP_intr(:,2)/dble(N_spin))**2!sc-polarization_y^2
              sum_scP_intr(:,6) = sum_scP_intr(:,6) + (scP_intr(:,3)/dble(N_spin))**2!sc-polarization_z^2
              sum_scP2_intr(:,1) = sum_scP2_intr(:,1) + abs(scP_intr(:,1)/dble(N_spin))!sc-polarization_x-abs
              sum_scP2_intr(:,2) = sum_scP2_intr(:,2) + abs(scP_intr(:,2)/dble(N_spin))!sc-polarization_y-abs
              sum_scP2_intr(:,3) = sum_scP2_intr(:,3) + abs(scP_intr(:,3)/dble(N_spin))!sc-polarization_y-abs
              !************

              !***para*****
              scP_para = scP - scP_intr
              sum_scP_para(:,1) = sum_scP_para(:,1) + scP_para(:,1)/dble(N_spin)!sc-polarization_x
              sum_scP_para(:,2) = sum_scP_para(:,2) + scP_para(:,2)/dble(N_spin)!sc-polarization_y
              sum_scP_para(:,3) = sum_scP_para(:,3) + scP_para(:,3)/dble(N_spin)!sc-polarization_z
              sum_scP_para(:,4) = sum_scP_para(:,4) + (scP_para(:,1)/dble(N_spin))**2!sc-polarization_x^2
              sum_scP_para(:,5) = sum_scP_para(:,5) + (scP_para(:,2)/dble(N_spin))**2!sc-polarization_y^2
              sum_scP_para(:,6) = sum_scP_para(:,6) + (scP_para(:,3)/dble(N_spin))**2!sc-polarization_z^2
              sum_scP2_para(:,1) = sum_scP2_para(:,1) + abs(scP_para(:,1)/dble(N_spin))!sc-polarization_x-abs
              sum_scP2_para(:,2) = sum_scP2_para(:,2) + abs(scP_para(:,2)/dble(N_spin))!sc-polarization_y-abs
              sum_scP2_para(:,3) = sum_scP2_para(:,3) + abs(scP_para(:,3)/dble(N_spin))!sc-polarization_y-abs
              !************

              sum_tmP(:,1) = sum_tmP(:,1) + tmP(:,1)/dble(N_spin*L)!tm-pola_x
              sum_tmP(:,2) = sum_tmP(:,2) + tmP(:,2)/dble(N_spin*L)!tm-pola_y
              sum_tmP(:,3) = sum_tmP(:,3) + tmP(:,3)/dble(N_spin*L)!tm-pola_z
              sum_tmP(:,4) = sum_tmP(:,4) + (tmP(:,1)/dble(N_spin*L))**2!tm-pola_x^2
              sum_tmP(:,5) = sum_tmP(:,5) + (tmP(:,2)/dble(N_spin*L))**2!tm-pola_y^2
              sum_tmP(:,6) = sum_tmP(:,6) + (tmP(:,3)/dble(N_spin*L))**2!tm-pola_z^2
              sum_tmP2(:,1) = sum_tmP(:,1) + abs(tmP(:,1)/dble(N_spin*L))!tm-pola_x-abs
              sum_tmP2(:,2) = sum_tmP(:,2) + abs(tmP(:,2)/dble(N_spin*L))!tm-pola_y-abs
              sum_tmP2(:,3) = sum_tmP(:,3) + abs(tmP(:,3)/dble(N_spin*L))!tm-pola_z-abs

              !*********spin space
              sum_SscP(:,1) = sum_SscP(:,1) + SscP(:,1)/dble(N_spin)!Ssc-polarization_x
              sum_SscP(:,2) = sum_SscP(:,2) + SscP(:,2)/dble(N_spin)!Ssc-polarization_y
              sum_SscP(:,3) = sum_SscP(:,3) + SscP(:,3)/dble(N_spin)!Ssc-polarization_z
              sum_SscP(:,4) = sum_SscP(:,4) + (SscP(:,1)/dble(N_spin))**2!Ssc-polarization_x^2
              sum_SscP(:,5) = sum_SscP(:,5) + (SscP(:,2)/dble(N_spin))**2!Ssc-polarization_y^2
              sum_SscP(:,6) = sum_SscP(:,6) + (SscP(:,3)/dble(N_spin))**2!Ssc-polarization_z^2
              sum_SscP2(:,1) = sum_SscP2(:,1) + abs(SscP(:,1)/dble(N_spin))!Ssc-polarization_x-abs
              sum_SscP2(:,2) = sum_SscP2(:,2) + abs(SscP(:,2)/dble(N_spin))!Ssc-polarization_y-abs
              sum_SscP2(:,3) = sum_SscP2(:,3) + abs(SscP(:,3)/dble(N_spin))!Ssc-polarization_y-abs

              sum_StmP(:,1) = sum_StmP(:,1) + StmP(:,1)/dble(N_spin*L)!Stm-pola_x
              sum_StmP(:,2) = sum_StmP(:,2) + StmP(:,2)/dble(N_spin*L)!Stm-pola_y
              sum_StmP(:,3) = sum_StmP(:,3) + StmP(:,3)/dble(N_spin*L)!Stm-pola_z
              sum_StmP(:,4) = sum_StmP(:,4) + (StmP(:,1)/dble(N_spin*L))**2!Stm-pola_x^2
              sum_StmP(:,5) = sum_StmP(:,5) + (StmP(:,2)/dble(N_spin*L))**2!Stm-pola_y^2
              sum_StmP(:,6) = sum_StmP(:,6) + (StmP(:,3)/dble(N_spin*L))**2!Stm-pola_z^2
              sum_StmP2(:,1) = sum_StmP(:,1) + abs(StmP(:,1)/dble(N_spin*L))!Stm-pola_x-abs
              sum_StmP2(:,2) = sum_StmP(:,2) + abs(StmP(:,2)/dble(N_spin*L))!Stm-pola_y-abs
              sum_StmP2(:,3) = sum_StmP(:,3) + abs(StmP(:,3)/dble(N_spin*L))!tm-pola_z-abs
              !***********

              !****intr****
              sum_SscP_intr(:,1) = sum_SscP_intr(:,1) + SscP_intr(:,1)/dble(N_spin)!sc-polarization_x
              sum_SscP_intr(:,2) = sum_SscP_intr(:,2) + SscP_intr(:,2)/dble(N_spin)!sc-polarization_y
              sum_SscP_intr(:,3) = sum_SscP_intr(:,3) + SscP_intr(:,3)/dble(N_spin)!sc-polarization_z
              sum_SscP_intr(:,4) = sum_SscP_intr(:,4) + (SscP_intr(:,1)/dble(N_spin))**2!sc-polarization_x^2
              sum_SscP_intr(:,5) = sum_SscP_intr(:,5) + (SscP_intr(:,2)/dble(N_spin))**2!sc-polarization_y^2
              sum_SscP_intr(:,6) = sum_SscP_intr(:,6) + (SscP_intr(:,3)/dble(N_spin))**2!sc-polarization_z^2
              sum_SscP2_intr(:,1) = sum_SscP2_intr(:,1) + abs(SscP_intr(:,1)/dble(N_spin))!sc-polarization_x-abs
              sum_SscP2_intr(:,2) = sum_SscP2_intr(:,2) + abs(SscP_intr(:,2)/dble(N_spin))!sc-polarization_y-abs
              sum_SscP2_intr(:,3) = sum_SscP2_intr(:,3) + abs(SscP_intr(:,3)/dble(N_spin))!sc-polarization_y-abs
              !************
              !***para*****
              SscP_para = SscP - SscP_intr
              sum_SscP_para(:,1) = sum_SscP_para(:,1) + SscP_para(:,1)/dble(N_spin)!sc-polarization_x
              sum_SscP_para(:,2) = sum_SscP_para(:,2) + SscP_para(:,2)/dble(N_spin)!sc-polarization_y
              sum_SscP_para(:,3) = sum_SscP_para(:,3) + SscP_para(:,3)/dble(N_spin)!sc-polarization_z
              sum_SscP_para(:,4) = sum_SscP_para(:,4) + (SscP_para(:,1)/dble(N_spin))**2!sc-polarization_x^2
              sum_SscP_para(:,5) = sum_SscP_para(:,5) + (SscP_para(:,2)/dble(N_spin))**2!sc-polarization_y^2
              sum_SscP_para(:,6) = sum_SscP_para(:,6) + (SscP_para(:,3)/dble(N_spin))**2!sc-polarization_z^2
              sum_SscP2_para(:,1) = sum_SscP2_para(:,1) + abs(SscP_para(:,1)/dble(N_spin))!sc-polarization_x-abs
              sum_SscP2_para(:,2) = sum_SscP2_para(:,2) + abs(SscP_para(:,2)/dble(N_spin))!sc-polarization_y-abs
              sum_SscP2_para(:,3) = sum_SscP2_para(:,3) + abs(SscP_para(:,3)/dble(N_spin))!sc-polarization_y-abs
              !************


              sum_mono = sum_mono + (mono/dble(N_spin*L))**2
              sum_quad(:,1) = sum_quad(:,1) + (quad(:,1)/dble(N_spin*L))**2
              sum_quad(:,2) = sum_quad(:,2) + (quad(:,2)/dble(N_spin*L))**2
              sum_quad(:,3) = sum_quad(:,3) + (quad(:,3)/dble(N_spin*L))**2
              sum_quad(:,4) = sum_quad(:,4) + (quad(:,4)/dble(N_spin*L))**2
              sum_quad(:,5) = sum_quad(:,5) + (quad(:,5)/dble(N_spin*L))**2
              sum_quad(:,6) = sum_quad(:,6) + (quad(:,6)/dble(N_spin*L))**2
              sum_quad(:,7) = sum_quad(:,7) + (quad(:,7)/dble(N_spin*L))**2
              sum_quad(:,8) = sum_quad(:,8) + (quad(:,8)/dble(N_spin*L))**2
              sum_quad(:,9) = sum_quad(:,9) + (quad(:,9)/dble(N_spin*L))**2

              sum_vmAF(:,:)=sum_vmAF(:,:)+AF/dble(N_spin)
              sum_vm = sum_vm +m/dble(N_spin)

              sum_k_norm(:,1) = sum_k_norm(:,1) + k(:,1,1)/dble(N_spin)*k(:,1,1)/dble(N_spin) + &
                   k(:,1,2)/dble(N_spin)*k(:,1,2)/dble(N_spin) + k(:,1,3)/dble(N_spin)*k(:,1,3)/dble(N_spin)
              sum_k_norm(:,2) = sum_k_norm(:,2) + k(:,2,1)/dble(N_spin)*k(:,2,1)/dble(N_spin) + &
                   k(:,2,2)/dble(N_spin)*k(:,2,2)/dble(N_spin) + k(:,2,3)/dble(N_spin)*k(:,2,3)/dble(N_spin)
              sum_k_norm(:,3) = sum_k_norm(:,3) + k(:,3,1)/dble(N_spin)*k(:,3,1)/dble(N_spin) + &
                   k(:,3,2)/dble(N_spin)*k(:,3,2)/dble(N_spin) + k(:,3,3)/dble(N_spin)*k(:,3,3)/dble(N_spin)
              sum_k_norm(:,4) = sum_k_norm(:,4) + k(:,4,1)/dble(N_spin)*k(:,4,1)/dble(N_spin) + &
                        k(:,4,2)/dble(N_spin)*k(:,4,2)/dble(N_spin) + k(:,4,3)/dble(N_spin)*k(:,4,3)/dble(N_spin)


           end if

           call stock_random(randnum)
           count = 1


           do T_i = 1, N_T
!            call HeatBath(rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
!            call over_relaxation(rg,spinA,spinB,T_i,J_ij)
             !$omp parallel do num_threads(N_core)
             do ii = 1,L
               do jj = 1,L
                 do kk = 1,L
                   !count = 2*(kk-1)*L*L+2*(jj-1)*L+2*ii-1
                   call HeatBathA(ii,jj,kk,rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
                 end do
               end do
             end do
             !$omp end parallel do

             !$omp parallel do num_threads(N_core)
             do ii = 1,L
               do jj = 1,L
                 do kk = 1,L
                   !count = 2*(kk-1)*L*L+2*(jj-1)*L+2*ii-1+2*L*L*L
                   call HeatBathB(ii,jj,kk,rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
                 end do
               end do
             end do
             !$omp end parallel do

             do R_i = 1,N_relax
               !$omp parallel do num_threads(N_core)
               do ii = 1,L
                 do jj = 1,L
                   do kk = 1,L
                     call over_relaxationA(ii,jj,kk,rg,spinA,spinB,T_i,J_ij)
                   end do
                 end do
               end do
               !$omp end parallel do
!             end do

!             do R_i = 1,N_relax
               !$omp parallel do num_threads(N_core)
               do ii = 1,L
                 do jj = 1,L
                   do kk = 1,L
                     call over_relaxationB(ii,jj,kk,rg,spinA,spinB,T_i,J_ij)
                   end do
                 end do
               end do
               !$omp end parallel do
             end do

           end do
        end do
        !************************end MC sweep********************************
        ! write(*,*)counter


        !********************thermal average*********************
        do T_i = 1,N_T

           tmp_i = T(N_T)

           rmsm(T_i) = sqrt(sum_rmsm(T_i)/dble(sweep2))
           C(T_i) = (dble(N_spin)/tmp_i**2)*(sum_e(T_i,2)/dble(sweep2) - (sum_e(T_i,1)/dble(sweep2))**2)
           E(T_i) = sum_e(T_i,1)/dble(sweep2) !energy-density
           sum_k_norm = sum_k_norm/dble(sweep2)
           AF(T_i,1) = sqrt(sum_AF(T_i,4)/dble(sweep2))!staggered_x
           AF(T_i,2) = sqrt(sum_AF(T_i,5)/dble(sweep2))!staggered_y
           AF(T_i,3) = sqrt(sum_AF(T_i,6)/dble(sweep2))!staggered_z
           xAF(T_i,1) = dble(N_spin)/tmp_i * (sum_AF(T_i,4)/dble(sweep2) - (sum_AF(T_i,1)/dble(sweep2))**2)!AFsuSsceptibity_x
           xAF(T_i,2) = dble(N_spin)/tmp_i * (sum_AF(T_i,5)/dble(sweep2) - (sum_AF(T_i,2)/dble(sweep2))**2)!AFsuSsceptibity_y
           xAF(T_i,3) = dble(N_spin)/tmp_i * (sum_AF(T_i,6)/dble(sweep2) - (sum_AF(T_i,3)/dble(sweep2))**2)!AFsusceptibity_z
           tm(T_i,1) = sqrt(sum_tm(T_i,1)/dble(sweep2))!troidalmoment-bond1
           tm(T_i,2) = sqrt(sum_tm(T_i,2)/dble(sweep2))!troidalmoment-bond2
           tm(T_i,3) = sqrt(sum_tm(T_i,3)/dble(sweep2))!troidalmoment-bond3
           scP(T_i,1) = sum_scP(T_i,1)/dble(sweep2)!spincurrent-polarization_x
           scP(T_i,2) = sum_scP(T_i,2)/dble(sweep2)!spincurrent-polarization_y
           absscP(T_i,1) = sqrt(sum_scP(T_i,4)/dble(sweep2))!spincurrent-polarization_x-abs
           absscP(T_i,2) = sqrt(sum_scP(T_i,5)/dble(sweep2))!spincurrent-polarization_y-abs
           tmP(T_i,1) = sum_tmP(T_i,1)/dble(sweep2)!troidalmoment-polarization_x
           tmP(T_i,2) = sum_tmP(T_i,2)/dble(sweep2)!troidalmoment-polarization_y
           tmP(T_i,3) = sum_tmP(T_i,3)/dble(sweep2)!troidalmoment-polarization_z
           abstmP(T_i,1) = sqrt(sum_tmP(T_i,4)/dble(sweep2))!troidalmoment-polarization_x-abs
           abstmP(T_i,2) = sqrt(sum_tmP(T_i,5)/dble(sweep2))!troidalmoment-polarization_y-abs
           abstmP(T_i,3) = sqrt(sum_tmP(T_i,6)/dble(sweep2))!troidalmoment-polarization_z-abs
           scperm(T_i,1) = dble(N_spin)/tmp_i * (absscP(T_i,1)**2 - (sum_scP2(T_i,1)/dble(sweep2))**2)!scpermittivity_x
           scperm(T_i,2) = dble(N_spin)/tmp_i * (absscP(T_i,2)**2 - (sum_scP2(T_i,2)/dble(sweep2))**2)!scpermittivity_y
           tmperm(T_i,1) = dble(N_spin*L)/tmp_i * ( abstmP(T_i,1)**2 - (sum_tmP2(T_i,1)/dble(sweep2))**2)!tmpermittivity_x
           tmperm(T_i,2) = dble(N_spin*L)/tmp_i * ( abstmP(T_i,2)**2 - (sum_tmP2(T_i,2)/dble(sweep2))**2)!tmpermittivity_y
           tmperm(T_i,3) = dble(N_spin*L)/tmp_i * ( abstmP(T_i,3)**2 - (sum_tmP2(T_i,3)/dble(sweep2))**2)!tmpermittivity_z
           X_SG_kmin(T_i,2) = sum_q_uv_kmin(T_i,1)/dble(sweep2) !coreration ratio guzai_L_xx
           X_SG_kmin(T_i,3) = sum_q_uv_kmin(T_i,4)/dble(sweep2) !coreration ratio guzai_L_yy
           X_SG(T_i,1) = X_SG(T_i,1)/dble(sweep2)!spinoverlap-tensor__xx
           X_SG(T_i,2) = X_SG(T_i,2)/dble(sweep2)!spinoverlap-tensor__xy
           X_SG(T_i,3) = X_SG(T_i,3)/dble(sweep2)!spinoverlap-tensor__yx
           X_SG(T_i,4) = X_SG(T_i,4)/dble(sweep2)!spinoverlap-tensor__yy
           X_SG(T_i,5) = sumX_SG(T_i,1)/dble(sweep2)!sum SG_X
           X_SG_kmin(T_i,1) = sumX_SG(T_i,2)/dble(sweep2)
           X_CG(T_i,1) = sum_q_k(T_i,1)/dble(sweep2)
           X_CG(T_i,2) = sum_q_k(T_i,2)/dble(sweep2)
           X_CG(T_i,3) = sum_q_k_kmin(T_i,1)/dble(sweep2)
           X_CG(T_i,4) = sum_q_k_kmin(T_i,2)/dble(sweep2)

           !**************************write physical quantities***************************
           write(seed,*)tmp_i,E(T_i),C(T_i),rmsm(T_i),sum_k_norm(T_i,1:3),scP(T_i,1),scP(T_i,2),&
                absscP(T_i,1),absscP(T_i,2),xAF(T_i,1),xAF(T_i,2),AF(T_i,1),AF(T_i,2),&
                tm(T_i,1),tm(T_i,2),tm(T_i,3),tmP(T_i,1),tmP(T_i,2),tmP(T_i,3),&
                abstmP(T_i,1),abstmP(T_i,2),abstmP(T_i,3),scperm(T_i,1),scperm(T_i,2),&
                tmperm(T_i,1),tmperm(T_i,2),tmperm(T_i,3),X_SG(T_i,1),X_SG(T_i,2),X_SG(T_i,3),X_SG(T_i,4),&
                X_SG(T_i,5), X_SG_kmin(T_i,1), X_SG_kmin(T_i,2), X_SG_kmin(T_i,3),&
                X_CG(T_i,1), X_CG(T_i,2), X_CG(T_i,3), X_CG(T_i,4)
           !******************************************************************************

           mono(T_i) = sqrt(sum_mono(T_i)/dble(sweep2))
           quad(T_i,1) = sqrt(sum_quad(T_i,1)/dble(sweep2))
           quad(T_i,2) = sqrt(sum_quad(T_i,2)/dble(sweep2))
           quad(T_i,3) = sqrt(sum_quad(T_i,3)/dble(sweep2))
           quad(T_i,4) = sqrt(sum_quad(T_i,4)/dble(sweep2))
           quad(T_i,5) = sqrt(sum_quad(T_i,5)/dble(sweep2))
           quad(T_i,6) = sqrt(sum_quad(T_i,6)/dble(sweep2))
           quad(T_i,7) = sqrt(sum_quad(T_i,7)/dble(sweep2))
           quad(T_i,8)  = sqrt(sum_quad(T_i,8)/dble(sweep2))
           quad(T_i,9) = sqrt(sum_quad(T_i,9)/dble(sweep2))
           AF(T_i,1) = (sum_vmAF(T_i,1)/dble(sweep2))!staggered_x
           AF(T_i,2) = (sum_vmAF(T_i,2)/dble(sweep2))!staggered_y
           m(T_i,1) = (sum_vm(T_i,1)/dble(sweep2))!staggered_x
           m(T_i,2) = (sum_vm(T_i,2)/dble(sweep2))!staggered_y
           sus = dble(N_spin)/tmp_i*(rmsm(T_i)**2 - (sum_rmsm2(T_i)/dble(sweep2))**2) !Susceptibility
           k(T_i,1,:) = sqrt(sum_k(T_i,1,:)/dble(sweep2))!chiral-bond1
           k(T_i,2,:) = sqrt(sum_k(T_i,2,:)/dble(sweep2))!chiral-bond2
           k(T_i,3,:) = sqrt(sum_k(T_i,3,:)/dble(sweep2))!chiral-bond3
           scP(T_i,3) = sum_scP(T_i,3)/dble(sweep2)!spincurrent-polarization_y
           absscP(T_i,3) = sqrt(sum_scP(T_i,6)/dble(sweep2))!spincurrent-polarization_x-abs
           scperm(T_i,3) = dble(N_spin)/tmp_i * (absscP(T_i,3)**2 - (sum_scP2(T_i,3)/dble(sweep2))**2)!scpermittivity_y
           rmsm(T_i) = sqrt(sum_rmsAF(T_i,1)/dble(sweep2))
           X_SG(T_i,1) = dble(N_spin)/tmp_i*(rmsm(T_i)**2 - (sum_rmsAF(T_i,2)/dble(sweep2))**2) !AFSusceptibility

           !**************************write physical quantities***************************
           write(fo2,*)tmp_i,mono(T_i),quad(T_i,1),quad(T_i,2),quad(T_i,3),quad(T_i,4),&
                quad(T_i,5),quad(T_i,6),quad(T_i,7),quad(T_i,8),quad(T_i,9),&
                AF(T_i,1),AF(T_i,2),m(T_i,1),m(T_i,2),sus,&
                k(T_i,1,:),k(T_i,2,:),k(T_i,3,:),scP(T_i,3),absscP(T_i,3),scperm(T_i,3),AF(T_i,3),xAF(T_i,3),rmsm(T_i),X_SG(T_i,1)
           !******************************************************************************


           sum_zig(T_i,1) = sqrt(sum_zig(T_i,1)/dble(sweep2))
           sum_zig(T_i,2) = sqrt(sum_zig(T_i,2)/dble(sweep2))
           sum_zig(T_i,3) = sqrt(sum_zig(T_i,3)/dble(sweep2))
           sum_zig(T_i,4) = sqrt(sum_zig(T_i,4)/dble(sweep2))
           sum_zig1(T_i,1) = sqrt(sum_zig1(T_i,1)/dble(sweep2))
           sum_zig1(T_i,2) = sqrt(sum_zig1(T_i,2)/dble(sweep2))
           sum_zig1(T_i,3) = sqrt(sum_zig1(T_i,3)/dble(sweep2))
           sum_zig1(T_i,4) = sqrt(sum_zig1(T_i,4)/dble(sweep2))
           sum_zig2(T_i,1) = sqrt(sum_zig2(T_i,1)/dble(sweep2))
           sum_zig2(T_i,2) = sqrt(sum_zig2(T_i,2)/dble(sweep2))
           sum_zig2(T_i,3) = sqrt(sum_zig2(T_i,3)/dble(sweep2))
           sum_zig2(T_i,4) = sqrt(sum_zig2(T_i,4)/dble(sweep2))
           sum_zig3(T_i,1) = sqrt(sum_zig3(T_i,1)/dble(sweep2))
           sum_zig3(T_i,2) = sqrt(sum_zig3(T_i,2)/dble(sweep2))
           sum_zig3(T_i,3) = sqrt(sum_zig3(T_i,3)/dble(sweep2))
           sum_zig3(T_i,4) = sqrt(sum_zig3(T_i,4)/dble(sweep2))

           sum_stripe(T_i,1) = sqrt(sum_stripe(T_i,1)/dble(sweep2))
           sum_stripe(T_i,2) = sqrt(sum_stripe(T_i,2)/dble(sweep2))
           sum_stripe(T_i,3) = sqrt(sum_stripe(T_i,3)/dble(sweep2))
           sum_stripe(T_i,4) = sqrt(sum_stripe(T_i,4)/dble(sweep2))
           sum_stripe1(T_i,1) = sqrt(sum_stripe1(T_i,1)/dble(sweep2))
           sum_stripe1(T_i,2) = sqrt(sum_stripe1(T_i,2)/dble(sweep2))
           sum_stripe1(T_i,3) = sqrt(sum_stripe1(T_i,3)/dble(sweep2))
           sum_stripe1(T_i,4) = sqrt(sum_stripe1(T_i,4)/dble(sweep2))
           sum_stripe2(T_i,1) = sqrt(sum_stripe2(T_i,1)/dble(sweep2))
           sum_stripe2(T_i,2) = sqrt(sum_stripe2(T_i,2)/dble(sweep2))
           sum_stripe2(T_i,3) = sqrt(sum_stripe2(T_i,3)/dble(sweep2))
           sum_stripe2(T_i,4) = sqrt(sum_stripe2(T_i,4)/dble(sweep2))
           sum_stripe3(T_i,1) = sqrt(sum_stripe3(T_i,1)/dble(sweep2))
           sum_stripe3(T_i,2) = sqrt(sum_stripe3(T_i,2)/dble(sweep2))
           sum_stripe3(T_i,3) = sqrt(sum_stripe3(T_i,3)/dble(sweep2))
           sum_stripe3(T_i,4) = sqrt(sum_stripe3(T_i,4)/dble(sweep2))
           sum_pm(T_i,1) =  sqrt(sum_pm(T_i,1)/dble(sweep2))
           sum_pm(T_i,2) =  sqrt(sum_pm(T_i,2) /dble(sweep2))
           sum_pm(T_i,3) =  sqrt(sum_pm(T_i,3)/dble(sweep2))
           sum_pm(T_i,4) =  sqrt(sum_pm(T_i,4)/dble(sweep2))
           sum_m(T_i,1) =  sqrt(sum_m(T_i,1) /dble(sweep2))
           sum_m(T_i,2) =  sqrt(sum_m(T_i,2)/dble(sweep2))
           sum_m(T_i,3) =  sqrt(sum_m(T_i,3)/dble(sweep2))
           k(T_i,4,:) = sqrt(sum_k(T_i,4,:)/dble(sweep2))!chiral-inter

           !**scPintrpara****
           absscP_intr(T_i,1) = sqrt(sum_scP_intr(T_i,4)/dble(sweep2))!spincurrent-polarization_x-abs
           absscP_intr(T_i,2) = sqrt(sum_scP_intr(T_i,5)/dble(sweep2))!spincurrent-polarization_y-abs
           absscP_intr(T_i,3) = sqrt(sum_scP_intr(T_i,6)/dble(sweep2))!spincurrent-polarization_z-abs
           scperm_intr(T_i,1) = dble(N_spin)/tmp_i * (absscP_intr(T_i,1)**2 - (sum_scP2_intr(T_i,1)/dble(sweep2))**2)!scpermittivity_x
           scperm_intr(T_i,2) = dble(N_spin)/tmp_i * (absscP_intr(T_i,2)**2 - (sum_scP2_intr(T_i,2)/dble(sweep2))**2)!scpermittivity_y
           scperm_intr(T_i,3) = dble(N_spin)/tmp_i * (absscP_intr(T_i,3)**2 - (sum_scP2_intr(T_i,3)/dble(sweep2))**2)!scpermittivity_z

           absscP_para(T_i,1) = sqrt(sum_scP_para(T_i,4)/dble(sweep2))!spincurrent-polarization_x-abs
           absscP_para(T_i,2) = sqrt(sum_scP_para(T_i,5)/dble(sweep2))!spincurrent-polarization_y-abs
           absscP_para(T_i,3) = sqrt(sum_scP_para(T_i,6)/dble(sweep2))!spincurrent-polarization_z-abs
           scperm_para(T_i,1) = dble(N_spin)/tmp_i * (absscP_para(T_i,1)**2 - (sum_scP2_para(T_i,1)/dble(sweep2))**2)!scpermittivity_x
           scperm_para(T_i,2) = dble(N_spin)/tmp_i * (absscP_para(T_i,2)**2 - (sum_scP2_para(T_i,2)/dble(sweep2))**2)!scpermittivity_y
           scperm_para(T_i,3) = dble(N_spin)/tmp_i * (absscP_para(T_i,3)**2 - (sum_scP2_para(T_i,3)/dble(sweep2))**2)!scpermittivity_z
           !*****************

           !**
           sum_scP_bond(T_i,:,:) = sqrt(sum_scP_bond(T_i,:,:)/dble(sweep2))
           !**

           !**************************write physical quantities***************************
           write(fo3,*)tmp_i,sum_zig(T_i,1),sum_zig(T_i,2),sum_zig(T_i,3),sum_zig(T_i,4)&
           ,sum_zig1(T_i,1),sum_zig1(T_i,2),sum_zig1(T_i,3),sum_zig1(T_i,4)&
           ,sum_zig2(T_i,1),sum_zig2(T_i,2),sum_zig2(T_i,3),sum_zig2(T_i,4)&
           ,sum_zig3(T_i,1),sum_zig3(T_i,2),sum_zig3(T_i,3),sum_zig3(T_i,4)&
           ,sum_stripe(T_i,1),sum_stripe(T_i,2),sum_stripe(T_i,3),sum_stripe(T_i,4)&
           ,sum_stripe1(T_i,1),sum_stripe1(T_i,2),sum_stripe1(T_i,3),sum_stripe1(T_i,4)&
           ,sum_stripe2(T_i,1),sum_stripe2(T_i,2),sum_stripe2(T_i,3),sum_stripe2(T_i,4)&
           ,sum_stripe3(T_i,1),sum_stripe3(T_i,2),sum_stripe3(T_i,3),sum_stripe3(T_i,4)&
           ,sum_pm(T_i,1),sum_pm(T_i,2),sum_pm(T_i,3),sum_pm(T_i,4)&
           ,sum_m(T_i,1),sum_m(T_i,2),sum_m(T_i,3),k(T_i,4,:),sum_k_norm(T_i,4)&
           ,absscP_intr(T_i,:),absscP_para(T_i,:),scperm_intr(T_i,:),scperm_para(T_i,:)&
           ,sum_scP_bond(T_i,1,1),sum_scP_bond(T_i,2,1),sum_scP_bond(T_i,3,1)&
           ,sum_scP_bond(T_i,1,2),sum_scP_bond(T_i,2,2),sum_scP_bond(T_i,3,2)&
           ,sum_scP_bond(T_i,1,3),sum_scP_bond(T_i,2,3),sum_scP_bond(T_i,3,3)
           !******************************************************************************


           sum_Szig(T_i,1) = sqrt(sum_Szig(T_i,1)/dble(sweep2))
           sum_Szig(T_i,2) = sqrt(sum_Szig(T_i,2)/dble(sweep2))
           sum_Szig(T_i,3) = sqrt(sum_Szig(T_i,3)/dble(sweep2))
           sum_Szig1(T_i,1) = sqrt(sum_Szig1(T_i,1)/dble(sweep2))
           sum_Szig1(T_i,2) = sqrt(sum_Szig1(T_i,2)/dble(sweep2))
           sum_Szig1(T_i,3) = sqrt(sum_Szig1(T_i,3)/dble(sweep2))
           sum_Szig2(T_i,1) = sqrt(sum_Szig2(T_i,1)/dble(sweep2))
           sum_Szig2(T_i,2) = sqrt(sum_Szig2(T_i,2)/dble(sweep2))
           sum_Szig2(T_i,3) = sqrt(sum_Szig2(T_i,3)/dble(sweep2))
           sum_Szig3(T_i,1) = sqrt(sum_Szig3(T_i,1)/dble(sweep2))
           sum_Szig3(T_i,2) = sqrt(sum_Szig3(T_i,2)/dble(sweep2))
           sum_Szig3(T_i,3) = sqrt(sum_Szig3(T_i,3)/dble(sweep2))

           sum_Sstripe(T_i,1) = sqrt(sum_Sstripe(T_i,1)/dble(sweep2))
           sum_Sstripe(T_i,2) = sqrt(sum_Sstripe(T_i,2)/dble(sweep2))
           sum_Sstripe(T_i,3) = sqrt(sum_Sstripe(T_i,3)/dble(sweep2))
           sum_Sstripe1(T_i,1) = sqrt(sum_Sstripe1(T_i,1)/dble(sweep2))
           sum_Sstripe1(T_i,2) = sqrt(sum_Sstripe1(T_i,2)/dble(sweep2))
           sum_Sstripe1(T_i,3) = sqrt(sum_Sstripe1(T_i,3)/dble(sweep2))
           sum_Sstripe2(T_i,1) = sqrt(sum_Sstripe2(T_i,1)/dble(sweep2))
           sum_Sstripe2(T_i,2) = sqrt(sum_Sstripe2(T_i,2)/dble(sweep2))
           sum_Sstripe2(T_i,3) = sqrt(sum_Sstripe2(T_i,3)/dble(sweep2))
           sum_Sstripe3(T_i,1) = sqrt(sum_Sstripe3(T_i,1)/dble(sweep2))
           sum_Sstripe3(T_i,2) = sqrt(sum_Sstripe3(T_i,2)/dble(sweep2))
           sum_Sstripe3(T_i,3) = sqrt(sum_Sstripe3(T_i,3)/dble(sweep2))
           sum_Spm(T_i,1) =  sqrt(sum_Spm(T_i,1)/dble(sweep2))
           sum_Spm(T_i,2) =  sqrt(sum_Spm(T_i,2) /dble(sweep2))
           sum_Spm(T_i,3) =  sqrt(sum_Spm(T_i,3)/dble(sweep2))
           sum_Sm(T_i,1) =  sqrt(sum_Sm(T_i,1) /dble(sweep2))
           sum_Sm(T_i,2) =  sqrt(sum_Sm(T_i,2)/dble(sweep2))
           sum_Sm(T_i,3) =  sqrt(sum_Sm(T_i,3)/dble(sweep2))

           Sk(T_i,1,:) = sum_Sk(T_i,1,:)/dble(sweep2)!chiral-bond1
           Sk(T_i,2,:) = sum_Sk(T_i,2,:)/dble(sweep2)!chiral-bond2
           Sk(T_i,3,:) = sum_Sk(T_i,3,:)/dble(sweep2)!chiral-bond3
           Sk(T_i,4,:) = sum_Sk(T_i,4,:)/dble(sweep2)!chiral-inter

           SAF(T_i,1) = sqrt(sum_SAF(T_i,4)/dble(sweep2))!staggered_x
           SAF(T_i,2) = sqrt(sum_SAF(T_i,5)/dble(sweep2))!staggered_y
           SAF(T_i,3) = sqrt(sum_SAF(T_i,6)/dble(sweep2))!staggered_z
           SxAF(T_i,1) = dble(N_spin)/tmp_i * (sum_SAF(T_i,4)/dble(sweep2) - (sum_SAF(T_i,1)/dble(sweep2))**2)!AFsuSsceptibity_x
           SxAF(T_i,2) = dble(N_spin)/tmp_i * (sum_SAF(T_i,5)/dble(sweep2) - (sum_SAF(T_i,2)/dble(sweep2))**2)!AFsuSsceptibity_y
           SxAF(T_i,3) = dble(N_spin)/tmp_i * (sum_SAF(T_i,6)/dble(sweep2) - (sum_SAF(T_i,3)/dble(sweep2))**2)!AFsusceptibity_z
           Stm(T_i,1) = sqrt(sum_Stm(T_i,1)/dble(sweep2))!troidalmoment-bond1
           Stm(T_i,2) = sqrt(sum_Stm(T_i,2)/dble(sweep2))!troidalmoment-bond2
           Stm(T_i,3) = sqrt(sum_Stm(T_i,3)/dble(sweep2))!troidalmoment-bond3
           SscP(T_i,1) = sum_SscP(T_i,1)/dble(sweep2)!spincurrent-polarization_x
           SscP(T_i,2) = sum_SscP(T_i,2)/dble(sweep2)!spincurrent-polarization_y
           SscP(T_i,3) = sum_SscP(T_i,3)/dble(sweep2)!spincurrent-polarization_y
           SabsscP(T_i,1) = sqrt(sum_SscP(T_i,4)/dble(sweep2))!spincurrent-polarization_x-abs
           SabsscP(T_i,2) = sqrt(sum_SscP(T_i,5)/dble(sweep2))!spincurrent-polarization_y-abs
           SabsscP(T_i,3) = sqrt(sum_SscP(T_i,6)/dble(sweep2))!spincurrent-polarization_y-abs
           StmP(T_i,1) = sum_StmP(T_i,1)/dble(sweep2)!troidalmoment-polarization_x
           StmP(T_i,2) = sum_StmP(T_i,2)/dble(sweep2)!troidalmoment-polarization_y
           StmP(T_i,3) = sum_StmP(T_i,3)/dble(sweep2)!troidalmoment-polarization_z
           SabstmP(T_i,1) = sqrt(sum_StmP(T_i,4)/dble(sweep2))!troidalmoment-polarization_x-abs
           SabstmP(T_i,2) = sqrt(sum_StmP(T_i,5)/dble(sweep2))!troidalmoment-polarization_y-abs
           SabstmP(T_i,3) = sqrt(sum_StmP(T_i,6)/dble(sweep2))!troidalmoment-polarization_z-abs
           Sscperm(T_i,1) = dble(N_spin)/tmp_i * (SabsscP(T_i,1)**2 - (sum_SscP2(T_i,1)/dble(sweep2))**2)!scpermittivity_x
           Sscperm(T_i,2) = dble(N_spin)/tmp_i * (SabsscP(T_i,2)**2 - (sum_SscP2(T_i,2)/dble(sweep2))**2)!scpermittivity_y
           Sscperm(T_i,3) = dble(N_spin)/tmp_i * (SabsscP(T_i,3)**2 - (sum_SscP2(T_i,3)/dble(sweep2))**2)!scpermittivity_y
           Stmperm(T_i,1) = dble(N_spin*L)/tmp_i * ( SabstmP(T_i,1)**2 - (sum_StmP2(T_i,1)/dble(sweep2))**2)!tmpermittivity_x
           Stmperm(T_i,2) = dble(N_spin*L)/tmp_i * ( SabstmP(T_i,2)**2 - (sum_StmP2(T_i,2)/dble(sweep2))**2)!tmpermittivity_y
           Stmperm(T_i,3) = dble(N_spin*L)/tmp_i * ( SabstmP(T_i,3)**2 - (sum_StmP2(T_i,3)/dble(sweep2))**2)!tmpermittivity_z

           !****scpparaintr***
           SabsscP_intr(T_i,1) = sqrt(sum_SscP_intr(T_i,4)/dble(sweep2))!spincurrent-polarization_x-abs
           SabsscP_intr(T_i,2) = sqrt(sum_SscP_intr(T_i,5)/dble(sweep2))!spincurrent-polarization_y-abs
           SabsscP_intr(T_i,3) = sqrt(sum_SscP_intr(T_i,6)/dble(sweep2))!spincurrent-polarization_y-abs
           Sscperm_intr(T_i,1) = dble(N_spin)/tmp_i * (SabsscP_intr(T_i,1)**2 - (sum_SscP2_intr(T_i,1)/dble(sweep2))**2)!scpermittivity_x
           Sscperm_intr(T_i,2) = dble(N_spin)/tmp_i * (SabsscP_intr(T_i,2)**2 - (sum_SscP2_intr(T_i,2)/dble(sweep2))**2)!scpermittivity_y
           Sscperm_intr(T_i,3) = dble(N_spin)/tmp_i * (SabsscP_intr(T_i,3)**2 - (sum_SscP2_intr(T_i,3)/dble(sweep2))**2)!scpermittivity_y

           SabsscP_para(T_i,1) = sqrt(sum_SscP_para(T_i,4)/dble(sweep2))!spincurrent-polarization_x-abs
           SabsscP_para(T_i,2) = sqrt(sum_SscP_para(T_i,5)/dble(sweep2))!spincurrent-polarization_y-abs
           SabsscP_para(T_i,3) = sqrt(sum_SscP_para(T_i,6)/dble(sweep2))!spincurrent-polarization_y-abs
           Sscperm_para(T_i,1) = dble(N_spin)/tmp_i * (SabsscP_para(T_i,1)**2 - (sum_SscP2_para(T_i,1)/dble(sweep2))**2)!scpermittivity_x
           Sscperm_para(T_i,2) = dble(N_spin)/tmp_i * (SabsscP_para(T_i,2)**2 - (sum_SscP2_para(T_i,2)/dble(sweep2))**2)!scpermittivity_y
           Sscperm_para(T_i,3) = dble(N_spin)/tmp_i * (SabsscP_para(T_i,3)**2 - (sum_SscP2_para(T_i,3)/dble(sweep2))**2)!scpermittivity_y
           !******************

           !*******write
           write(fo3+1000,*)T,sum_Szig(T_i,1),sum_Szig(T_i,2),sum_Szig(T_i,3),sum_Szig1(T_i,1),sum_Szig1(T_i,2)&
           ,sum_Szig1(T_i,3),sum_Szig2(T_i,1),sum_Szig2(T_i,2),sum_Szig2(T_i,3),sum_Szig3(T_i,1),sum_Szig3(T_i,2)&
           ,sum_Szig3(T_i,3),sum_Sstripe(T_i,1),sum_Sstripe(T_i,2),sum_Sstripe(T_i,3),sum_Sstripe1(T_i,1)&
           ,sum_Sstripe1(T_i,2),sum_Sstripe1(T_i,3),sum_Sstripe2(T_i,1),sum_Sstripe2(T_i,2),sum_Sstripe2(T_i,3)&
           ,sum_Sstripe3(T_i,1),sum_Sstripe3(T_i,2),sum_Sstripe3(T_i,3),Sk(T_i,1,:),Sk(T_i,2,:),Sk(T_i,3,:)&
           ,SscP(T_i,:),SabsscP(T_i,:),SxAF(T_i,:),SAF(T_i,:) ,Stm(T_i,:),StmP(T_i,:),SabstmP(T_i,:),Sscperm(T_i,:)&
           ,Stmperm(T_i,:),sum_Spm(T_i,:),sum_Sm(T_i,:),Sk(T_i,4,:),SabsscP_intr(T_i,:),SabsscP_para(T_i,:)&
           ,Sscperm_intr(T_i,:),Sscperm_para(T_i,:)
           !*******



        end do


     end do

     close(seed)
     close(fo2)
     close(fo3)
     close(fo3+1000)

!!$     allocate(spinA0(2,L,L,L))
!!$     allocate(spinB0(2,L,L,L))
!!$
!!$     spinA0 = spinA(:,:,:,:,T_ex(1))
!!$     spinB0 = spinB(:,:,:,:,T_ex(1))
!!$     m = 0.0d0
!!$     call calc_vec(spinA,spinB,1,T_ex,m)
!!$     m0 = m(1,:)
!!$
!!$     Ht = 0.0d0
!!$     H=0.0d0
!!$     H(2)=0.1d0
!!$     h0 = 1.0d0
!!$     !***********************file****************************
!!$     write(fname3(seed),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a)') 'HsweepL',L,'x',int(H_ini(1)),'y',int(H_ini(2)),'seed',seed,'K.d' !output format(system size,magnetic field,seed)
!!$     fname3(seed) =  adjustl(fname3(seed))
!!$     open(seed, file = fname3(seed) ,action ='write')
!!$
!!$     write(fname4(seed),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a)') 'Hsweep2L',L,'x',int(H_ini(1)),'y',int(H_ini(2)),'seed',seed,'K.d' !output format(system size,magnetic field,seed)
!!$     fname4(seed) =  adjustl(fname4(seed))
!!$     open(fo2, file = fname4(seed) ,action ='write')
!!$
!!$     open(fo3, file = 'all_tmPx.d' ,action='write')
!!$     !*******************************************************
!!$
!!$     !***********************H_sweep*******************************
!!$     do i = 0,H_mesh*H_speed
!!$
!!$
!!$        !***********totoidal************************
!!$        tmP=0.0d0
!!$        tm=0.0d0
!!$        m=0.0d0
!!$        AF=0.0d0
!!$        E=0.0d0
!!$        T_i=1
!!$        call sysE(rg,spinA,spinB,T_i,T_ex,J_ij,E)
!!$        call calc_vec(spinA,spinB,T_i,T_ex,m)
!!$        call troidal(rg,spinA,spinB,T_i,T_ex,m,tm,tmP)
!!$        call staggered(spinA,spinB,T_i,T_ex,AF)
!!$        call ac(rg,spinA,spinB,spinA0,spinB0,m0,T_i,T_ex,m,spin_ac,toro_ac)
!!$        write(*,'(100f)')H,dble(i),tm,AF,abs(tm),E
!!$        write(fo3,'(100f)')dble(i),abs(spin_ac),abs(toro_ac)
!!$        !*******************************************
!!$
!!$        if(mod(i,H_speed)==0.0d0) then
!!$          !***********totoidal************************
!!$           tmP=0.0d0
!!$        tm=0.0d0
!!$        m=0.0d0
!!$        AF=0.0d0
!!$        E=0.0d0
!!$        T_i=1
!!$        call calc_vec(spinA,spinB,T_i,T_ex,m)
!!$        call troidal(rg,spinA,spinB,T_i,T_ex,m,tm,tmP)
!!$
!!$        !*******************************************
!!$        end if
!!$
!!$        if(mod(i,H_speed)==0)then
!!$
!!$           !if(i/=1) then
!!$           ! H(1)=H_ini(1) - (i/H_speed-1)*(2.0d0*H_ini(1)/dble(H_mesh))
!!$           ! H(2)=H_ini(2) - (i/H_speed-1)*(2.0d0*H_ini(2)/dble(H_mesh))
!!$           !end if
!!$           !*************average and write****************************************
!!$           rmsm(T_i) = sqrt(sum_rmsm(T_i)/dble(sweep2))
!!$           C(T_i) = (dble(N_spin)/tmp_i**2)*(sum_e(T_i,2)/dble(sweep2) - (sum_e(T_i,1)/dble(sweep2))**2)
!!$           E(T_i) = sum_e(T_i,1)/dble(sweep2) !energy-density
!!$           k(T_i,1) = sum_k(T_i,1)/dble(sweep2)!chiral-bond1
!!$           k(T_i,2) = sum_k(T_i,2)/dble(sweep2)!chiral-bond2
!!$           k(T_i,3) = dble(N_spin)/tmp_i*(rmsm(T_i)**2 - (sum_rmsm2(T_i)/dble(sweep2))**2) !Susceptibility
!!$           AF(T_i,1) = sqrt(sum_AF(T_i,3)/dble(sweep2))!staggered_x
!!$           AF(T_i,2) = sqrt(sum_AF(T_i,4)/dble(sweep2))!staggered_y
!!$           xAF(T_i,1) = dble(N_spin)/tmp_i * (sum_AF(T_i,3)/dble(sweep2) - (sum_AF(T_i,1)/dble(sweep2))**2)!AFsusceptibity_x
!!$           xAF(T_i,2) = dble(N_spin)/tmp_i * (sum_AF(T_i,4)/dble(sweep2) - (sum_AF(T_i,2)/dble(sweep2))**2)!AFsusceptibity_y
!!$           tm(T_i,1) = sum_tm(T_i,1)/dble(sweep2)!troidalmoment-bond1
!!$           tm(T_i,2) = sum_tm(T_i,2)/dble(sweep2)!troidalmoment-bond2
!!$           tm(T_i,3) = sum_tm(T_i,3)/dble(sweep2)!troidalmoment-bond3
!!$           scP(T_i,1) = sum_scP(T_i,1)/dble(sweep2)!spincurrent-polarization_x
!!$           scP(T_i,2) = sum_scP(T_i,2)/dble(sweep2)!spincurrent-polarization_y
!!$           absscP(T_i,1) = sqrt(sum_scP(T_i,3)/dble(sweep2))!spincurrent-polarization_x-abs
!!$           absscP(T_i,2) = sqrt(sum_scP(T_i,4)/dble(sweep2))!spincurrent-polarization_y-abs
!!$           tmP(T_i,1) = sum_tmP(T_i,1)/dble(sweep2)!troidalmoment-polarization_x
!!$           tmP(T_i,2) = sum_tmP(T_i,2)/dble(sweep2)!troidalmoment-polarization_y
!!$           tmP(T_i,3) = sum_tmP(T_i,3)/dble(sweep2)!troidalmoment-polarization_z
!!$           abstmP(T_i,1) = sqrt(sum_tmP(T_i,4)/dble(sweep2))!troidalmoment-polarization_x-abs
!!$           abstmP(T_i,2) = sqrt(sum_tmP(T_i,5)/dble(sweep2))!troidalmoment-polarization_y-abs
!!$           abstmP(T_i,3) = sqrt(sum_tmP(T_i,6)/dble(sweep2))!troidalmoment-polarization_z-abs
!!$           scperm(T_i,1) = dble(N_spin)/tmp_i * (absscP(T_i,1)**2 - (sum_scP2(T_i,1)/dble(sweep2))**2)!scpermittivity_x
!!$           scperm(T_i,2) = dble(N_spin)/tmp_i * (absscP(T_i,2)**2 - (sum_scP2(T_i,2)/dble(sweep2))**2)!scpermittivity_y
!!$           tmperm(T_i,1) = dble(N_spin)/tmp_i * ( abstmP(T_i,1)**2 - (sum_tmP2(T_i,1)/dble(sweep2))**2)!tmpermittivity_x
!!$           tmperm(T_i,2) = dble(N_spin)/tmp_i * ( abstmP(T_i,2)**2 - (sum_tmP2(T_i,2)/dble(sweep2))**2)!tmpermittivity_y
!!$           tmperm(T_i,3) = dble(N_spin)/tmp_i * ( abstmP(T_i,3)**2 - (sum_tmP2(T_i,3)/dble(sweep2))**2)!tmpermittivity_z
!!$           X_SG_kmin(T_i,2) = sum_q_uv_kmin(T_i,1)/dble(sweep2) !coreration ratio guzai_L_xx
!!$           X_SG_kmin(T_i,3) = sum_q_uv_kmin(T_i,4)/dble(sweep2) !coreration ratio guzai_L_yy
!!$           X_SG(T_i,1) = X_SG(T_i,1)/dble(sweep2)!spinoverlap-tensor__xx
!!$           X_SG(T_i,2) = X_SG(T_i,2)/dble(sweep2)!spinoverlap-tensor__xy
!!$           X_SG(T_i,3) = X_SG(T_i,3)/dble(sweep2)!spinoverlap-tensor__yx
!!$           X_SG(T_i,4) = X_SG(T_i,4)/dble(sweep2)!spinoverlap-tensor__yy
!!$           X_SG(T_i,5) = sumX_SG(T_i,1)/dble(sweep2)!sum SG_X
!!$           X_SG_kmin(T_i,1) = sumX_SG(T_i,2)/dble(sweep2)
!!$           X_CG(T_i,1) = sum_q_k(T_i,1)/dble(sweep2)
!!$           X_CG(T_i,2) = sum_q_k(T_i,2)/dble(sweep2)
!!$           X_CG(T_i,3) = sum_q_k_kmin(T_i,1)/dble(sweep2)
!!$           X_CG(T_i,4) = sum_q_k_kmin(T_i,2)/dble(sweep2)
!!$
!!$           !**************************write physical quantities***************************
!!$           write(seed,*)H,E(T_i),C(T_i),rmsm(T_i),k(T_i,1),k(T_i,2),k(T_i,3),scP(T_i,1),scP(T_i,2),&
!!$                absscP(T_i,1),absscP(T_i,2),xAF(T_i,1),xAF(T_i,2),AF(T_i,1),AF(T_i,2),&
!!$                tm(T_i,1),tm(T_i,2),tm(T_i,3),tmP(T_i,1),tmP(T_i,2),tmP(T_i,3),&
!!$                abstmP(T_i,1),abstmP(T_i,2),abstmP(T_i,3),scperm(T_i,1),scperm(T_i,2),&
!!$                tmperm(T_i,1),tmperm(T_i,2),tmperm(T_i,3),X_SG(T_i,1),X_SG(T_i,2),X_SG(T_i,3),X_SG(T_i,4),&
!!$                X_SG(T_i,5), X_SG_kmin(T_i,1), X_SG_kmin(T_i,2), X_SG_kmin(T_i,3),&
!!$                X_CG(T_i,1), X_CG(T_i,2), X_CG(T_i,3), X_CG(T_i,4)
!!$           !******************************************************************************
!!$
!!$           mono(T_i) = sqrt(sum_mono(T_i)/dble(sweep2))
!!$           quad(T_i,1) = sqrt(sum_quad(T_i,1)/dble(sweep2))
!!$           quad(T_i,2) = sqrt(sum_quad(T_i,2)/dble(sweep2))
!!$           quad(T_i,3) = sqrt(sum_quad(T_i,3)/dble(sweep2))
!!$           quad(T_i,4) = sqrt(sum_quad(T_i,4)/dble(sweep2))
!!$           quad(T_i,5) = sqrt(sum_quad(T_i,5)/dble(sweep2))
!!$           quad(T_i,6) = sqrt(sum_quad(T_i,6)/dble(sweep2))
!!$           quad(T_i,7) = sqrt(sum_quad(T_i,7)/dble(sweep2))
!!$           quad(T_i,8)  = sqrt(sum_quad(T_i,8)/dble(sweep2))
!!$           quad(T_i,9) = sqrt(sum_quad(T_i,9)/dble(sweep2))
!!$           AF(T_i,1) = (sum_vmAF(T_i,1)/dble(sweep2))!staggered_x
!!$           AF(T_i,2) = (sum_vmAF(T_i,2)/dble(sweep2))!staggered_y
!!$           m(T_i,1) = (sum_vm(T_i,1)/dble(sweep2))!staggered_x
!!$           m(T_i,2) = (sum_vm(T_i,2)/dble(sweep2))!staggered_y
!!$
!!$           !**************************write physical quantities***************************
!!$           write(fo2,*)H,mono(T_i),quad(T_i,1),quad(T_i,2),quad(T_i,3),quad(T_i,4),&
!!$                quad(T_i,5),quad(T_i,6),quad(T_i,7),quad(T_i,8),quad(T_i,9),&
!!$                AF(T_i,1),AF(T_i,2),m(T_i,1),m(T_i,2)
!!$           !******************************************************************************
!!$
!!$
!!$
!!$           !****************************************************************************
!!$
!!$           !*******************zero shift****************
!!$           sweep2 = 0
!!$           tmp_max = 1
!!$           tmp_min = N_T
!!$           counter = 0
!!$           sum_e =0.0d0
!!$           sum_rmsm = 0.0d0
!!$           sum_AF = 0.0d0
!!$           sum_k = 0.0d0
!!$           sum_tm = 0.0d0
!!$           sum_scP = 0.0d0
!!$           sum_scP2 = 0.0d0
!!$           sum_tmP = 0.0d0
!!$           sum_tmP2 = 0.0d0
!!$           sum_q_uv_kmin = 0.0d0
!!$           sumX_SG = 0.0d0
!!$           X_SG = 0.0d0
!!$           X_CG = 0.0d0
!!$           sum_q_k = 0.0d0
!!$           sum_q_k_kmin = 0.0d0
!!$           q_uv = 0.0d0
!!$           q_uv_kmin = 0.0d0
!!$           sum_rmsm2 = 0.0d0
!!$           sum_mono=0.0d0
!!$           sum_quad=0.0d0
!!$           sum_vmAF=0.0d0
!!$           sum_vm=0.0d0
!!$           !*********************************************
!!$        end if
!!$
!!$        if(mod(int(2*i/H_speed)+1,2)==0)then !dvide by N_sweep/2
!!$
!!$           !**zero shift**   to inistiarize observals
!!$           m = 0.0d0
!!$           AF = 0.0d0
!!$           k = 0.0d0
!!$           scP = 0.0d0
!!$           tm = 0.0d0
!!$           tmP = 0.0d0
!!$           q_uv = 0.0d0
!!$           q_uv_kmin = 0.0d0
!!$           q_k = 0.0d0
!!$           q_k_kmin = 0.0d0
!!$           mono=0.0d0
!!$           quad=0.0d0
!!$           !*************
!!$
!!$           sweep2 = sweep2 + 1
!!$
!!$           !$omp parallel do num_threads(N_core)
!!$           do T_i = 1,N_T
!!$              call calc_vec(spinA,spinB,T_i,T_ex,m)
!!$              call sysE(rg,spinA,spinB,T_i,T_ex,J_ij,E)
!!$              call staggered(spinA,spinB,T_i,T_ex,AF)
!!$              call calc_chila_pola(spinA,spinB,T_i,T_ex,k,scP)
!!$              ! call spin_overlap(spinA,spin_repA,spinB,spin_repB,T_i,T_ex,T_ex_rep,q_uv,q_uv_kmin)
!!$              !call chiral_overlap(spinA,spin_repA,spinB,spin_repB,T_i,T_ex,T_ex_rep,q_k,q_k_kmin)
!!$              call troidal(rg,spinA,spinB,T_i,T_ex,m,tm,tmP)
!!$              call monpole(rg,spinA,spinB,T_i,T_ex,m,mono)
!!$              call quadpole(rg,spinA,spinB,T_i,T_ex,m,quad)
!!$           end do
!!$           !$omp end parallel do
!!$
!!$           sum_e(1:N_T,1) = sum_e(1:N_T,1) + E(1:N_T)/dble(N_spin)
!!$           sum_e(1:N_T,2) = sum_e(1:N_T,2) + E(1:N_T)/dble(N_spin)*E/dble(N_spin)
!!$
!!$           do T_i = 1,N_T
!!$              sum_rmsm(T_i) = dot_product(m(T_i,:),m(T_i,:)) + sum_rmsm(T_i)
!!$              sum_rmsm2(T_i) =  sum_rmsm2(T_i) + sqrt(dot_product(m(T_i,:),m(T_i,:)))
!!$           end do
!!$
!!$           sum_k(:,1) = sum_k(:,1) + k(:,1)/dble(N_spin)!chiral-bond1
!!$           sum_k(:,2) = sum_k(:,2) + k(:,2)/dble(N_spin)!chiral-bond2
!!$           sum_k(:,3) = sum_k(:,3) + k(:,3)/dble(N_spin)!chiral-bond3
!!$           sum_AF(:,1) =  sum_AF(:,1) + abs(AF(:,1)/dble(N_spin))!staggered_x
!!$           sum_AF(:,2) =  sum_AF(:,2) + abs(AF(:,2)/dble(N_spin))!staggered_y
!!$           sum_AF(:,3) =  sum_AF(:,3) + (AF(:,1)/dble(N_spin))**2!staggered_x^2
!!$           sum_AF(:,4) =  sum_AF(:,4) + (AF(:,2)/dble(N_spin))**2!staggered_y^2
!!$           sum_tm(:,1) = sum_tm(:,1) + tm(:,1)/dble(N_spin)!troidal-moment_x
!!$           sum_tm(:,2) = sum_tm(:,2) + tm(:,2)/dble(N_spin)!troidal-moment_y
!!$           sum_tm(:,3) = sum_tm(:,3) + tm(:,3)/dble(N_spin)!troidal-moment_z
!!$           sum_scP(:,1) = sum_scP(:,1) + scP(:,1)/dble(N_spin)!sc-polarization_x
!!$           sum_scP(:,2) = sum_scP(:,2) + scP(:,2)/dble(N_spin)!sc-polarization_y
!!$           sum_scP(:,3) = sum_scP(:,3) + (scP(:,1)/dble(N_spin))**2!sc-polarization_x^2
!!$           sum_scP(:,4) = sum_scP(:,4) + (scP(:,2)/dble(N_spin))**2!sc-polarization_y^2
!!$           sum_scP2(:,1) = sum_scP2(:,1) + abs(scP(:,1)/dble(N_spin))!sc-polarization_x-abs
!!$           sum_scP2(:,2) = sum_scP2(:,2) + abs(scP(:,2)/dble(N_spin))!sc-polarization_y-abs
!!$           sum_tmP(:,1) = sum_tmP(:,1) + tmP(:,1)/dble(N_spin)!tm-pola_x
!!$           sum_tmP(:,2) = sum_tmP(:,2) + tmP(:,2)/dble(N_spin)!tm-pola_y
!!$           sum_tmP(:,3) = sum_tmP(:,3) + tmP(:,3)/dble(N_spin)!tm-pola_z
!!$           sum_tmP(:,4) = sum_tmP(:,4) + (tmP(:,1)/dble(N_spin))**2!tm-pola_x^2
!!$           sum_tmP(:,5) = sum_tmP(:,5) + (tmP(:,2)/dble(N_spin))**2!tm-pola_y^2
!!$           sum_tmP(:,6) = sum_tmP(:,6) + (tmP(:,3)/dble(N_spin))**2!tm-pola_z^2
!!$           sum_tmP2(:,1) = sum_tmP(:,1) + abs(tmP(:,1)/dble(N_spin))!tm-pola_x-abs
!!$           sum_tmP2(:,2) = sum_tmP(:,2) + abs(tmP(:,2)/dble(N_spin))!tm-pola_y-abs
!!$           sum_tmP2(:,3) = sum_tmP(:,3) + abs(tmP(:,3)/dble(N_spin))!tm-pola_z-abs
!!$
!!$
!!$
!!$           sum_mono = sum_mono + (mono/dble(N_spin*L))**2
!!$           sum_quad(:,1) = sum_quad(:,1) + (quad(:,1)/dble(N_spin))**2
!!$           sum_quad(:,2) = sum_quad(:,2) + (quad(:,2)/dble(N_spin))**2
!!$           sum_quad(:,3) = sum_quad(:,3) + (quad(:,3)/dble(N_spin))**2
!!$           sum_quad(:,4) = sum_quad(:,4) + (quad(:,4)/dble(N_spin))**2
!!$           sum_quad(:,5) = sum_quad(:,5) + (quad(:,5)/dble(N_spin))**2
!!$           sum_quad(:,6) = sum_quad(:,6) + (quad(:,6)/dble(N_spin))**2
!!$           sum_quad(:,7) = sum_quad(:,7) + (quad(:,7)/dble(N_spin))**2
!!$           sum_quad(:,8) = sum_quad(:,8) + (quad(:,8)/dble(N_spin))**2
!!$           sum_quad(:,9) = sum_quad(:,9) + (quad(:,9)/dble(N_spin))**2
!!$
!!$           sum_vmAF(:,:)=sum_vmAF(:,:)+AF/dble(N_spin)
!!$           sum_vm = sum_vm +m/dble(N_spin)
!!$
!!$        end if
!!$
!!$        call stock_random(randnum)
!!$
!!$
!!$        !$omp parallel do num_threads(N_core)
!!$        do T_i = 1, N_T
!!$
!!$           call metropolis(rg,spinA,spinB,J_ij,T_i,T,alpha,randnum)
!!$           call over_relaxation(rg,spinA,spinB,T_i,J_ij)
!!$
!!$        end do
!!$
!!$     end do
!!$
!!$
!!$     !******************end H_asweep*******************************

     deallocate(spinA)
     deallocate(spinB)
     !deallocate(spin_repA)
     !deallocate(spin_repB)
     deallocate(J_ij)
     deallocate(randnum)
     !deallocate(randnum_rep)
     !deallocate(spinA0)
     !deallocate(spinB0)


  end do
  !****************************end sample loop****************************

  !**********MPI*************
  !if(MPIswitch) then
  !   call mpi_finalize(ierr)
  !end if
  !**************************

  call statistics(dataname,fname,mesh,N_sample)
  !call boot_strap(fname,dataname2,mesh,N_sample,30,L,alpha,T_min)
  call statistics(dataname3,fname2,mesh,N_sample)
  call statistics(dataname4,fname3,mesh,N_sample)
  call statistics(dataname5,fname5,mesh,N_sample)


end program main

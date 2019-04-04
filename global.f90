module global
  implicit none
  !***********constant*************
  real(8),parameter :: pi = acos(-1.0d0)
  real(8),parameter :: twopi = 2.0d0*pi
  !********************************

  !**********switch**********************
  logical, parameter ::MPIswitch=.false.
  logical, parameter ::joreiswitch=.false.
  logical, parameter ::nitiswitch=.false.
  !**************************************

  !*********************parameter*********************
  integer,parameter :: L = 24 !system size
  integer,parameter :: N_bond = 3
  integer,parameter :: N_sweep = 50000 ,N_relax = L
  integer,parameter :: mesh = 18
  integer,parameter :: H_speed = 1 !montecarlo_sweep
  integer,parameter :: H_mesh = 1
  integer,parameter :: N_spin = 2*L*L*L !number of spin
  integer,parameter :: N_sample = 5 !number of sample
  integer,parameter :: N_sampn = 1!sample per node
  integer,parameter :: N_core=24
  integer,parameter :: N_node=1
  integer,parameter :: N_T = 1
  integer,parameter :: warm = 5000
  real(8),parameter :: HB_limit = 0d0
  real(8),parameter :: ave = 0.0d0,Niratio=0.5
  real(8),parameter :: T_max = 1.05d0 , T_min = 0.03d0
  real(8)  H(3) !magnetic field
  real(8)  Ht(3)  !electric field
  real(8) h0
  real(8),parameter:: cnst=1.0d0/sqrt(6.0d0)
  real(8),parameter ::  H_ini(3) = (/1.0d0*cnst , 1.0d0*cnst , -2.0d0*cnst/) !magnetic field
  real(8),parameter ::  Ht_ini(3) = (/0.0d0 , 0.0d0 , 0.0d0/) !electric field
  real(8),parameter :: Jp = 0.3d0
  real(8),parameter :: regJ = 1.0d0
  real(8),parameter :: J_0 = 1.0d0 !Dispersion
  real(8),parameter :: K_0 = 1.0d0 !kitaev
  real(8),parameter :: XY_0 = 0.0d0 !XYiho 0>:XY ani , 0<:ising ani
  real(8),parameter :: Z_0 = 1.0d0 !Ziho
  real(8),parameter :: gam = 0.0d0!gamma
  !***************************************************

  !********************variable***********************
  integer,save :: fo !output file
  integer,save :: seed
  integer,save :: samething
  !***************************************************

end module global

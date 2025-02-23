
   module tqsubs_module

      implicit none

   contains

! subroutine tild03                all systems                92/11/01
! portability : all systems
! 92/11/01 lu : original version
!
! purpose :
!  initial values of the variables for constrained minimization.
!  linearly constrained and dense version.
!
! parameters :
!  io  nf  number of variables.
!  io  nc  number of constraints.
!  ii  ncl  number of linear constraints.
!  ro  x(nf)  vector of variables.
!  io  ix(nf)  vector containing types of bounds.
!  ro  xl(nf)  vector containing lower bounds for variables.
!  ro  xu(nf)  vector containing upper bounds for variables.
!  io  ic(nc)  vector containing types of constraints.
!  ro  cl(nc)  vector containing lower bounds for constraint functions.
!  ro  cu(nc)  vector containing upper bounds for constraint functions.
!  ro  cg(nf*nc) matrix whose columns are normals of the linear
!         constraints.
!  ii  next  number of the test problem.
!  io  ierr  error indicator.
!
      subroutine tild03(nf,nc,x,ix,xl,xu,ic,cl,cu,cg,fmin,xmax,next, &
       ierr)
      integer nf,nc,ix(nf),ic(nc),next,iext,ierr
      double precision x(nf),xl(nf),xu(nf),cl(nc),cu(nc),cg(nf*nc), &
       fmin,xmax
      integer i,j,ncl
      double precision y(235)
      common /empr03/ y
      double precision eta9
      parameter (eta9=1.0d60)
      fmin=0.0d0
      xmax=1.0d3
      iext=0
      ierr=0
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150) next
   10 nf=2
      nc=1
      ncl=1
      x(1)=0.0d0
      x(2)=0.0d0
      ix(1)=0
      ix(2)=0
      ic(1)=5
      cl(1)=0.0d0
      cu(1)=0.0d0
      cg(1)=4.0d0
      cg(2)=-3.0d0
      fmin=-eta9
      return
   20 nf=2
      nc=2
      ncl=2
      x(1)=1.0d0
      x(2)=0.5d0
      ix(1)=1
      ix(2)=1
      xl(1)=0.0d0
      xl(2)=0.0d0
      ic(1)=1
      ic(2)=3
      cl(1)=0.0d0
      cl(2)=0.0d0
      cu(2)=6.0d0
      cg(1)=1.0d0/sqrt(3.0d0)
      cg(2)=-1.0d0
      cg(3)=1.0d0
      cg(4)=sqrt(3.0d0)
      fmin=-eta9
      return
   30 nf=3
      nc=1
      ncl=1
      x(1)=0.7d0
      x(2)=0.2d0
      x(3)=0.1d0
      do 31 i=1,nf
      ix(i)=3
      xl(i)=0.0d0
      xu(i)=1.0d0
      cg(i)=1.0d0
   31 continue
      ic(1)=5
      cl(1)=1.0d0
      cu(1)=1.0d0
      fmin=-eta9
      return
   40 nf=3
      nc=1
      ncl=1
      do 41 i=1,nf
      x(i)=1.0d1
      ix(i)=3
      xl(i)=0.0d0
      xu(i)=4.2d1
   41 continue
      ic(1)=3
      cl(1)=0.0d0
      cu(1)=7.2d1
      cg(1)=1.0d0
      cg(2)=2.0d0
      cg(3)=2.0d0
      fmin=-eta9
      return
   50 nf=4
      nc=1
      ncl=1
      do 51 i=1,nf
      x(i)=2.0d0
      ix(i)=3
      xl(i)=0.0d0
      xu(i)=1.0d0
   51 continue
      xu(4)=2.0d0
      ic(1)=5
      cl(1)=0.0d0
      cu(1)=0.0d0
      cg(1)=1.0d0
      cg(2)=2.0d0
      cg(3)=2.0d0
      cg(4)=-1.0d0
      return
   60 nf=4
      nc=1
      ncl=1
      x(1)=3.0d0
      x(2)=-1.0d0
      x(3)=0.0d0
      x(4)=1.0d0
      do 61 i=1,nf
      ix(i)=2
      xu(i)=2.0d1
      cg(i)=1.0d0
   61 continue
      ic(1)=1
      cl(1)=1.0d0
      return
   70 nf=5
      nc=2
      ncl=2
      x(1)=1.0d1
      x(2)=7.0d0
      x(3)=2.0d0
      x(4)=-3.0d0
      x(5)=0.8d0
      do 71 i=1,nf
      ix(i)=0
      cg(i)=1.0d0
      cg(i+nf)=0.0d0
   71 continue
      ic(1)=5
      ic(2)=5
      cl(1)=7.0d0
      cu(1)=7.0d0
      cl(2)=6.0d0
      cu(2)=6.0d0
      cg(4)=4.0d0
      cg(5)=0.0d0
      cg(8)=1.0d0
      cg(10)=5.0d0
      return
   80 nf=5
      nc=3
      ncl=3
      x(1)=3.5d1
      x(2)=-3.1d1
      x(3)=1.1d1
      x(4)=5.0d0
      x(5)=-5.0d0
      do 82 i=1,nf
      ix(i)=0
      do 81 j=1,nc
      cg((j-1)*nf+i)=0.0d0
   81 continue
   82 continue
      do 83 j=1,nc
      ic(j)=5
      cl(j)=6.0d0
      cu(j)=6.0d0
      cg((j-1)*nf+j)=1.0d0
      cg((j-1)*nf+j+1)=2.0d0
      cg((j-1)*nf+j+2)=3.0d0
   83 continue
      return
   90 nf=5
      nc=3
      ncl=3
      do 91 j=1,5
      x(j)=2.0d0
      ix(j)=0
   91 continue
      do 92 j=1,3
      ic(j)=5
   92 continue
      do 93 j=1,3
      cl(j)=0.0d0
      cu(j)=0.0d0
   93 continue
      do 94 j=1,15
      cg(j)=0.0d0
   94 continue
      cg(1)=1.0d0
      cg(2)=3.0d0
      cg(8)=1.0d0
      cg(9)=1.0d0
      cg(10)=-2.0d0
      cg(12)=1.0d0
      cg(15)=-1.0d0
      return
  100 nf=5
      nc=10
      ncl=10
      do 101 j=1,5
      x(j)=0.0d0
      ix(j)=1
      xl(j)=0.0d0
  101 continue
      x(5)=1.0d0
      do 102 j=1,10
      ic(j)=1
  102 continue
      do 103 j=1,50
      cg(j)=0.0d0
  103 continue
      cg(1)=-16.0d0
      cg(2)=2.0d0
      cg(4)=1.0d0
      cg(7)=-2.0d0
      cg(9)=4.0d0
      cg(10)=2.0d0
      cg(11)=-3.5d0
      cg(13)=2.0d0
      cg(17)=-2.0d0
      cg(19)=-4.0d0
      cg(20)=-1.0d0
      cg(22)=-9.0d0
      cg(23)=-2.0d0
      cg(24)=1.0d0
      cg(25)=-2.8d0
      cg(26)=2.0d0
      cg(28)=-4.0d0
      do 104 j=1,5
      cg(30+j)=-1.0d0
  104 continue
      cg(36)=-1.0d0
      cg(37)=-2.0d0
      cg(38)=-3.0d0
      cg(39)=-2.0d0
      cg(40)=-1.0d0
      cg(41)=1.0d0
      cg(42)=2.0d0
      cg(43)=3.0d0
      cg(44)=4.0d0
      cg(45)=5.0d0
      do 105 j=1,5
      cg(45+j)=1.0d0
  105 continue
      cl(1)=-40.0d0
      cl(2)=-2.0d0
      cl(3)=-0.25d0
      cl(4)=-4.0d0
      cl(5)=-4.0d0
      cl(6)=-1.0d0
      cl(7)=-40.0d0
      cl(8)=-60.0d0
      cl(9)=5.0d0
      cl(10)=1.0d0
      fmin=-eta9
      return
  110 nf=6
      nc=1
      ncl=1
      x(1)=6.0d3
      x(2)=1.5d0
      x(3)=4.0d6
      x(4)=2.0d0
      x(5)=3.0d-3
      x(6)=5.0d7
      do 111 j=1,6
      ix(j)=3
  111 continue
      do 112 j=1,6
      xl(j)=0.0d0
  112 continue
      xl(2)=-10.0d0
      xl(5)=-1.0d0
      xu(1)=2.0d4
      xu(2)=10.0d0
      xu(3)=1.0d7
      xu(4)=20.0d0
      xu(5)=1.0d0
      xu(6)=2.0d8
      ic(1)=5
      cl(1)=1.76d4
      cu(1)=cl(1)
      do 113 j=1,6
      cg(j)=0.0d0
  113 continue
      cg(1)=1.0d0
      cg(2)=4.0d3
      fmin=-eta9
      return
  120 nf=6
      nc=6
      ncl=6
      do 121 j=1,6
      x(j)=0.0d0
      ix(j)=1
      xl(j)=0.0d0
  121 continue
      ix(1)=3
      ix(4)=3
      xu(1)=1.0d0
      xu(4)=1.0d0
      x(1)=1.0d0
      x(2)=2.0d0
      x(6)=2.0d0
      do 122 j=1,6
      ic(j)=5
  122 continue
      cl(1)=6.0d0
      cl(2)=3.0d0
      cl(3)=2.0d0
      cl(4)=1.0d0
      cl(5)=2.0d0
      cl(6)=2.0d0
      do 123 j=1,6
      cu(j)=cl(j)
  123 continue
      do 124 j=1,36
      cg(j)=0.0d0
  124 continue
      cg(1)=1.0d0
      cg(2)=2.0d0
      cg(5)=5.0d0
      cg(7)=1.0d0
      cg(8)=1.0d0
      cg(9)=1.0d0
      cg(16)=1.0d0
      cg(17)=1.0d0
      cg(18)=1.0d0
      cg(19)=1.0d0
      cg(22)=1.0d0
      cg(26)=1.0d0
      cg(29)=1.0d0
      cg(33)=1.0d0
      cg(36)=1.0d0
      return
  130 nf=8
      nc=1
      ncl=1
      x(1)=0.1d0
      x(2)=0.2d0
      x(3)=100.0d0
      x(4)=125.0d0
      x(5)=175.0d0
      x(6)=11.2d0
      x(7)=13.2d0
      x(8)=15.8d0
      do 131 j=1,8
      ix(j)=3
  131 continue
      xl(1)=0.001d0
      xu(1)=0.499d0
      xl(2)=xl(1)
      xu(2)=xu(1)
      xl(3)=100.0d0
      xu(3)=180.0d0
      xl(4)=130.0d0
      xu(4)=210.0d0
      xl(5)=170.0d0
      xu(5)=240.0d0
      do 132 j=6,8
      xl(j)=5.0d0
      xu(j)=25.0d0
  132 continue
      ic(1)=2
      cu(1)=1.0d0
      do 133 j=3,8
      cg(j)=0.0d0
  133 continue
      cg(1)=1.0d0
      cg(2)=1.0d0
      y(1)=95.0d0
      y(2)=105.0d0
      do 1301  j=3,6
      y(j)=110.0d0
 1301 continue
      do 1302 j=7,10
      y(j)=115.0d0
 1302 continue
      do 1303 j=11,25
      y(j)=120.0d0
 1303 continue
      do 1304 j=26,40
      y(j)=125.0d0
 1304 continue
      do 1305 j=41,55
      y(j)=130.0d0
 1305 continue
      do 1306 j=56,68
 1306 y(j)=135.0d0
      continue
      do 1307 j=69,89
      y(j)=140.0d0
 1307 continue
      do 1308 j=90,101
      y(j)=145.0d0
 1308 continue
      do 1309 j=102,118
      y(j)=150.0d0
 1309 continue
      do 1310 j=119,122
      y(j)=155.0d0
 1310 continue
      do 1311 j=123,142
      y(j)=160.0d0
 1311 continue
      do 1312 j=143,150
      y(j)=165.0d0
 1312 continue
      do 1313 j=151,167
      y(j)=170.0d0
 1313 continue
      do 1314 j=168,175
      y(j)=175.0d0
 1314 continue
      do 1315 j=176,181
      y(j)=180.0d0
 1315 continue
      do 1316 j=182,187
      y(j)=185.0d0
 1316 continue
      do 1317 j=188,194
      y(j)=190.0d0
 1317 continue
      do 1318 j=195,198
      y(j)=195.0d0
 1318 continue
      do 1319 j=199,201
      y(j)=200.0d0
 1319 continue
      do 1320 j=202,204
      y(j)=205.0d0
 1320 continue
      do 1321 j=205,212
      y(j)=210.0d0
 1321 continue
      y(213)=215.0d0
      do 1322 j=214,219
      y(j)=220.0d0
 1322 continue
      do 1323 j=220,224
      y(j)=230.0d0
 1323 continue
      y(225)=235.0d0
      do 1324 j=226,232
      y(j)=240.0d0
 1324 continue
      y(233)=245.0d0
      do 1325 j=234,235
      y(j)=250.0d0
 1325 continue
      return
  140 nf=10
      nc=3
      ncl=3
      do 141 j=1,10
      x(j)=0.1d0
      ix(j)=1
      xl(j)=1.0d-6
  141 continue
      ic(1)=5
      ic(2)=5
      ic(3)=5
      cl(1)=2.0d0
      cl(2)=1.0d0
      cl(3)=1.0d0
      do 142 j=1,3
      cu(j)=cl(j)
  142 continue
      do 143 j=1,30
      cg(j)=0.0d0
  143 continue
      cg(1)=1.0d0
      cg(2)=2.0d0
      cg(3)=2.0d0
      cg(6)=1.0d0
      cg(10)=1.0d0
      cg(14)=1.0d0
      cg(15)=2.0d0
      cg(16)=1.0d0
      cg(17)=1.0d0
      cg(23)=1.0d0
      cg(27)=1.0d0
      cg(28)=1.0d0
      cg(29)=2.0d0
      cg(30)=1.0d0
      fmin=-eta9
      return
  150 nf=16
      nc=8
      ncl=8
      do 151 j=1,16
      x(j)=10.0d0
      ix(j)=3
      xl(j)=0.0d0
      xu(j)=5.0d0
  151 continue
      do 152 j=1,8
      ic(j)=5
  152 continue
      do 153 j=1,128
      cg(j)=0.0d0
  153 continue
      cg(1)=0.22d0
      cg(2)=0.20d0
      cg(3)=0.19d0
      cg(4)=0.25d0
      cg(5)=0.15d0
      cg(6)=0.11d0
      cg(7)=0.12d0
      cg(8)=0.13d0
      cg(9)=1.0d0
      cg(17)=-1.46d0
      cg(19)=-1.30d0
      cg(20)=1.82d0
      cg(21)=-1.15d0
      cg(23)=0.80d0
      cg(26)=1.0d0
      cg(33)=1.29d0
      cg(34)=-0.89d0
      cg(37)=-1.16d0
      cg(38)=-0.96d0
      cg(40)=-0.49d0
      cg(43)=1.0d0
      cg(49)=-1.10d0
      cg(50)=-1.06d0
      cg(51)=0.95d0
      cg(52)=-0.54d0
      cg(54)=-1.78d0
      cg(55)=-0.41d0
      cg(60)=1.0d0
      cg(68)=-1.43d0
      cg(69)=1.51d0
      cg(70)=0.59d0
      cg(71)=-0.33d0
      cg(72)=-0.43d0
      cg(77)=1.0d0
      cg(82)=-1.72d0
      cg(83)=-0.33d0
      cg(85)=1.62d0
      cg(86)=1.24d0
      cg(87)=0.21d0
      cg(88)=-0.26d0
      cg(94)=1.0d0
      cg(97)=1.12d0
      cg(100)=0.31d0
      cg(103)=1.12d0
      cg(105)=-0.36d0
      cg(111)=1.0d0
      cg(114)=0.45d0
      cg(115)=0.26d0
      cg(116)=-1.10d0
      cg(117)=0.58d0
      cg(119)=-1.03d0
      cg(120)=0.10d0
      cg(128)=1.0d0
      cl(1)=2.5d0
      cl(2)=1.1d0
      cl(3)=-3.1d0
      cl(4)=-3.5d0
      cl(5)=1.3d0
      cl(6)=2.1d0
      cl(7)=2.3d0
      cl(8)=-1.5d0
      do 154 j=1,8
      cu(j)=cl(j)
  154 continue
      return
      end subroutine tild03

! subroutine tffu03                all systems                92/11/01
! portability : all systems
! 92/11/01 lu : original v0ersion
!
! purpose :
!  values of model functions for constrained minimization.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ri  x(n)  vector of variables.
!  ro  f  value of the model function.
!  ii  next  number of the test problem.
!
      subroutine tffu03(n,x,f,next)
      integer n,next
      double precision x(n),f
      double precision a,b,c,h,souc,sq,hk,pi
      data pi /3.14159265358979323d0/
      integer i,j
      double precision cc(5,5),aa(16,16),ev(5),dv(5),y(235),cv(10)
      common /empr03/ y
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150) next
   10 a=pi
      f=sin(a*x(1)/1.2d1)*cos(a*x(2)/1.6d1)
      return
   20 f=((x(1)-3.0d0)**2-9.0d0)*x(2)**3/(2.7d1*sqrt(3.0d0))
      return
   30 a=x(3)+3.0d-2
      b=a+x(2)
      c=b+x(1)
      f=-32.174d0*(255.0d0*log(c/(b+9.0d-2*x(1)))+&
      0280.0d0*log(b/(a+7.0d-2*x(2)))+&
      0290.0d0*log(a/(3.0d-2+1.3d-1*x(3))))
      return
   40 f=-x(1)*x(2)*x(3)
      return
   50 f=2.0d0-x(1)*x(2)*x(3)
      return
   60 f=(x(1)+1.0d1*x(2))**2+5.0d0*(x(3)-x(4))**2+ &
        (x(2)-2.0d0*x(3))**4+1.0d1*(x(1)-x(4))**4
      return
   70 f=(x(1)-x(2))**2+(x(3)-1.0d0)**2+(x(4)-1.0d0)**4+ &
        (x(5)-1.0d0)**6
      return
   80 f=(x(1)-x(2))**2+(x(2)-x(3))**2+(x(3)-x(4))**4+(x(4)-x(5))**2
      return
   90 f=(x(1)-x(2))**2+(x(2)+x(3)-2.0d0)**2+ &
         (x(4)-1.0d0)**2+(x(5)-1.0d0)**2
      return
  100 ev(1)=-15.0d0
      ev(2)=-27.0d0
      ev(3)=-36.0d0
      ev(4)=-18.0d0
      ev(5)=-12.0d0
      cc(1,1)=30.0d0
      cc(1,2)=-20.0d0
      cc(1,3)=-10.0d0
      cc(1,4)=32.0d0
      cc(1,5)=-10.0d0
      cc(2,1)=-20.0d0
      cc(2,2)=39.0d0
      cc(2,3)=-6.0d0
      cc(2,4)=-31.0d0
      cc(2,5)=32.0d0
      cc(3,1)=-10.0d0
      cc(3,2)=-6.0d0
      cc(3,3)=10.0d0
      cc(3,4)=-6.0d0
      cc(3,5)=-10.0d0
      cc(4,1)=32.0d0
      cc(4,2)=-31.0d0
      cc(4,3)=-6.0d0
      cc(4,4)=39.0d0
      cc(4,5)=-20.0d0
      cc(5,1)=-10.0d0
      cc(5,2)=32.0d0
      cc(5,3)=-10.0d0
      cc(5,4)=-20.0d0
      cc(5,5)=30.0d0
      dv(1)=4.0d0
      dv(2)=8.0d0
      dv(3)=10.0d0
      dv(4)=6.0d0
      dv(5)=2.0d0
      f=0.0d0
      do 101 j=1,5
      do 102 i=1,5
      f=f+cc(i,j)*x(i)*x(j)
  102 continue
      f=f+ev(j)*x(j)+dv(j)*x(j)**3
  101 continue
      return
  110 hk=0.96d0*4.9d13
      h=((x(1)-1.0d4)**2/6.4d7+(x(1)-1.0d4)*(x(2)-1.0d0)/2.0d4+ &
        (x(2)-1.0d0)**2)*(x(3)-2.0d6)**2/hk+ &
        (x(4)-10.0d0)**2/2.5d3+(x(5)-1.0d-3)**2/2.5d-3+ &
        (x(6)-1.0d8)**2/2.5d17
      f=-exp(-h/2.0d0)
      return
  120 f=x(1)+2.0d0*x(2)+4.0d0*x(5)+exp(x(1)*x(4))
      return
  130 f=0.0d0
      sq=sqrt(2.0d0*pi)
      do 131 i=1,235
      a=x(1)/x(6)*exp(-(y(i)-x(3))**2/(2.0d0*x(6)**2))
      b=x(2)/x(7)*exp(-(y(i)-x(4))**2/(2.0d0*x(7)**2))
      c=(1.0d0-x(2)-x(1))/x(8)* &
           exp(-(y(i)-x(5))**2/(2.0d0*x(8)**2))
      f=f-log((a+b+c)/sq)
  131 continue
      return
  140 f=0.0d0
      cv(1)=-6.089d0
      cv(2)=-17.164d0
      cv(3)=-34.054d0
      cv(4)=-5.914d0
      cv(5)=-24.721d0
      cv(6)=-14.986d0
      cv(7)=-24.100d0
      cv(8)=-10.708d0
      cv(9)=-26.662d0
      cv(10)=-22.179d0
      souc=0.0d0
      do 141 j=1,10
      souc=souc+x(j)
  141 continue
      do 142 j=1,10
      f=f+x(j)*(cv(j)+log(x(j)/souc))
  142 continue
      return
  150 f=0.0d0
      do 151 i=1,16
      do 152 j=1,16
      aa(i,j)=0.0d0
  152 continue
  151 continue
      do 153 i=1,16
      aa(i,i)=1.0d0
  153 continue
      aa(1,4)=1.0d0
      aa(1,7)=1.0d0
      aa(1,8)=1.0d0
      aa(1,16)=1.0d0
      aa(2,3)=1.0d0
      aa(2,7)=1.0d0
      aa(2,10)=1.0d0
      aa(3,7)=1.0d0
      aa(3,9)=1.0d0
      aa(3,10)=1.0d0
      aa(3,14)=1.0d0
      aa(4,7)=1.0d0
      aa(4,11)=1.0d0
      aa(4,15)=1.0d0
      aa(5,6)=1.0d0
      aa(5,10)=1.0d0
      aa(5,12)=1.0d0
      aa(5,16)=1.0d0
      aa(6,8)=1.0d0
      aa(6,15)=1.0d0
      aa(7,11)=1.0d0
      aa(7,13)=1.0d0
      aa(8,10)=1.0d0
      aa(8,15)=1.0d0
      aa(9,12)=1.0d0
      aa(9,16)=1.0d0
      aa(10,14)=1.0d0
      aa(11,13)=1.0d0
      aa(12,14)=1.0d0
      aa(13,14)=1.0d0
      do 154 i=1,16
      do 155 j=1,16
      f=f+aa(i,j)*(x(i)**2+x(i)+1.0d0)*(x(j)**2+x(j)+1.0d0)
  155 continue
  154 continue
      return
      end subroutine tffu03

! subroutine tfgu03                all systems                92/11/01
! portability : all systems
! 92/11/01 lu : original version
!
! purpose :
!  gradients of model functions for constrained minimization.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ri  x(n)  vector of variables.
!  ro  g(n)  gradient of the model function.
!  ii  next  number of the test problem.
!
      subroutine tfgu03(n,x,g,next)
      integer n,next
      double precision x(n),g(n)
      double precision a,b,c,d,e,f,sou,aexp,bexp,cexp,hi,sq,h,hk,pi
      data pi /3.14159265358979323d0/
      integer i,j
      double precision cc(5,5),ev(5),dv(5),cv(10),y(235)
      double precision aa(16,16)
      common /empr03/ y
      go to (10,20,30,40,40,60,70,80,90,100,110,120,130,140,150) next
   10 a=pi
      b=a/1.2d1
      c=a/1.6d1
      g(1)=b*cos(b*x(1))*cos(c*x(2))
      g(2)=-c*sin(b*x(1))*sin(c*x(2))
      return
   20 a=2.7d1*sqrt(3.0d0)
      g(1)=2.0d0*(x(1)-3.0d0)*x(2)**3/a
      g(2)=((x(1)-3.0d0)**2-9.0d0)*3.0d0*x(2)**2/a
      return
   30 a=x(3)+3.0d-2
      b=a+x(2)
      c=b+x(1)
      d=1.3d-1*x(3)+3.0d-2
      e=a+7.0d-2*x(2)
      f=b+9.0d-2*x(1)
      g(1)=-32.174d0*255.0d0*(f-9.0d-2*c)/(f*c)
      g(2)=-32.174d0*(255.0d0*(f-c)/(f*c)+&
      0280.0d0*(e-7.0d-2*b)/(e*b))
      g(3)=-32.174d0*(255.0d0*(f-c)/(f*c)+ &
        280.0d0*(e-b)/(e*b)+290.0d0*(d-1.3d-1*a)/(d*a))
      return
   40 g(1)=-x(2)*x(3)
      g(2)=-x(1)*x(3)
      g(3)=-x(1)*x(2)
      return
   60 a=x(1)+1.0d1*x(2)
      b=x(3)-x(4)
      c=(x(2)-2.0d0*x(3))**3
      d=(x(1)-x(4))**3
      g(1)=2.0d0*a+4.0d1*d
      g(2)=2.0d1*a+4.0d0*c
      g(3)=1.0d1*b-8.0d0*c
      g(4)=-1.0d1*b-4.0d1*d
      return
   70 a=2.0d0*(x(1)-x(2))
      b=2.0d0*(x(3)-1.0d0)
      c=4.0d0*(x(4)-1.0d0)**3
      d=6.0d0*(x(5)-1.0d0)**5
      g(1)=a
      g(2)=-a
      g(3)=b
      g(4)=c
      g(5)=d
      return
   80 a=2.0d0*(x(1)-x(2))
      b=2.0d0*(x(2)-x(3))
      c=4.0d0*(x(3)-x(4))**3
      d=2.0d0*(x(4)-x(5))
      g(1)=a
      g(2)=b-a
      g(3)=c-b
      g(4)=d-c
      g(5)=-d
      return
   90 a=2.0d0*(x(1)-x(2))
      b=2.0d0*(x(2)+x(3)-2.0d0)
      g(1)=a
      g(2)=-a+b
      g(3)=b
      g(4)=2.0d0*(x(4)-1.0d0)
      g(5)=2.0d0*(x(5)-1.0d0)
      return
  100 ev(1)=-15.0d0
      ev(2)=-27.0d0
      ev(3)=-36.0d0
      ev(4)=-18.0d0
      ev(5)=-12.0d0
      cc(1,1)=30.0d0
      cc(1,2)=-20.0d0
      cc(1,3)=-10.0d0
      cc(1,4)=32.0d0
      cc(1,5)=-10.0d0
      cc(2,1)=-20.0d0
      cc(2,2)=39.0d0
      cc(2,3)=-6.0d0
      cc(2,4)=-31.0d0
      cc(2,5)=32.0d0
      cc(3,1)=-10.0d0
      cc(3,2)=-6.0d0
      cc(3,3)=10.0d0
      cc(3,4)=-6.0d0
      cc(3,5)=-10.0d0
      cc(4,1)=32.0d0
      cc(4,2)=-31.0d0
      cc(4,3)=-6.0d0
      cc(4,4)=39.0d0
      cc(4,5)=-20.0d0
      cc(5,1)=-10.0d0
      cc(5,2)=32.0d0
      cc(5,3)=-10.0d0
      cc(5,4)=-20.0d0
      cc(5,5)=30.0d0
      dv(1)=4.0d0
      dv(2)=8.0d0
      dv(3)=10.0d0
      dv(4)=6.0d0
      dv(5)=2.0d0
      do 105 j=1,5
      g(j)=0.0d0
  105 continue
      do 106 j=1,5
      do 107 i=1,5
      g(j)=g(j)+(cc(i,j)+cc(j,i))*x(i)
  107 continue
      g(j)=g(j)+ev(j)+3.0d0*dv(j)*x(j)**2
  106 continue
      return
  110 hk=0.96d0*4.9d13
      h=((x(1)-1.0d4)**2/6.4d7+(x(1)-1.0d4)*(x(2)-1.0d0)/2.0d4+ &
        (x(2)-1.0d0)**2)*(x(3)-2.0d6)**2/hk+ &
        (x(4)-10.0d0)**2/2.5d3+(x(5)-1.0d-3)**2/2.5d-3+ &
        (x(6)-1.0d8)**2/2.5d17
      d=exp(-h/2.0d0)/2.0d0
      g(1)=d*(2.0d0*(x(1)-1.0d4)/6.4d7+(x(2)-1.0d0)/2.0d4)* &
           (x(3)-2.0d6)**2/hk
      g(2)=d*((x(1)-1.0d4)/2.0d4+ &
           2.0d0*(x(2)-1.0d0))*(x(3)-2.0d6)**2/hk
      g(3)=d*((x(1)-1.0d4)**2/6.4d7+ &
           (x(1)-1.0d4)*(x(2)-1.0d0)/2.0d4+ &
           (x(2)-1.0d0)**2)*2.0d0*(x(3)-2.0d6)/hk
      g(4)=d*2.0d0*(x(4)-10.0d0)/2.5d3
      g(5)=d*2.0d0*(x(5)-1.0d-3)/2.5d-3
      g(6)=d*2.0d0*(x(6)-1.0d8)/2.5d17
      return
  120 g(1)=1.0d0+x(4)*exp(x(1)*x(4))
      g(2)=2.0d0
      g(3)=0.0d0
      g(4)=x(1)*exp(x(1)*x(4))
      g(5)=4.0d0
      g(6)=0.0d0
      return
  130 do 131 j=1,8
      g(j)=0.0d0
  131 continue
      sq=sqrt(2.0d0*pi)
      do 132 i=1,235
      aexp=exp(-(y(i)-x(3))**2/(2.0d0*x(6)**2))
      a=x(1)/x(6)*aexp
      bexp=exp(-(y(i)-x(4))**2/(2.0d0*x(7)**2))
      b=x(2)/x(7)*bexp
      cexp=exp(-(y(i)-x(5))**2/(2.0d0*x(8)**2))
      c=(1.0d0-x(2)-x(1))/x(8)*cexp
      hi=-sq/(a+b+c)
      g(1)=g(1)+hi*(aexp/x(6)-cexp/x(8))
      g(2)=g(2)+hi*(bexp/x(7)-cexp/x(8))
      g(3)=g(3)+hi*x(1)/x(6)*aexp*(y(i)-x(3))/(x(6)**2)
      g(4)=g(4)+hi*x(2)/x(7)*bexp*(y(i)-x(4))/(x(7)**2)
      g(5)=g(5)+hi*(1.0d0-x(2)-x(1))/x(8)*cexp*(y(i)-x(5))/(x(8)**2)
      g(6)=g(6)+hi*((-x(1)/(x(6)**2))*aexp+ &
         x(1)/x(6)*aexp*(y(i)-x(3))**2/(x(6)**3))
      g(7)=g(7)+hi*((-x(2)/(x(7)**2))*bexp+ &
         x(2)/x(7)*bexp*(y(i)-x(4))**2/(x(7)**3))
      g(8)=g(8)+hi*((x(1)+x(2)-1.0d0)/(x(8)**2)*cexp+ &
         (1.0d0-x(2)-x(1))/x(8)*cexp*(y(i)-x(5))**2/(x(8)**3))
  132 continue
      return
  140 cv(1)=-6.089d0
      cv(2)=-17.164d0
      cv(3)=-34.054d0
      cv(4)=-5.914d0
      cv(5)=-24.721d0
      cv(6)=-14.986d0
      cv(7)=-24.100d0
      cv(8)=-10.708d0
      cv(9)=-26.662d0
      cv(10)=-22.179d0
      sou=0.0d0
      do 142 j=1,10
      sou=sou+x(j)
  142 continue
      do 143 j=1,10
      g(j)=cv(j)+log(x(j)/sou)
  143 continue
      return
  150 do 151 i=1,16
      do 152 j=1,16
      aa(i,j)=0.0d0
  152 continue
  151 continue
      do 153 i=1,16
      aa(i,i)=1.0d0
  153 continue
      aa(1,4)=1.0d0
      aa(1,7)=1.0d0
      aa(1,8)=1.0d0
      aa(1,16)=1.0d0
      aa(2,3)=1.0d0
      aa(2,7)=1.0d0
      aa(2,10)=1.0d0
      aa(3,7)=1.0d0
      aa(3,9)=1.0d0
      aa(3,10)=1.0d0
      aa(3,14)=1.0d0
      aa(4,7)=1.0d0
      aa(4,11)=1.0d0
      aa(4,15)=1.0d0
      aa(5,6)=1.0d0
      aa(5,10)=1.0d0
      aa(5,12)=1.0d0
      aa(5,16)=1.0d0
      aa(6,8)=1.0d0
      aa(6,15)=1.0d0
      aa(7,11)=1.0d0
      aa(7,13)=1.0d0
      aa(8,10)=1.0d0
      aa(8,15)=1.0d0
      aa(9,12)=1.0d0
      aa(9,16)=1.0d0
      aa(10,14)=1.0d0
      aa(11,13)=1.0d0
      aa(12,14)=1.0d0
      aa(13,14)=1.0d0
      do 154 j=1,16
      g(j)=0.0d0
  154 continue
      do 155 j=1,16
      do 156 i=1,16
      g(j)=g(j)+ &
        (aa(i,j)+aa(j,i))*(x(i)**2+x(i)+1.0d0)*(2.0d0*x(j)+1.0d0)
  156 continue
  155 continue
      return
      end subroutine tfgu03

!
!     test subroutines for constrained optimization
!
! subroutine tind07             all systems                 90/12/01
! portability : all systems
! 90/12/01 lu : original version
!
! purpose :
!  initial values of the variables for nonlinear programming.
!  dense version.
!
! parameters :
!  io  n  number of variables.
!  io  nc  number of constraints.
!  io  ncl number of linear constraints.
!  ro  x(n)  vector of variables.
!  io  ix(nf)  vector containing types of bounds.
!  ro  xl(nf)  vector containing lower bounds for variables.
!  ro  xu(nf)  vector containing upper bounds for variables.
!  io  ic(nc)  types of constraints.
!  ro  cl(nc)  lower bounds of constraints.
!  ro  cu(nc)  upper bounds of constraints.
!  ro  fmin  lower bound for value of the objective function.
!  ro  xmax  maximum stepsize.
!  io  next  number of the test problem.
!  io  ierr  error indicator.
!
      subroutine tind07(n,nc,x,ix,xl,xu,ic,cl,cu,fmin,xmax,next, &
       ierr)
      integer i,n,nc,ncl,ix(n),ic(nc),next,ierr
      double precision x(n),xl(n),xu(n),cl(nc),cu(nc),fmin,xmax
      double precision y(128)
      common /empr07/ y
      ncl=0
      fmin=-1.0d60
      xmax=1.0d3
      ierr=0
      do 1 i=1,n
      ix(i)=0
      xl(i)=0.0d0
      xu(i)=0.0d0
    1 continue
      do 2 i=1,nc
      ic(i)=1
      cl(i)=0.0d0
      cu(i)=0.0d0
    2 continue
      go to(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170, &
       180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330, &
       340) next
   10 if(n.ge.2.and.nc.ge.1) then
      n=2
      nc=1
      x(1)=2.0d0
      x(2)=0.0d0
      ic(1)=2
      cu(1)=1.0d0
      else
      ierr=1
      endif
      return
   20 if(n.ge.2.and.nc.ge.1) then
      n=2
      nc=1
      x(1)=-10.0d0
      x(2)= 10.0d0
      ic(1)=2
      cu(1)=1.0d0
      else
      ierr=1
      endif
      return
   30 if(n.ge.2.and.nc.ge.1) then
      n=2
      nc=1
      x(1)=0.0d0
      x(2)=0.0d0
      ic(1)=2
      cu(1)=25.0d0
      else
      ierr=1
      endif
      return
   40 if(n.ge.2.and.nc.ge.1) then
      n=2
      nc=1
      x(1)=-2.0d0
      x(2)=-2.0d0
      ix(1)=1
      ix(2)=1
      else
      ierr=1
      endif
      return
   50 if(n.ge.2.and.nc.ge.1) then
      n=2
      nc=1
      x(1)=4.2d-1
      x(2)=5.0d0
      ix(1)=1
      ix(2)=1
      xl(1)= 0.4d0
      xl(2)=-4.0d0
      cl(1)=0.09d0
      y(1)=8.0d0
      y(2)=8.0d0
      y(3)=10.0d0
      y(4)=10.0d0
      y(5)=10.0d0
      y(6)=10.0d0
      y(7)=12.0d0
      y(8)=12.0d0
      y(9)=12.0d0
      y(10)=12.0d0
      y(11)=14.0d0
      y(12)=14.0d0
      y(13)=14.0d0
      y(14)=16.0d0
      y(15)=16.0d0
      y(16)=16.0d0
      y(17)=18.0d0
      y(18)=18.0d0
      y(19)=20.0d0
      y(20)=20.0d0
      y(21)=20.0d0
      y(22)=22.0d0
      y(23)=22.0d0
      y(24)=22.0d0
      y(25)=24.0d0
      y(26)=24.0d0
      y(27)=24.0d0
      y(28)=26.0d0
      y(29)=26.0d0
      y(30)=26.0d0
      y(31)=28.0d0
      y(32)=28.0d0
      y(33)=30.0d0
      y(34)=30.0d0
      y(35)=30.0d0
      y(36)=32.0d0
      y(37)=32.0d0
      y(38)=34.0d0
      y(39)=36.0d0
      y(40)=36.0d0
      y(41)=38.0d0
      y(42)=38.0d0
      y(43)=40.0d0
      y(44)=42.0d0
      y(45)=0.49d0
      y(46)=0.49d0
      y(47)=0.48d0
      y(48)=0.47d0
      y(49)=0.48d0
      y(50)=0.47d0
      y(51)=0.46d0
      y(52)=0.46d0
      y(53)=0.45d0
      y(54)=0.43d0
      y(55)=0.45d0
      y(56)=0.43d0
      y(57)=0.43d0
      y(58)=0.44d0
      y(59)=0.43d0
      y(60)=0.43d0
      y(61)=0.46d0
      y(62)=0.45d0
      y(63)=0.42d0
      y(64)=0.42d0
      y(65)=0.43d0
      y(66)=0.41d0
      y(67)=0.41d0
      y(68)=0.40d0
      y(69)=0.42d0
      y(70)=0.40d0
      y(71)=0.40d0
      y(72)=0.41d0
      y(73)=0.40d0
      y(74)=0.41d0
      y(75)=0.41d0
      y(76)=0.40d0
      y(77)=0.40d0
      y(78)=0.40d0
      y(79)=0.38d0
      y(80)=0.41d0
      y(81)=0.40d0
      y(82)=0.40d0
      y(83)=0.41d0
      y(84)=0.38d0
      y(85)=0.40d0
      y(86)=0.40d0
      y(87)=0.39d0
      y(88)=0.39d0
      else
      ierr=1
      endif
      return
   60 if(n.ge.2.and.nc.ge.2) then
      n=2
      nc=2
      x(1)=-2.0d0
      x(2)= 1.0d0
      ix(1)=3
      ix(2)=2
      xl(1)=-0.5d0
      xu(1)= 0.5d0
      xu(2)= 1.0d0
      else
      ierr=1
      endif
      return
   70 if(n.ge.2.and.nc.ge.2) then
      n=2
      nc=2
      x(1)=20.1d0
      x(2)=5.84d0
      ix(1)=3
      ix(2)=3
      xl(1)=1.3d1
      xu(1)=1.0d2
      xl(2)=0.0d0
      xu(2)=1.0d2
      ic(2)=2
      cl(1)=1.0d2
      cu(2)=82.81d0
      else
      ierr=1
      endif
      return
   80 if(n.ge.2.and.nc.ge.3) then
      n=2
      nc=3
      x(1)=-2.0d0
      x(2)= 1.0d0
      ix(1)=3
      xl(1)=-0.5d0
      xu(1)= 0.5d0
      cl(3)= 1.0d0
      else
      ierr=1
      endif
      return
   90 if(n.ge.2.and.nc.ge.3) then
      n=2
      nc=3
      x(1)=90.0d0
      x(2)=10.0d0
      ix(1)=3
      ix(2)=3
      xl(1)=0.0d0
      xu(1)=75.0d0
      xl(2)=0.0d0
      xu(2)=65.0d0
      cl(1)=700.0d0
      else
      ierr=1
      endif
      return
  100 if(n.ge.3.and.nc.ge.1) then
      n=3
      nc=1
      x(1)=1.0d0
      x(2)=1.0d0
      x(3)=1.0d0
      ic(1)=2
      cu(1)=48.0d0
      else
      ierr=1
      endif
      return
  110 if(n.ge.3.and.nc.ge.1) then
      n=3
      nc=1
      x(1)=1.0d0
      x(2)=1.0d0
      x(3)=1.0d0
      ix(1)=3
      ix(2)=3
      ix(3)=3
      xl(1)= 1.0d0
      xu(1)= 1.0d1
      xl(2)=-1.0d1
      xu(2)= 1.0d1
      xl(3)=-1.0d1
      xu(3)= 1.0d1
      cl(1)= 1.0d0
      else
      ierr=1
      endif
      return
  120 if(n.ge.3.and.nc.ge.1) then
      n=3
      nc=1
      x(1)=1.0d0
      x(2)=1.0d0
      x(3)=1.0d0
      ix(1)=3
      ix(2)=3
      ix(3)=3
      xl(1)=-1.0d1
      xu(1)= 1.0d1
      xl(2)= 1.0d0
      xu(2)= 1.0d1
      xl(3)=-1.0d1
      xu(3)= 1.0d0
      cl(1)= 1.0d0
      else
      ierr=1
      endif
      return
  130 if(n.ge.3.and.nc.ge.1) then
      n=3
      nc=1
      x(1)=1.0d0
      x(2)=1.0d0
      x(3)=1.0d0
      ix(1)=1
      ix(2)=1
      ix(3)=1
      xl(1)=1.0d-5
      xl(2)=1.0d-5
      xl(3)=1.0d-5
      ic(1)=2
      cu(1)=1.0d0
      else
      ierr=1
      endif
      return
  140 if(n.ge.3.and.nc.ge.2) then
      n=3
      nc=2
      x(1)=0.0d0
      x(2)=1.05d0
      x(3)=2.9d0
      ix(1)=3
      ix(2)=3
      ix(3)=3
      xl(1)=0.0d0
      xu(1)=1.0d2
      xl(2)=0.0d0
      xu(2)=1.0d2
      xl(3)=0.0d0
      xu(3)=1.0d1
      ic(1)=2
      ic(2)=2
      else
      ierr=1
      endif
      return
  150 if(n.ge.3.and.nc.ge.7) then
      n=3
      nc=7
      x(1)=1745.0d0
      x(2)=12000.0d0
      x(3)=110.0d0
      ix(1)=3
      ix(2)=3
      ix(3)=3
      xl(1)=1.0d-5
      xu(1)=2.0d3
      xl(2)=1.0d-5
      xu(2)=1.6d4
      xl(3)=1.0d-5
      xu(3)=1.2d2
      do 151 i=1,7
      ic(i)=3
  151 continue
      cl(1)=0.0d0
      cl(2)=0.0d0
      cl(3)=8.5d1
      cl(4)=9.0d1
      cl(5)=3.0d0
      cl(6)=1.0d-2
      cl(7)=1.45d2
      cu(1)=5.0d3
      cu(2)=2.0d3
      cu(3)=9.3d1
      cu(4)=9.5d1
      cu(5)=1.2d1
      cu(6)=4.0d0
      cu(7)=1.62d2
      else
      ierr=1
      endif
      return
  160 if(n.ge.4.and.nc.ge.1) then
      n=4
      nc=1
      x(1)=2.0d0
      x(2)=4.0d0
      x(3)=4.0d-2
      x(4)=2.0d0
      do 161 i=1,4
      ix(i)=3
      xl(i)=1.0d-5
      xu(i)=1.0d2
  161 continue
      xu(3)=1.0d0
      do 162 i=2,19
      y(i)=dble(i-1)
      y(i+38)=log(y(i))
  162 continue
      y(1)=0.1d0
      y(1+38)=log(y(1))
      y(20)=0.00189d0
      y(21)=0.1038d0
      y(22)=0.268d0
      y(23)=0.506d0
      y(24)=0.577d0
      y(25)=0.604d0
      y(26)=0.725d0
      y(27)=0.898d0
      y(28)=0.947d0
      y(29)=0.845d0
      y(30)=0.702d0
      y(31)=0.528d0
      y(32)=0.385d0
      y(33)=0.257d0
      y(34)=0.159d0
      y(35)=0.0869d0
      y(36)=0.0453d0
      y(37)=0.01509d0
      y(38)=0.00189d0
      else
      ierr=1
      endif
      return
  170 if(n.ge.4.and.nc.ge.2) then
      n=4
      nc=2
      x(1)=1.0d0
      x(2)=1.0d0
      x(3)=1.0d0
      x(4)=1.0d0
      do 171 i=1,4
      ix(i)=3
      xl(i)=1.0d-3
      xu(i)=dble(5-i)*1.0d5
  171 continue
      ic(1)=2
      ic(2)=2
      cu(1)=0.0401d0
      cu(2)=0.010085d0
      else
      ierr=1
      endif
      return
  180 if(n.ge.4.and.nc.ge.3) then
      n=4
      nc=3
      do 181 i=1,n
      x(i)=0.0d0
  181 continue
      ic(1)=2
      ic(2)=2
      ic(3)=2
      cu(1)=8.0d0
      cu(2)=1.0d1
      cu(3)=5.0d0
      else
      ierr=1
      endif
      return
  190 if(n.ge.5.and.nc.ge.3) then
      n=5
      nc=3
      x(1)=1.0d0
      x(2)=1.0d0
      x(3)=1.0d0
      x(4)=1.0d0
      x(5)=1.0d0
      ic(1)=2
      cu(1)= 2.0d1
      cl(2)=-2.0d0
      cl(3)= 5.0d0
      else
      ierr=1
      endif
      return
  200 if(n.ge.5.and.nc.ge.3) then
      n=5
      nc=3
      x(1)=78.0d0
      x(2)=33.0d0
      x(3)=27.0d0
      x(4)=27.0d0
      x(5)=27.0d0
      do 201 i=1,5
      ix(i)=3
      xl(i)=27.0d0
      xu(i)=45.0d0
  201 continue
      xl(1)=78.0d0
      xl(2)=33.0d0
      xu(1)=10.2d1
      ic(1)=3
      ic(2)=3
      ic(3)=3
      cl(1)=0.0d0
      cu(1)=9.2d1
      cl(2)=9.0d1
      cu(2)=1.1d2
      cl(3)=2.0d1
      cu(3)=2.5d1
      else
      ierr=1
      endif
      return
  210 if(n.ge.5.and.nc.ge.3) then
      n=5
      nc=3
      x(1)=2.52d0
      x(2)=2.0d0
      x(3)=37.5d0
      x(4)=9.25d0
      x(5)=6.8d0
      do 211 i=1,5
      ix(i)=3
  211 continue
      xl(1)=0.0d0
      xu(1)=1.0d3
      xl(2)=1.2d0
      xu(2)=2.4d0
      xl(3)=2.0d1
      xu(3)=6.0d1
      xl(4)=9.0d0
      xu(4)=9.3d0
      xl(5)=6.5d0
      xu(5)=7.0d0
      ic(1)=3
      ic(2)=3
      ic(3)=3
      cl(1)=0.0d0
      cu(1)=294.0d3
      cl(2)=0.0d0
      cu(2)=294.0d3
      cl(3)=0.0d0
      cu(3)=277.2d3
      else
      ierr=1
      endif
      return
  220 if(n.ge.5.and.nc.ge.21) then
      n=5
      nc=21
      x(1)=900.0d0
      x(2)=80.0d0
      x(3)=115.0d0
      x(4)=267.0d0
      x(5)=27.0d0
      do 221 i=1,5
      ix(i)=3
  221 continue
      xl(1)=704.4148d0
      xu(1)=906.3855d0
      xl(2)=68.6d0
      xu(2)=288.88d0
      xl(3)=0.0d0
      xu(3)=134.75d0
      xl(4)=193.0d0
      xu(4)=287.0966d0
      xl(5)=25.0d0
      xu(5)=84.1988d0
      ic(3)=2
      cu(3)=21.0d0
      cl(4)=110.6d0
      do 222 i=5,21
      ic(i)=3
  222 continue
      cl(5)= 213.1d0
      cl(6)= 17.505d0
      cl(7)= 11.275d0
      cl(8)= 214.228d0
      cl(9)= 7.458d0
      cl(10)= 0.961d0
      cl(11)= 1.612d0
      cl(12)= 0.146d0
      cl(13)= 107.99d0
      cl(14)= 922.693d0
      cl(15)= 926.832d0
      cl(16)= 18.766d0
      cl(17)= 1072.163d0
      cl(18)= 8961.448d0
      cl(19)= 0.063d0
      cl(20)= 71084.33d0
      cl(21)= 2.802713d6
      cu(5)= 405.23d0
      cu(6)= 1053.6667d0
      cu(7)= 35.03d0
      cu(8)= 665.585d0
      cu(9)= 584.463d0
      cu(10)= 265.916d0
      cu(11)= 7.046d0
      cu(12)= 0.222d0
      cu(13)= 273.366d0
      cu(14)= 1286.105d0
      cu(15)= 1444.046d0
      cu(16)= 537.141d0
      cu(17)= 3247.039d0
      cu(18)= 26844.086d0
      cu(19)= 0.386d0
      cu(20)= 1.4d5
      cu(21)= 1.2146108d7
      else
      ierr=1
      endif
      return
  230 if(n.ge.6.and.nc.ge.2) then
      n=6
      nc=2
      x(1)=1.0d0
      x(2)=1.0d0
      x(3)=1.0d0
      x(4)=3.0d0
      x(5)=0.0d0
      x(6)=0.5d0
      ix(6)=3
      xl(6)=4.0d0
      xu(6)=8.0d0
      ic(1)=2
      ic(2)=2
      cu(1)=5.0d0
      cu(2)=1.0d0
      else
      ierr=1
      endif
      return
  240 if(n.ge.6.and.nc.ge.2) then
      n=6
      nc=2
      x(1)=5.54d0
      x(2)=4.4d0
      x(3)=12.02d0
      x(4)=11.82d0
      x(5)=0.702d0
      x(6)=0.852d0
      do 241 i=1,6
      ix(i)=1
  241 continue
      ic(2)=2
      cl(1)=2.07d0
      cu(2)=1.0d0
      else
      ierr=1
      endif
      return
  250 if(n.ge.6.and.nc.ge.4) then
      n=6
      nc=4
      do 251 i=1,6
      x(i)=0.0d0
      ix(i)=3
  251 continue
      xu(1)=0.31d0
      xu(2)=0.046d0
      xu(3)=0.068d0
      xu(4)=0.042d0
      xu(5)=0.028d0
      xu(6)=0.134d-1
      cl(1)=32.97d0
      cl(2)=25.12d0
      cl(3)=-124.08d0
      cl(4)=-173.02d0
      else
      ierr=1
      endif
      return
  260 if(n.ge.7.and.nc.ge.4) then
      n=7
      nc=4
      x(1)=1.0d0
      x(2)=2.0d0
      x(3)=0.0d0
      x(4)=4.0d0
      x(5)=0.0d0
      x(6)=1.0d0
      x(7)=1.0d0
      ic(1)=2
      ic(2)=2
      ic(3)=2
      ic(4)=2
      cu(1)=127.0d0
      cu(2)=282.0d0
      cu(3)=196.0d0
      else
      ierr=1
      endif
      return
  270 if(n.ge.7.and.nc.ge.5) then
      n=7
      nc=5
      do 271 i=1,7
      x(i)=6.0d0
      ix(i)=3
      xl(i)=1.0d-1
      xu(i)=1.0d1
  271 continue
      xl(7)=1.0d-2
      do 272 i=1,4
      ic(i)=2
      cu(i)=1.0d0
  272 continue
      ic(5)=3
      cl(5)=1.0d2
      cu(5)=3.0d3
      else
      ierr=1
      endif
      return
  280 if(n.ge.8.and.nc.ge.5) then
      n=8
      nc=5
      x(1)=6.0d0
      x(2)=3.0d0
      x(3)=0.4d0
      x(4)=0.2d0
      x(5)=6.0d0
      x(6)=6.0d0
      x(7)=1.0d0
      x(8)=0.5d0
      do 281 i=1,8
      ix(i)=3
      xl(i)=1.0d-1
      xu(i)=1.0d1
  281 continue
      do 282 i=1,4
      ic(i)=2
      cu(i)=1.0d0
  282 continue
      ic(5)=3
      cl(5)=1.0d0
      cu(5)=4.2d0
      else
      ierr=1
      endif
      return
  290 if(n.ge.8.and.nc.ge.6) then
      n=8
      nc=6
      do 291 i=1,8
      x(i)=5.0d3
      ix(i)=3
      xl(i)=1.0d1
      xu(i)=1.0d3
  291 continue
      x(4)=200.0d0
      x(5)=350.0d0
      x(6)=150.0d0
      x(7)=225.0d0
      x(8)=425.0d0
      xl(1)=1.0d2
      xl(2)=1.0d3
      xl(3)=1.0d3
      xu(1)=1.0d4
      xu(2)=1.0d4
      xu(3)=1.0d4
      do 292 i=1,3
      ic(i)=2
      cu(i)=1.0d0
  292 continue
      else
      ierr=1
      endif
      return
  300 if(n.ge.9.and.nc.ge.13) then
      n=9
      nc=13
      do 301 i=1,9
      x(i)=1.0d0
      ix(i)=1
      xl(i)=0.0d0
      ic(i)=2
      cu(i)=1.0d0
  301 continue
      ix(9)=1
      xl(9)=0.0d0
      ic(13)=2
      else
      ierr=1
      endif
      return
  310 if(n.ge.10.and.nc.ge.8) then
      n=10
      nc=8
      x(1)=2.0d0
      x(2)=3.0d0
      x(3)=5.0d0
      x(4)=5.0d0
      x(5)=1.0d0
      x(6)=2.0d0
      x(7)=7.0d0
      x(8)=3.0d0
      x(9)=6.0d0
      x(10)=1.0d1
      do 311 i=1,8
      ic(i)=2
  311 continue
      cu(1)=120.0d0
      cu(2)=40.0d0
      cu(3)=30.0d0
      cu(4)=0.0d0
      cu(5)=105.0d0
      cu(6)=0.0d0
      cu(7)=0.0d0
      cu(8)=12.0d0
      else
      ierr=1
      endif
      return
  320 if(n.ge.15.and.nc.ge.5) then
      n=15
      nc=5
      do 321 i=1,15
      x(i)=0.001d0
      ix(i)=1
      xl(i)=0.0d0
  321 continue
      x(7)=60.0d0
      cl(1)=15.0d0
      cl(2)=27.0d0
      cl(3)=36.0d0
      cl(4)=18.0d0
      cl(5)=12.0d0
      do 322 i=1,90
      y(i)=0.0d0
  322 continue
      y(1)=-16.0d0
      y(3)=-3.5d0
      y(6)= 2.0d0
      y(7)=-1.0d0
      y(8)=-1.0d0
      y(9)= 1.0d0
      y(10)= 1.0d0
      y(11)= 2.0d0
      y(12)=-2.0d0
      y(14)=-2.0d0
      y(15)=-9.0d0
      y(17)=-1.0d0
      y(18)=-2.0d0
      y(19)= 2.0d0
      y(20)= 1.0d0
      y(23)= 2.0d0
      y(25)=-2.0d0
      y(26)=-4.0d0
      y(27)=-1.0d0
      y(28)=-3.0d0
      y(29)= 3.0d0
      y(30)= 1.0d0
      y(31)= 1.0d0
      y(32)= 4.0d0
      y(34)=-4.0d0
      y(35)= 1.0d0
      y(37)=-1.0d0
      y(38)=-2.0d0
      y(39)= 4.0d0
      y(40)= 1.0d0
      y(42)= 2.0d0
      y(44)=-1.0d0
      y(45)=-2.8d0
      y(47)=-1.0d0
      y(48)=-1.0d0
      y(49)= 5.0d0
      y(50)= 1.0d0
      y(51)=-4.0d1
      y(52)=-2.0d0
      y(53)=-2.5d-1
      y(54)=-4.0d0
      y(55)=-4.0d0
      y(56)=-1.0d0
      y(57)=-4.0d1
      y(58)=-6.0d1
      y(59)= 5.0d0
      y(60)= 1.0d0
      y(61)= 3.0d1
      y(62)=-2.0d1
      y(63)=-1.0d1
      y(64)= 3.2d1
      y(65)=-1.0d1
      y(66)=-2.0d1
      y(67)= 3.9d1
      y(68)=-6.0d0
      y(69)=-3.1d1
      y(70)= 3.2d1
      y(71)=-1.0d1
      y(72)=-6.0d0
      y(73)= 1.0d1
      y(74)=-6.0d0
      y(75)=-1.0d1
      y(76)= 3.2d1
      y(77)=-3.1d1
      y(78)=-6.0d0
      y(79)= 3.9d1
      y(80)=-2.0d1
      y(81)=-1.0d1
      y(82)= 3.2d1
      y(83)=-1.0d1
      y(84)=-2.0d1
      y(85)= 3.0d1
      y(86)= 4.0d0
      y(87)= 8.0d0
      y(88)= 1.0d1
      y(89)= 6.0d0
      y(90)= 2.0d0
      else
      ierr=1
      endif
      return
  330 if(n.ge.16.and.nc.ge.19) then
      n=16
      nc=19
      x(1)=0.8d0
      x(2)=0.83d0
      x(3)=0.85d0
      x(4)=0.87d0
      x(5)=0.9d0
      x(6)=0.1d0
      x(7)=0.12d0
      x(8)=0.19d0
      x(9)=0.25d0
      x(10)=0.29d0
      x(11)=512.0d0
      x(12)=13.1d0
      x(13)=71.8d0
      x(14)=640.0d0
      x(15)=650.0d0
      x(16)=5.7d0
      do 331 i=1,16
      ix(i)=3
  331 continue
      do 332 i=1,10
      xl(i)=0.1d0
      xu(i)=0.9d0
  332 continue
      xl(5)=0.9d0
      xu(5)=1.0d0
      xl(6)=1.0d-4
      xu(6)=0.1d0
      xl(11)=1.0d0
      xu(11)=1.0d3
      xl(12)=1.0d-6
      xu(12)=5.0d2
      xl(13)=1.0d0
      xu(13)=5.0d2
      do 333 i=14,15
      xl(i)=5.0d2
      xu(i)=1.0d3
  333 continue
      xl(16)=1.0d-6
      xu(16)=5.0d2
      do 334 i=1,19
      ic(i)=2
      cu(i)=1.0d0
  334 continue
      do 335 i=1,5
      y(i)=1.262626d0
      y(i+5)=-1.231060d0
      y(i+25)=1.0d0
      y(i+34)=0.002d0
      y(i+43)=1.d0
      y(i+51)=1.d0
      y(i+56)=1.d0
  335 continue
      do 336 i=11,23,3
      y(i)=0.03475d0
      y(i+1)=0.975d0
      y(i+2)=-0.00975d0
  336 continue
      y(28)=-y(28)
      y(30)=0.002d0
      y(31)=y(30)
      y(32)=-y(30)
      y(33)=-y(30)
      y(34)=1.0d0
      y(37)=1.0d0
      y(38)=-y(38)
      y(39)=-y(39)
      y(40)=1.0d0
      y(41)=y(40)
      y(42)=5.0d2
      y(43)=-y(42)
      y(44)=-y(44)
      y(47)=y(42)
      y(48)=-y(48)
      y(49)=y(43)
      y(50)=0.9d0
      y(51)=0.002d0
      y(52)=-y(51)
      y(53)=y(51)
      y(54)=y(52)
      else
      ierr=1
      endif
      return
  340 if(n.ge.20.and.nc.ge.17) then
      n=20
      nc=17
      x(1)=2.0d0
      x(2)=3.0d0
      x(3)=5.0d0
      x(4)=5.0d0
      x(5)=1.0d0
      x(6)=2.0d0
      x(7)=7.0d0
      x(8)=3.0d0
      x(9)=6.0d0
      x(10)=1.0d1
      x(11)=2.0d0
      x(12)=2.0d0
      x(13)=6.0d0
      x(14)=1.5d1
      x(15)=1.0d0
      x(16)=2.0d0
      x(17)=1.0d0
      x(18)=2.0d0
      x(19)=1.0d0
      x(20)=3.0d0
      do 341 i=1,17
      ic(i)=2
  341 continue
      cu(1)=120.0d0
      cu(2)=40.0d0
      cu(3)=30.0d0
      cu(4)=0.0d0
      cu(5)=105.0d0
      cu(6)=0.0d0
      cu(7)=0.0d0
      cu(8)=12.0d0
      cu(9)=0.0d0
      cu(10)=28.0d0
      cu(11)=87.0d0
      cu(12)=10.0d0
      cu(13)=92.0d0
      cu(14)=54.0d0
      cu(15)=68.0d0
      cu(16)=-19.0d0
      cu(17)=0.0d0
      else
      ierr=1
      endif
      return
      end subroutine tind07

! subroutine tcfu07             all systems                 91/12/01
! portability : all systems
! 91/12/01 lu : original version
!
! purpose :
!  values of test functions for nonlinear approximation.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ii  kc  index of the constraint.
!  ri  x(n)  vector of variables.
!  ro  fc value of the constraint at the selected point.
!  ii  next  number of the test problem.
!
      subroutine tcfu07(n,kc,x,fc,next)
      integer n,kc,next
      double precision x(n),fc
      double precision x1,x2,x3,x4,x5,x6,x7,x8
      integer i,j,k
      double precision y(128)
      common /empr07/ y
      go to(10,20,30,40,50,60,70,80,90,100,10,91,130,140,150,160,170, &
       180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330, &
       340),next
   10 fc=x(1)**2+x(2)**2
      return
   20 fc=3.0d0*x(1)**2-2.0d0*x(1)*x(2)+x(2)**2
      return
   30 fc=4.0d0*x(1)**2+x(2)**2
      return
   40 fc=(1.0d0-x(1))**3-x(2)
      return
   50 fc=0.49d0*x(2)-x(1)*x(2)
      return
   60 go to (61,62),kc
   61 fc=x(2)**2-x(1)
      return
   62 fc=x(1)**2-x(2)
      return
   70 go to (71,72),kc
   71 fc=(x(1)-5.0d0)**2+(x(2)-5.0d0)**2
      return
   72 fc=(x(1)-6.0d0)**2+(x(2)-5.0d0)**2
      return
   80 go to (81,82,10),kc
   81 fc=x(1)+x(2)**2
      return
   82 fc=x(1)**2+x(2)
      return
   90 go to (91,92,93),kc
   91 fc=x(1)*x(2)
      return
   92 fc=x(2)-x(1)**2/125.0d0
      return
   93 fc=(x(2)-5.0d1)**2-5.0d0*(x(1)-5.5d1)
      return
  100 fc=x(1)**2+2.0d0*x(2)**2+4.0d0*x(3)**2
      return
  130 fc=4.0d0/x(1)+32.0d0/x(2)+120.0d0/x(3)
      return
  140 go to (141,142),kc
  141 fc=exp(x(1))-x(2)
      return
  142 fc=exp(x(2))-x(3)
      return
  150 x2=1.6d0*x(1)
  151 x3=1.22d0*x2-x(1)
      x6=(x(2)+x3)/x(1)
      x1=0.01d0*x(1)*(112.0d0+13.167d0*x6-0.6667d0*x6**2)
      if (abs(x1-x2).gt.0.001d0) then
      x2=x1
      go to 151
      endif
      x4=93.0d0
  152 x5=86.35d0+1.098d0*x6-0.038d0*x6**2+0.325d0*(x4-89.0d0)
      x8=3.0d0*x5-133.0d0
      x7=35.82d0-0.222d0*x8
      x1=98000.0d0*x(3)/(x2*x7+1000.0d0*x(3))
      if (abs(x1-x4).gt.0.001d0) then
      x4=x1
      go to 152
      endif
      go to (153,154,155,156,157,158,159),kc
  153 fc=x2
      return
  154 fc=x3
      return
  155 fc=x4
      return
  156 fc=x5
      return
  157 fc=x6
      return
  158 fc=x7
      return
  159 fc=x8
      return
  160 fc=x(3)+(1.0d0-x(3))*x(4)
      return
  170 go to (171,172),kc
  171 fc=4.0d0/x(1)+2.25d0/x(2)+1.0d0/x(3)+0.25d0/x(4)
      return
  172 fc=0.16d0/x(1)+0.36d0/x(2)+0.64d0/x(3)+0.64d0/x(4)
      return
  180 x1=x(1)*x(1)
      x2=x(2)*x(2)
      x3=x(3)*x(3)
      x4=x(4)*x(4)
      x5=x(1)+x(1)
      x6=x(2)+x(2)
      x7=x(3)+x(3)
      x8=x(4)+x(4)
      go to (181,182,183),kc
  181 fc=x1+x2+x3+x4+x(1)-x(2)+x(3)-x(4)
      return
  182 fc=x1+x2+x3+x4+x4-x(1)-x(4)
      return
  183 fc=x1+x1+x2+x3+x5-x(2)-x(4)
      return
  190 go to (191,192,193),kc
  191 fc=x(1)**2+x(2)**2+x(3)**2+x(4)**2+x(5)**2
      return
  192 fc=x(1)**2*x(3)+x(4)*x(5)
      return
  193 fc=x(2)**2*x(4)+10.0d0*x(1)*x(5)
      return
  200 go to (201,202,203),kc
  201 fc=85.334407d0+0.0056868d0*x(2)*x(5)+0.0006262d0*x(1)*x(4)- &
       0.0022053d0*x(3)*x(5)
      return
  202 fc=80.51249d0+0.0071317d0*x(2)*x(5)+0.0029955d0*x(1)*x(2)+ &
       0.0021813d0*x(3)**2
      return
  203 fc=9.300961d0+0.0047026d0*x(3)*x(5)+0.0012547d0*x(1)*x(3)+ &
       0.0019085d0*x(3)*x(4)
      return
  210 go to (211,212,213),kc
  211 fc=-145421.4020d0*x(1)+2931.15060d0*x(1)*x(2)- &
       40.4279320d0*x(1)*x(3)+5106.1920d0*x(1)*x(4)+ &
       15711.36d0*x(1)*x(5)
      return
  212 fc=-155011.1084d0*x(1)+4360.53352d0*x(1)*x(2)+ &
       12.9492344d0*x(1)*x(3)+10236.884d0*x(1)*x(4)+ &
       13176.786d0*x(1)*x(5)
      return
  213 fc=-326669.5104d0*x(1)+7390.68412d0*x(1)*x(2)- &
       27.8986976d0*x(1)*x(3)+16643.076d0*x(1)*x(4)+ &
       30988.146d0*x(1)*x(5)
      return
  220 continue
      x1=146.312d3
      y(1)=x(2)+x(3)+41.6d0
      y(18)=0.024d0*x(4)-4.62d0
      y(2)=12.5d0/y(18)+12.0d0
      y(19)=0.0003535d0*x(1)**2+0.5311d0*x(1)+0.08705d0*y(2)*x(1)
      y(20)=0.052d0*x(1)+78.0d0+0.002377d0*y(2)*x(1)
      y(3)=y(19)/y(20)
      y(4)=19.0d0*y(3)
      y(21)=0.04782d0*(x(1)-y(3))+0.1956d0*(x(1)-y(3))**2/x(2)+ &
       0.6376d0*y(4)+1.594d0*y(3)
      y(22)=100.0d0*x(2)
      y(23)=x(1)-y(3)-y(4)
      y(24)=0.95d0-y(21)/y(22)
      y(5)=y(23)*y(24)
      y(6)=x(1)-y(5)-y(4)-y(3)
      y(25)=(y(5)+y(4))*0.995d0
      y(7)=y(25)/y(1)
      y(8)=y(25)/3798.0d0
      y(26)=y(7)-0.0663d0*y(7)/y(8)-0.3153d0
      y(9)=96.82d0/y(26)+0.321d0*y(1)
      y(10)=1.29d0*y(5)+1.258d0*y(4)+2.29d0*y(3)+1.71d0*y(6)
      y(11)=1.71d0*x(1)-0.452d0*y(4)+0.58d0*y(3)
      y(27)=12.3d0/752.3d0
      y(28)=1.75d0*y(2)*0.995d0*x(1)
      y(29)=0.995d0*y(10)+1998.0d0
      y(12)=y(27)*x(1)+y(28)/y(29)
      y(13)=y(29)-1.75d0*y(2)
      y(14)=3623.0d0+64.4d0*x(2)+58.4d0*x(3)+x1/(y(9)+ &
       x(5))
      y(30)=0.995d0*y(10)+60.8d0*x(2)+48.0d0*x(4)-0.1121d0*y(14)- &
       5095.0d0
      y(15)=y(13)/y(30)
      y(16)=148.0d3-331.0d3*y(15)+40.0d0*y(13)-61.0d0*y(15)*y(13)
      y(31)=2324.0d0*y(10)-2.874d7*y(2)
      y(17)=1.413d7-1328.0d0*y(10)-531.0d0*y(11)+y(31)/y(29)
      y(32)=y(13)/y(15)-y(13)/0.52d0
      y(33)=1.104d0-0.72d0*y(15)
      y(34)=y(9)+x(5)
      if(kc.gt.4) go to 225
      go to (221,222,223,224),kc
  221 fc=y(4)-0.28d0*y(5)/0.72d0
      return
  222 fc=1.5d0*x(2)-x(3)
      return
  223 fc=3496.0d0*y(2)/y(29)
      return
  224 fc=62212.0d0/y(34)-y(1)
      return
  225 fc=y(kc-4)
      return
  230 go to (231,232),kc
  231 fc=x(1)**2+x(2)**2+x(3)**2
      return
  232 fc=x(5)**2+(x(4)-3.0d0)**2
      return
  240 go to (241,242),kc
  241 fc=0.001d0*x(1)*x(2)*x(3)*x(4)*x(5)*x(6)
      return
  242 fc=0.00062d0*x(1)*x(4)*x(5)**2*(x(1)+x(2)+x(3))+ &
       0.00058d0*x(2)*x(3)*x(6)**2*(x(1)+1.57d0*x(2)+x(4))
      return
  250 go to (251,252,253,254),kc
  251 fc=17.1d0*x(1)+38.2d0*x(2)+204.2d0*x(3)+212.3d0*x(4)+ &
       623.4d0*x(5)+1495.5d0*x(6)-169.0d0*x(1)*x(3)- &
       3580.0d0*x(3)*x(5)-3810.0d0*x(4)*x(5)-18500.0d0*x(4)*x(6)- &
       24300.0d0*x(5)*x(6)
      return
  252 fc=17.9d0*x(1)+36.8d0*x(2)+113.9d0*x(3)+169.7d0*x(4)+ &
       337.8d0*x(5)+1385.2d0*x(6)-139.0d0*x(1)*x(3)- &
       2450.0d0*x(4)*x(5)-16600.0d0*x(4)*x(6)-17200.0d0*x(5)*x(6)
      return
  253 fc=-273.0d0*x(2)-70.0d0*x(4)-819.0d0*x(5)+ &
       26000.0d0*x(4)*x(5)
      return
  254 fc=159.9d0*x(1)-311.0d0*x(2)+587.0d0*x(4)+391.0d0*x(5)+ &
       2198.0d0*x(6)-14000.0d0*x(1)*x(6)
      return
  260 go to (261,262,263,264),kc
  261 fc=2.0d0*x(1)**2+3.0d0*x(2)**4+x(3)+4.0d0*x(4)**2+ &
      5.0d0*x(5)
      return
  262 fc=7.0d0*x(1)+3.0d0*x(2)+1.0d1*x(3)**2+x(4)-x(5)
      return
  263 fc=2.3d1*x(1)+x(2)**2+6.0d0*x(6)**2-8.0d0*x(7)
      return
  264 fc=4.0d0*x(1)**2+x(2)**2-3.0d0*x(1)*x(2)+2.0d0*x(3)**2+ &
      5.0d0*x(6)-1.1d1*x(7)
      return
  270 go to (271,272,273,274,275),kc
  271 fc=0.5d0*x(1)**0.5d0*x(7)/(x(3)*x(6)**2)+ &
       0.7d0*x(1)**3*x(2)*x(6)*x(7)**0.5d0/x(3)**2+ &
       0.2d0*x(3)*x(6)**(2.0d0/3.0d0)*x(7)**0.25d0/ &
       (x(2)*x(4)**0.5d0)
      return
  272 fc=1.3d0*x(2)*x(6)/(x(1)**0.5d0*x(3)*x(5))+ &
       0.8d0*x(3)*x(6)**2/(x(4)*x(5))+ &
       3.1d0*x(2)**0.5d0*x(6)**(1.0d0/3.0d0)/ &
       (x(1)*x(4)**2*x(5))
      return
  273 fc=2.0d0*x(1)*x(5)*x(7)**(4.0d0/3.0d0)/(x(3)**1.5d0*x(6))+ &
       0.1d0*x(2)*x(5)/(x(3)**0.5d0*x(6)*x(7)**0.5d0)+ &
       x(2)*x(3)**0.5d0*x(5)/x(1)+ &
       0.65d0*x(3)*x(5)*x(7)/(x(2)**2*x(6))
      return
  274 fc=0.2d0*x(2)*x(5)**0.5d0*x(7)**(1.0d0/3.0d0)/ &
       (x(1)**2*x(4))+ &
       0.3d0*x(1)**0.5d0*x(2)**2*x(3)*x(4)**(1.0d0/3.0d0)* &
       x(7)**0.25d0/x(5)**(2.0d0/3.0d0)+ &
       0.4d0*x(3)*x(5)*x(7)**0.75d0/(x(1)**3*x(2)**2)+ &
       0.5d0*x(4)*x(7)**0.5d0/x(3)**2
      return
  275 fc=10.0d0*x(1)*x(4)**2/(x(2)*x(6)**3*x(7)**0.25d0)+ &
       15.0d0*x(3)*x(4)/(x(1)*x(2)**2*x(5)*x(7)**0.5d0)+ &
       20.0d0*x(2)*x(6)/(x(1)**2*x(4)*x(5)**2)+ &
       25.0d0*x(1)**2*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6)**2)
      return
  280 go to (281,282,283,284,285),kc
  281 fc=0.0588d0*x(5)*x(7)+0.1d0*x(1)
      return
  282 fc=0.0588d0*x(6)*x(8)+0.1d0*(x(1)+x(2))
      return
  283 fc=4.0d0*x(3)/x(5)+2.0d0/(x(3)**0.71d0*x(5))+ &
       0.0588d0*x(7)/x(3)**1.3d0
      return
  284 fc=4.0d0*x(4)/x(6)+2.0d0/(x(4)**0.71d0*x(6))+ &
       0.0588d0*x(8)/x(4)**1.3d0
      return
  285 fc=0.4d0*(x(1)/x(7))**0.67d0+0.4d0*(x(2)/x(8))**0.67d0+ &
       10.0d0-x(1)-x(2)
      return
  290 go to (291,292,293,294,295,296),kc
  291 fc=0.0025d0*(x(4)+x(6))
      return
  292 fc=0.0025d0*(x(5)+x(7)-x(4))
      return
  293 fc=0.01d0*(x(8)-x(5))
      return
  294 fc=x(1)*x(6)-833.33252d0*x(4)-100.0d0*x(1)+83333.333d0
      return
  295 fc=x(2)*x(7)-1250.0d0*x(5)-x(2)*x(4)+1250.0d0*x(4)
      return
  296 fc=x(3)*x(8)-1250000.0d0-x(3)*x(5)+2500.0d0*x(5)
      return
  300 go to (301,302,303,304,305,306,307,308,309,110,111,112,113),kc
  301 fc=x(3)**2+x(4)**2
      return
  302 fc=x(5)**2+x(6)**2
      return
  303 fc=(x(1)-x(5))**2+(x(2)-x(6))**2
      return
  304 fc=(x(1)-x(7))**2+(x(2)-x(8))**2
      return
  305 fc=(x(3)-x(5))**2+(x(4)-x(6))**2
      return
  306 fc=(x(3)-x(7))**2+(x(4)-x(8))**2
      return
  307 fc=x(7)**2+(x(8)-x(9))**2
      return
  308 fc=x(1)**2+(x(2)-x(9))**2
      return
  309 fc=x(9)**2
      return
  110 fc=x(3)*x(9)
      return
  111 fc=x(5)*x(8)-x(6)*x(7)
      return
  112 fc=x(1)*x(4)-x(2)*x(3)
      return
  113 fc=x(5)*x(9)
      return
  310 go to (311,312,313,314,315,316,317,318),kc
  311 fc=3.0d0*(x(1)-2.0d0)**2+4.0d0*(x(2)- &
      3.0d0)**2+2.0d0*x(3)**2-7.0d0*x(4)
      return
  312 fc=5.0d0*x(1)**2+8.0d0*x(2)+(x(3)-6.0d0)**2- &
      2.0d0*x(4)
      return
  313 fc=0.5d0*(x(1)-8.0d0)**2+2.0d0*(x(2)- &
      4.0d0)**2+3.0d0*x(5)**2-x(6)
      return
  314 fc=x(1)**2+2.0d0*(x(2)-2.0d0)**2- &
      2.0d0*x(1)*x(2)+1.4d1*x(5)-6.0d0*x(6)
      return
  315 fc=4.0d0*x(1)+5.0d0*x(2)-3.0d0*x(7)+ &
      9.0d0*x(8)
      return
  316 fc=1.0d1*x(1)-8.0d0*x(2)-1.7d1*x(7)+ &
      2.0d0*x(8)
      return
  317 fc=6.0d0*x(2)-3.0d0*x(1)+1.2d1*(x(9)- &
      8.0d0)**2-7.0d0*x(10)
      return
  318 fc=2.0d0*x(2)-8.0d0*x(1)+5.0d0*x(9)- &
      2.0d0*x(10)
      return
  320 fc=3.0d0*y(kc+85)*x(kc+10)**2
      k=(kc-1)*5
      do 321 j=1,5
      fc=fc+2.0d0*y(k+j+60)*x(j+10)
  321 continue
      k=(kc-1)*10
      do 322 i=1,10
      fc=fc-y(k+i)*x(i)
  322 continue
      return
  330 if (kc.lt.6) then
      k=11+3*(kc-1)
      fc=(y(k)+y(k+1)*x(kc+5)+y(k+2)*x(kc))*x(kc)/x(kc+5)
      return
      else if (kc.lt.14) then
      go to (331,332,333,334,335,336,337,338),kc-5
  331 fc=(y(26)*x(6)+(y(27)*x(1)+y(28)*x(6))*x(12)/x(11))/x(7)
      return
  332 fc=(y(29)*x(7)+(y(30)*x(7)+y(33)*x(1))*x(12)+y(31)*x(2)*x(13))/ &
       x(8)+y(32)*x(13)
      return
  333 fc=y(34)*x(8)+y(37)*x(9)+(y(35)*x(8)+y(38)*x(2))*x(13)+ &
       (y(36)*x(3)+y(39)*x(9))*x(14)
      return
  334 fc=((y(40)*x(14)+y(43))*x(9)+(y(41)*x(4)+y(44)*x(8))*x(15)+ &
       y(42)*x(10))/(x(3)*x(14))
      return
  335 fc=(y(45)*x(5)*x(16)+y(49)*x(10))/(x(4)*x(15))+ &
       y(46)*x(10)/x(4)+(y(47)+y(48)*x(16))/x(15)
      return
  336 fc=(y(50)+y(52)*x(5)*x(16))/x(4)+y(51)*x(16)
      return
  337 fc=y(53)*x(11)+y(54)*x(12)
      return
  338 fc=y(55)*x(12)/x(11)
      return
      else if (kc.lt.18) then
      i=kc-13
      k=5-i
      fc=y(55+i)*x(k)/x(k+1)
      return
      else if (kc.eq.18) then
      fc=y(60)*x(9)/x(10)
      return
      else
      fc=y(61)*x(8)/x(9)
      return
      endif
  340 go to (311,312,313,314,315,316,317,318,341,342,343,344,345,346, &
       347,348,349),kc
  341 fc=x(1)+x(2)+4.0d0*x(11)-2.1d1*x(12)
      return
  342 fc=x(1)**2+1.5d1*x(11)-8.0d0*x(12)
      return
  343 fc=4.0d0*x(1)+9.0d0*x(2)+5.0d0*x(13)**2-9.0d0*x(14)
      return
  344 fc=3.0d0*x(1)+4.0d0*x(2)+3.0d0*(x(13)- &
      6.0d0)**2-1.4d1*x(14)
      return
  345 fc=1.4d1*x(1)**2+3.5d1*x(15)-7.9d1*x(16)
      return
  346 fc=1.5d1*x(2)**2+1.1d1*x(15)-6.1d1*x(16)
      return
  347 fc=5.0d0*x(1)**2+2.0d0*x(2)+9.0d0*x(17)**4-x(18)
      return
  348 fc=x(1)**2-x(2)+1.9d1*x(19)-2.0d1*x(20)
      return
  349 fc=7.0d0*x(1)**2+5.0d0*x(2)**2+x(19)**2-3.0d1*x(20)
      return
      end subroutine tcfu07

! subroutine tcgu07             all systems                 90/12/01
! portability : all systems
! 90/12/01 lu : original version
!
! purpose :
!  gradients of test functions for nonlinear programming.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ii  kc  index of the constraint.
!  ri  x(n)  vector of variables.
!  ro  gc gradient of the constraint function.
!  ii  next  number of the test problem.
!
      subroutine tcgu07(n,kc,x,gc,next)
      integer n,kc,next
      double precision x(n),gc(n)
      double precision x1,x2,x3,x4,x5,x6,x7,x8
      integer i,j,k
      double precision y(128)
      common /empr07/ y
      go to(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170, &
       180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330, &
       340),next
   10 gc(1)=2.0d0*x(1)
      gc(2)=2.0d0*x(2)
      return
   20 gc(1)=6.0d0*x(1)-2.0d0*x(2)
      gc(2)=2.0d0*(x(2)-x(1))
      return
   30 gc(1)=8.0d0*x(1)
      gc(2)=2.0d0*x(2)
      return
   40 gc(1)=-3.0d0*(1.0d0-x(1))**2
      gc(2)=-1.0d0
      return
   50 gc(1)=-x(2)
      gc(2)=0.49d0-x(1)
      return
   60 go to (61,62),kc
   61 gc(1)=-1.0d0
      gc(2)= 2.0d0*x(2)
      return
   62 gc(1)=2.0d0*x(1)
      gc(2)=-1.0d0
      return
   70 go to (71,72),kc
   71 gc(1)=2.0d0*(x(1)-5.0d0)
      gc(2)=2.0d0*(x(2)-5.0d0)
      return
   72 gc(1)=2.0d0*(x(1)-6.0d0)
      gc(2)=2.0d0*(x(2)-5.0d0)
      return
   80 go to (81,82,10),kc
   81 gc(1)=1.0d0
      gc(2)=2.0d0*x(2)
      return
   82 gc(1)=2.0d0*x(1)
      gc(2)=1.0d0
      return
   90 go to (91,92,93),kc
   91 gc(1)=x(2)
      gc(2)=x(1)
      return
   92 gc(1)=-2.0d0*x(1)/125.0d0
      gc(2)= 1.0d0
      return
   93 gc(1)=-5.0d0
      gc(2)= 2.0d0*(x(2)-5.0d1)
      return
  100 gc(1)= 2.0d0*x(1)
      gc(2)= 4.0d0*x(2)
      gc(3)= 8.0d0*x(3)
      return
  110 gc(1)=2.0d0*x(1)
      gc(2)=2.0d0*x(2)
      gc(3)=0.0d0
      return
  120 gc(1)=x(2)
      gc(2)=x(1)
      gc(3)=0.0d0
      return
  130 gc(1)=-4.0d0/x(1)**2
      gc(2)=-3.2d1/x(2)**2
      gc(3)=-1.2d2/x(3)**2
      return
  140 go to (141,142),kc
  141 gc(1)=exp(x(1))
      gc(2)=-1.0d0
      gc(3)= 0.0d0
      return
  142 gc(1)= 0.0d0
      gc(2)=exp(x(2))
      gc(3)=-1.0d0
      return
  150 x2=1.6d0*x(1)
      y(4)=1.6d0
      y(5)=0.0d0
      y(6)=0.0d0
  151 x3=1.22d0*x2-x(1)
      y(7)=1.22d0*y(4)-1.0d0
      y(8)=1.22d0*y(5)
      y(9)=1.22d0*y(6)
      x6=(x(2)+x3)/x(1)
      y(16)=y(7)/x(1)-(x(2)+x3)/x(1)**2
      y(17)=(1.0d0+y(8))/x(1)
      y(18)=y(9)/x(1)
      x1=0.01d0*x(1)*(112.0d0+13.167d0*x6-0.6667d0*x6**2)
      y(1)=0.01d0*(112.0d0+13.167d0*x6-0.6667d0*x6**2)+ &
       0.01d0*x(1)*(13.167*y(16)-1.3334d0*x6*y(16))
      y(2)=0.01d0*x(1)*(13.167d0*y(17)-1.3334d0*x6*y(17))
      y(3)=0.01d0*x(1)*(13.167d0*y(18)-1.3334d0*x6*y(18))
      if (abs(x1-x2).gt.0.001d0) then
      x2=x1
      y(4)=y(1)
      y(5)=y(2)
      y(6)=y(3)
      go to 151
      endif
      x4=93.0d0
      y(10)=0.0d0
      y(11)=0.0d0
      y(12)=0.0d0
  152 x5=86.35d0+1.098d0*x6-0.038d0*x6**2+0.325d0*(x4-89.0d0)
      y(13)=1.098d0*y(16)-0.076d0*x6*y(16)+0.325d0*y(10)
      y(14)=1.098d0*y(17)-0.076d0*x6*y(17)+0.325d0*y(11)
      y(15)=1.098d0*y(18)-0.076d0*x6*y(18)+0.325d0*y(12)
      x8=3.0d0*x5-133.0d0
      y(22)=3.0d0*y(13)
      y(23)=3.0d0*y(14)
      y(24)=3.0d0*y(15)
      x7=35.82d0-0.222d0*x8
      y(19)=-0.222d0*y(22)
      y(20)=-0.222d0*y(23)
      y(21)=-0.222d0*y(24)
      y(89)=x2*x7+1000.0d0*x(3)
      y(90)=y(89)**2
      x1=98000.0d0*x(3)/y(89)
      y(1)=-98000.0d0*x(3)*(y(4)*x7+x2*y(19))/y(90)
      y(2)=-98000.0d0*x(3)*(y(5)*x7+x2*y(20))/y(90)
      y(3)=-98000.0d0*x(3)*(y(6)*x7+x2*y(21)+1000.0d0)/y(90)+ &
       98000.0d0/y(89)
      if (abs(x1-x4).gt.0.001d0) then
      x4=x1
      y(10)=y(1)
      y(11)=y(2)
      y(12)=y(3)
      go to 152
      endif
      go to (153,154,155,156,157,158,159),kc
  153 continue
      gc(1)=y(4)
      gc(2)=y(5)
      gc(3)=y(6)
      return
  154 continue
      gc(1)=y(7)
      gc(2)=y(8)
      gc(3)=y(9)
      return
  155 continue
      gc(1)=y(10)
      gc(2)=y(11)
      gc(3)=y(12)
      return
  156 continue
      gc(1)=y(13)
      gc(2)=y(14)
      gc(3)=y(15)
      return
  157 continue
      gc(1)=y(16)
      gc(2)=y(17)
      gc(3)=y(18)
      return
  158 continue
      gc(1)=y(19)
      gc(2)=y(20)
      gc(3)=y(21)
      return
  159 continue
      gc(1)=y(22)
      gc(2)=y(23)
      gc(3)=y(24)
      return
  160 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=1.0d0-x(4)
      gc(4)=1.0d0-x(3)
      return
  170 go to (171,172),kc
  171 gc(1)=-4.00d0/x(1)**2
      gc(2)=-2.25d0/x(2)**2
      gc(3)=-1.00d0/x(3)**2
      gc(4)=-0.25d0/x(4)**2
      return
  172 gc(1)=-0.16d0/x(1)**2
      gc(2)=-0.36d0/x(2)**2
      gc(3)=-0.64d0/x(3)**2
      gc(4)=-0.64d0/x(4)**2
      return
  180 x5=x(1)+x(1)
      x6=x(2)+x(2)
      x7=x(3)+x(3)
      x8=x(4)+x(4)
  181 go to (182,183,184),kc
  182 gc(1)=x5+1.0d0
      gc(2)=x6-1.0d0
      gc(3)=x7+1.0d0
      gc(4)=x8-1.0d0
      return
  183 gc(1)=x5-1.0d0
      gc(2)=x6+x6
      gc(3)=x7
      gc(4)=x8+x8-1.0d0
      return
  184 gc(1)=x5+x5+2.0d0
      gc(2)=x6-1.0d0
      gc(3)=x7
      gc(4)=-1.0d0
      return
  190 go to (191,192,193),kc
  191 gc(1)=2.0d0*x(1)
      gc(2)=2.0d0*x(2)
      gc(3)=2.0d0*x(3)
      gc(4)=2.0d0*x(4)
      gc(5)=2.0d0*x(5)
      return
  192 gc(1)=2.0d0*x(1)*x(3)
      gc(2)=0.0d0
      gc(3)=x(1)**2
      gc(4)=x(5)
      gc(5)=x(4)
      return
  193 gc(1)=10.0d0*x(5)
      gc(2)=2.0d0*x(2)*x(4)
      gc(3)=0.0d0
      gc(4)=x(2)**2
      gc(5)=10.0d0*x(1)
      return
  200 go to (201,202,203),kc
  201 gc(1)= 0.0006262d0*x(4)
      gc(2)= 0.0056868d0*x(5)
      gc(3)=-0.0022053d0*x(5)
      gc(4)= 0.0006262d0*x(1)
      gc(5)= 0.0056868d0*x(2)-0.0022053d0*x(3)
      return
  202 gc(1)=0.0029955d0*x(2)
      gc(2)=0.0071317d0*x(5)+0.0029955d0*x(1)
      gc(3)=0.0043626d0*x(3)
      gc(4)=0.0d0
      gc(5)=0.0071317d0*x(2)
      return
  203 gc(1)=0.0012547d0*x(3)
      gc(2)=0.0d0
      gc(3)=0.0047026d0*x(5)+0.0012547d0*x(1)+0.0019085d0*x(4)
      gc(4)=0.0019085d0*x(3)
      gc(5)=0.0047026d0*x(3)
      return
  210 go to (211,212,213),kc
  211 gc(1)=-145421.4020d0+2931.15060d0*x(2)- &
       40.4279320d0*x(3)+5106.1920d0*x(4)+15711.36d0*x(5)
      gc(2)=2931.15060d0*x(1)
      gc(3)=-40.4279320d0*x(1)
      gc(4)=5106.1920d0*x(1)
      gc(5)=15711.36d0*x(1)
      return
  212 gc(1)=-155011.1084d0+4360.53352d0*x(2)+ &
       12.9492344d0*x(3)+10236.884d0*x(4)+13176.786d0*x(5)
      gc(2)=4360.53352d0*x(1)
      gc(3)=12.9492344d0*x(1)
      gc(4)=10236.884d0*x(1)
      gc(5)=13176.786d0*x(1)
      return
  213 gc(1)=-326669.5104d0+7390.68412d0*x(2)- &
       27.8986976d0*x(3)+16643.076d0*x(4)+30988.146d0*x(5)
      gc(2)=7390.68412d0*x(1)
      gc(3)=-27.8986976d0*x(1)
      gc(4)=16643.076d0*x(1)
      gc(5)=30988.146d0*x(1)
      return
  220 if(kc.eq.2) then
      gc(1)=0.0d0
      gc(2)=1.5d0
      gc(3)=-1.0d0
      gc(4)=0.0d0
      gc(5)=0.0d0
      return
      endif
      x1=146.312d3
      y(1)=x(2)+x(3)+41.6d0
      y(18)=0.024d0*x(4)-4.62d0
      y(2)=12.5d0/y(18)+12.0d0
      y(19)=0.0003535d0*x(1)**2+0.5311d0*x(1)+0.08705d0*y(2)*x(1)
      y(20)=0.052d0*x(1)+78.0d0+0.002377d0*y(2)*x(1)
      y(3)=y(19)/y(20)
      y(4)=19.0d0*y(3)
      y(21)=0.04782d0*(x(1)-y(3))+0.1956d0*(x(1)-y(3))**2/x(2)+ &
       0.6376d0*y(4)+1.594d0*y(3)
      y(22)=100.0d0*x(2)
      y(23)=x(1)-y(3)-y(4)
      y(24)=0.95d0-y(21)/y(22)
      y(5)=y(23)*y(24)
      y(6)=x(1)-y(5)-y(4)-y(3)
      y(25)=(y(5)+y(4))*0.995d0
      y(7)=y(25)/y(1)
      y(8)=y(25)/3798.0d0
      y(26)=y(7)-0.0663d0*y(7)/y(8)-0.3153d0
      y(9)=96.82d0/y(26)+0.321d0*y(1)
      y(10)=1.29d0*y(5)+1.258d0*y(4)+2.29d0*y(3)+1.71d0*y(6)
      y(11)=1.71d0*x(1)-0.452d0*y(4)+0.58d0*y(3)
      y(27)=12.3d0/752.3d0
      y(28)=1.75d0*y(2)*0.995d0*x(1)
      y(29)=0.995d0*y(10)+1998.0d0
      y(12)=y(27)*x(1)+y(28)/y(29)
      y(13)=y(29)-1.75d0*y(2)
      y(14)=3623.0d0+64.4d0*x(2)+58.4d0*x(3)+x1/(y(9)+ &
       x(5))
      y(30)=0.995d0*y(10)+60.8d0*x(2)+48.0d0*x(4)-0.1121d0*y(14)- &
       5095.0d0
      y(15)=y(13)/y(30)
      y(16)=148.0d3-331.0d3*y(15)+40.0d0*y(13)-61.0d0*y(15)*y(13)
      y(31)=2324.0d0*y(10)-2.874d7*y(2)
      y(17)=1.413d7-1328.0d0*y(10)-531.0d0*y(11)+y(31)/y(29)
      y(32)=y(13)/y(15)-y(13)/0.52d0
      y(33)=1.104d0-0.72d0*y(15)
      y(34)=y(9)+x(5)
      y(35)=1.0d0
      y(36)=0.024d0
      y(37)=-12.5d0*y(36)/y(18)**2
      y(38)=0.000707d0*x(1)+0.5311d0+0.08705d0*y(2)
      y(39)=0.08705d0*x(1)*y(37)
      y(40)=0.052d0+0.002377d0*y(2)
      y(41)=0.002377d0*x(1)*y(37)
      y(42)=(y(38)*y(20)-y(19)*y(40))/y(20)**2
      y(43)=(y(39)*y(20)-y(19)*y(41))/y(20)**2
      y(44)=19.0d0*y(42)
      y(45)=19.0d0*y(43)
      y(46)=0.04782d0*(1.0d0-y(42))+0.3912d0*(x(1)-y(3))*(1.0d0- &
       y(42))/x(2)+0.6376d0*y(44)+1.594d0*y(42)
      y(47)=-0.1956d0*(x(1)-y(3))**2/x(2)**2
      y(48)=-0.04782d0*y(43)-0.3912d0*(x(1)-y(3))*y(43)/x(2)+ &
       0.6376d0*y(45)+1.594d0*y(43)
      y(49)=100.0d0
      y(50)=1.0d0-y(42)-y(44)
      y(51)=-y(43)-y(45)
      y(52)=-y(46)/y(22)
      y(53)=(y(21)*y(49)-y(47)*y(22))/y(22)**2
      y(54)=-y(48)/y(22)
      y(55)=y(50)*y(24)+y(23)*y(52)
      y(56)=y(23)*y(53)
      y(57)=y(51)*y(24)+y(23)*y(54)
      y(58)=1.0d0-y(55)-y(44)-y(42)
      y(59)=-y(56)
      y(60)=-y(57)-y(45)-y(43)
      y(61)=0.995d0*(y(55)+y(44))
      y(62)=0.995d0*y(56)
      y(63)=0.995d0*(y(57)+y(45))
      y(64)=y(61)/y(1)
      y(65)=(y(62)*y(1)-y(25)*y(35))/y(1)**2
      y(66)=-y(25)*y(35)/y(1)**2
      y(67)=y(63)/y(1)
      y(68)=y(61)/3798.0d0
      y(69)=y(62)/3798.0d0
      y(70)=y(63)/3798.0d0
      y(71)=y(64)-0.0663d0*(y(64)*y(8)-y(7)*y(68))/y(8)**2
      y(72)=y(65)-0.0663d0*(y(65)*y(8)-y(7)*y(69))/y(8)**2
      y(73)=y(66)-0.0663d0*y(66)/y(8)
      y(74)=y(67)-0.0663d0*(y(67)*y(8)-y(7)*y(70))/y(8)**2
      y(75)=-96.82d0*y(71)/y(26)**2
      y(76)=-96.82d0*y(72)/y(26)**2+0.321d0*y(35)
      y(77)=-96.82d0*y(73)/y(26)**2+0.321d0*y(35)
      y(78)=-96.82d0*y(74)/y(26)**2
      y(79)=1.29d0*y(55)+1.258d0*y(44)+2.29d0*y(42)+1.71d0*y(58)
      y(80)=1.29d0*y(56)+1.71d0*y(59)
      y(81)=1.29d0*y(57)+1.258d0*y(45)+2.29d0*y(43)+1.71d0*y(60)
      y(82)=1.71d0-0.452d0*y(44)+0.58d0*y(42)
      y(83)=-0.452d0*y(45)+0.58d0*y(43)
      y(84)=1.75d0*y(2)*0.995d0
      y(85)=1.75d0*y(37)*0.995d0*x(1)
      y(86)=0.995d0*y(79)
      y(87)=0.995d0*y(80)
      y(88)=0.995d0*y(81)
      y(89)=y(27)+(y(84)*y(29)-y(28)*y(86))/y(29)**2
      y(90)=-y(28)*y(87)/y(29)**2
      y(91)=(y(85)*y(29)-y(28)*y(88))/y(29)**2
      y(92)=y(88)-1.75d0*y(37)
      y(93)=-x1*y(75)/(y(9)+x(5))**2
      y(94)=64.4d0-x1*y(76)/(y(9)+x(5))**2
      y(95)=58.4d0-x1*y(77)/(y(9)+x(5))**2
      y(96)=-x1*y(78)/(y(9)+x(5))**2
      y(97)=-x1/(y(9)+x(5))**2
      y(98)=0.995d0*y(79)-0.1121d0*y(93)
      y(99)=0.995d0*y(80)+60.8d0-0.1121d0*y(94)
      y(100)=-0.1121d0*y(95)
      y(101)=0.995d0*y(81)+48.0d0-0.1121d0*y(96)
      y(102)=-0.1121d0*y(97)
      y(103)=(y(86)*y(30)-y(13)*y(98))/y(30)**2
      y(104)=(y(87)*y(30)-y(13)*y(99))/y(30)**2
      y(105)=-y(13)*y(100)/y(30)**2
      y(106)=(y(92)*y(30)-y(13)*y(101))/y(30)**2
      y(107)=-y(13)*y(102)/y(30)**2
      y(108)=-3.31d5*y(103)+40.0d0*y(86)-61.0d0*(y(103)*y(13)+ &
       y(15)*y(86))
      y(109)=-3.31d5*y(104)+40.0d0*y(87)-61.0d0*(y(104)*y(13)+ &
       y(15)*y(87))
      y(110)=-3.31d5*y(105)-61.0d0*y(105)*y(13)
      y(111)=-3.31d5*y(106)+40.0d0*y(92)-61.0d0*(y(106)*y(13)+ &
       y(15)*y(92))
      y(112)=-3.31d5*y(107)-61.0d0*y(107)*y(13)
      y(113)=2.324d3*y(79)
      y(114)=2.324d3*y(80)
      y(115)=2.324d3*y(81)-2.874d7*y(37)
      y(116)=-1.328d3*y(79)-531.0d0*y(82)+(y(113)*y(29)-y(31)* &
       y(86))/y(29)**2
      y(117)=-1.328d3*y(80)+(y(114)*y(29)-y(31)*y(87))/y(29)**2
      y(118)=-1.328d3*y(81)-531.0d0*y(83)+(y(115)*y(29)-y(31)* &
       y(88))/y(29)**2
      y(119)=(y(86)*y(15)-y(13)*y(103))/y(15)**2-y(86)/0.52d0
      y(120)=(y(87)*y(15)-y(13)*y(104))/y(15)**2-y(87)/0.52d0
      y(121)=-y(13)*y(105)/y(15)**2
      y(122)=(y(92)*y(15)-y(13)*y(106))/y(15)**2-y(92)/0.52d0
      y(123)=-y(13)*y(107)/y(15)**2
      y(124)=-0.72d0*y(103)
      y(125)=-0.72d0*y(104)
      y(126)=-0.72d0*y(105)
      y(127)=-0.72d0*y(106)
      y(128)=-0.72d0*y(107)
      if (kc.eq.1) then
      gc(1)=y(44)-0.28d0*y(55)/0.72d0
      gc(2)=-0.28d0*y(56)/0.72d0
      gc(3)=0.0d0
      gc(4)=y(45)-0.28d0*y(57)/0.72d0
      gc(5)=0.0d0
      return
      endif
      go to (503,504,505,506,507,508,509,510,511,512,513,514,515,516, &
       517,518,519,520,521),kc-2
  503 gc(1)=-3496.0d0*y(2)*y(86)/y(29)**2
      gc(2)=-3496.0d0*y(2)*y(87)/y(29)**2
      gc(3)=0.0d0
      gc(4)=3496.0d0*(y(37)*y(29)-y(2)*y(88))/y(29)**2
      gc(5)=0.0d0
      return
  504 gc(1)=-62212.0d0*y(75)/y(34)**2
      gc(2)=-62212.0d0*y(76)/y(34)**2-y(35)
      gc(3)=-62212.0d0*y(77)/y(34)**2-y(35)
      gc(4)=-62212.0d0*y(78)/y(34)**2
      gc(5)=-62212.0d0*y(35)/y(34)**2
      return
  505 gc(1)=0.0d0
      gc(2)=y(35)
      gc(3)=y(35)
      gc(4)=0.0d0
      gc(5)=0.0d0
      return
  506 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=y(37)
      gc(5)=0.0d0
      return
  507 gc(1)=y(42)
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=y(43)
      gc(5)=0.0d0
      return
  508 gc(1)=y(44)
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=y(45)
      gc(5)=0.0d0
      return
  509 gc(1)=y(55)
      gc(2)=y(56)
      gc(3)=0.0d0
      gc(4)=y(57)
      gc(5)=0.0d0
      return
  510 gc(1)=y(58)
      gc(2)=y(59)
      gc(3)=0.0d0
      gc(4)=y(60)
      gc(5)=0.0d0
      return
  511 gc(1)=y(64)
      gc(2)=y(65)
      gc(3)=y(66)
      gc(4)=y(67)
      gc(5)=0.0d0
      return
  512 gc(1)=y(68)
      gc(2)=y(69)
      gc(3)=0.0d0
      gc(4)=y(70)
      gc(5)=0.0d0
      return
  513 gc(1)=y(75)
      gc(2)=y(76)
      gc(3)=y(77)
      gc(4)=y(78)
      gc(5)=0.0d0
      return
  514 gc(1)=y(79)
      gc(2)=y(80)
      gc(3)=0.0d0
      gc(4)=y(81)
      gc(5)=0.0d0
      return
  515 gc(1)=y(82)
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=y(83)
      gc(5)=0.0d0
      return
  516 gc(1)=y(89)
      gc(2)=y(90)
      gc(3)=0.0d0
      gc(4)=y(91)
      gc(5)=0.0d0
      return
  517 gc(1)=y(86)
      gc(2)=y(87)
      gc(3)=0.0d0
      gc(4)=y(92)
      gc(5)=0.0d0
      return
  518 gc(1)=y(93)
      gc(2)=y(94)
      gc(3)=y(95)
      gc(4)=y(96)
      gc(5)=y(97)
      return
  519 gc(1)=y(103)
      gc(2)=y(104)
      gc(3)=y(105)
      gc(4)=y(106)
      gc(5)=y(107)
      return
  520 gc(1)=y(108)
      gc(2)=y(109)
      gc(3)=y(110)
      gc(4)=y(111)
      gc(5)=y(112)
      return
  521 gc(1)=y(116)
      gc(2)=y(117)
      gc(3)=0.0d0
      gc(4)=y(118)
      gc(5)=0.0d0
      return
  230 go to (231,232),kc
  231 gc(1)=2.0d0*x(1)
      gc(2)=2.0d0*x(2)
      gc(3)=2.0d0*x(3)
      gc(4)=0.0d0
      gc(5)=0.0d0
      gc(6)=0.0d0
      return
  232 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=2.0d0*(x(4)-3.0d0)
      gc(5)=2.0d0*x(5)
      gc(6)=0.0d0
      return
  240 go to (241,242),kc
  241 gc(1)=0.001d0*x(2)*x(3)*x(4)*x(5)*x(6)
      gc(2)=0.001d0*x(1)*x(3)*x(4)*x(5)*x(6)
      gc(3)=0.001d0*x(1)*x(2)*x(4)*x(5)*x(6)
      gc(4)=0.001d0*x(1)*x(2)*x(3)*x(5)*x(6)
      gc(5)=0.001d0*x(1)*x(2)*x(3)*x(4)*x(6)
      gc(6)=0.001d0*x(1)*x(2)*x(3)*x(4)*x(5)
      return
  242 gc(1)=0.00124d0*x(1)*x(4)*x(5)**2+ &
       0.00062d0*x(4)*x(5)**2*(x(2)+x(3))+ &
       0.00058d0*x(2)*x(3)*x(6)**2
      gc(2)=0.00062d0*x(1)*x(4)*x(5)**2+ &
       0.00116d0*1.57d0*x(2)*x(3)*x(6)**2+ &
       0.00058d0*x(3)*x(6)**2*(x(1)+x(4))
      gc(3)=0.00062d0*x(1)*x(4)*x(5)**2+ &
       0.00058d0*x(2)*x(6)**2*(x(1)+1.57d0*x(2)+x(4))
      gc(4)=0.00062d0*x(1)*x(5)**2*(x(1)+x(2)+x(3))+ &
       0.00058d0*x(2)*x(3)*x(6)**2
      gc(5)=0.00124d0*x(1)*x(4)*x(5)*(x(1)+x(2)+x(3))
      gc(6)=0.00116d0*x(2)*x(3)*x(6)*(x(1)+1.57d0*x(2)+x(4))
      return
  250 go to (251,252,253,254),kc
  251 gc(1)=17.1d0-169.0d0*x(3)
      gc(2)=38.2d0
      gc(3)=204.2d0-169.0d0*x(1)-3580.0d0*x(5)
      gc(4)=212.3d0-3810.0d0*x(5)-18500.0d0*x(6)
      gc(5)=623.4d0-3580.0d0*x(3)-3810.0d0*x(4)-24300.0d0*x(6)
      gc(6)=1495.5d0-18500.0d0*x(4)-24300.0d0*x(5)
      return
  252 gc(1)=17.9d0-139.0d0*x(3)
      gc(2)=36.8d0
      gc(3)=113.9d0-139.0d0*x(1)
      gc(4)=169.7d0-2450.0d0*x(5)-16600.0d0*x(6)
      gc(5)=337.8d0-2450.0d0*x(4)-17200.0d0*x(6)
      gc(6)=1385.2d0-16600.0d0*x(4)-17200.0d0*x(5)
      return
  253 gc(1)=0.0d0
      gc(2)=-273.0d0
      gc(3)=0.0d0
      gc(4)=-70.0d0+26000.0d0*x(5)
      gc(5)=-819.0d0+26000.0d0*x(4)
      gc(6)=0d0
      return
  254 gc(1)=159.9d0-14000.0d0*x(6)
      gc(2)=-311.0d0
      gc(3)=0.0d0
      gc(4)=587.0d0
      gc(5)=391.0d0
      gc(6)=2198.0d0-14000.0d0*x(1)
      return
  260 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=0.0d0
      gc(5)=0.0d0
      gc(6)=0.0d0
      gc(7)=0.0d0
      go to (261,262,263,264),kc
  261 gc(1)=4.0d0*x(1)
      gc(2)=1.2d1*x(2)**3
      gc(3)=1.0d0
      gc(4)=8.0d0*x(4)
      gc(5)=5.0d0
      return
  262 gc(1)=7.0d0
      gc(2)=3.0d0
      gc(3)=2.0d1*x(3)
      gc(4)=1.0d0
      gc(5)=-1.0d0
      return
  263 gc(1)=2.3d1
      gc(2)=2.0d0*x(2)
      gc(6)=1.2d1*x(6)
      gc(7)=-8.0d0
      return
  264 gc(1)=8.0d0*x(1)-3.0d0*x(2)
      gc(2)=2.0d0*x(2)-3.0d0*x(1)
      gc(3)=4.0d0*x(3)
      gc(6)=5.0d0
      gc(7)=-1.1d1
      return
  270 go to (271,272,273,274,275),kc
  271 gc(1)=0.25d0*x(7)/(x(1)**0.5d0*x(3)*x(6)**2)+ &
       2.1d0*x(1)**2*x(2)*x(6)*x(7)**0.5d0/x(3)**2
      gc(2)=0.7d0*x(1)**3*x(6)*x(7)**0.5d0/x(3)**2- &
       0.2d0*x(3)*x(6)**(2.0d0/3.0d0)*x(7)**(1.0d0/4.0d0)/ &
       (x(2)**2*x(4)**0.5d0)
      gc(3)=-0.5d0*x(1)**0.5d0*x(7)/(x(3)*x(6))**2- &
       1.4d0*x(1)**3*x(2)*x(6)*x(7)**0.5d0/x(3)**3+ &
       0.2d0*x(6)**(2.0d0/3.0d0)*x(7)**(1.0d0/4.0d0)/ &
       (x(2)*x(4)**0.5d0)
      gc(4)=-0.1d0*x(3)*x(6)**(2.0d0/3.0d0)* &
       x(7)**(1.0d0/4.0d0)/(x(2)*x(4)**1.5d0)
      gc(5)=0.0d0
      gc(6)=-1.0d0*x(1)**0.5d0*x(7)/(x(3)*x(6)**3)+ &
       0.7d0*x(1)**3*x(2)*x(7)**0.5d0/x(3)**2+ &
       (0.4d0/3.0d0)*x(3)*x(7)**(1.0d0/4.0d0)/ &
       (x(2)*x(4)**0.5d0*x(6)**(1.0d0/3.0d0))
      gc(7)=0.5d0*x(1)**0.5d0/(x(3)*x(6)**2)+ &
       0.35d0*x(1)**3*x(2)*x(6)/(x(3)**2*x(7)**0.5d0)+ &
       0.05d0*x(3)*x(6)**(2.0d0/3.0d0)/ &
       (x(2)*x(4)**0.5d0*x(7)**0.75d0)
      return
  272 gc(1)=-0.65d0*x(2)*x(6)/(x(1)**1.5d0*x(3)*x(5))- &
       3.1d0*x(2)**0.5d0*x(6)**(1.0d0/3.0d0)/ &
       ((x(1)*x(4))**2*x(5))
      gc(2)=1.3d0*x(6)/(x(1)**0.5d0*x(3)*x(5))+ &
       1.55d0*x(6)**(1.0d0/3.0d0)/ &
       (x(1)*x(2)**0.5d0*x(4)**2*x(5))
      gc(3)=-1.3d0*x(2)*x(6)/(x(1)**0.5d0*x(3)**2*x(5))+ &
       0.8d0*x(6)**2/(x(4)*x(5))
      gc(4)=-0.8d0*x(3)*x(6)**2/(x(4)**2*x(5))- &
       6.2d0*x(2)**0.5d0*x(6)**(1.0d0/3.0d0)/ &
       (x(1)*x(4)**3*x(5))
      gc(5)=-1.3d0*x(2)*x(6)/(x(1)**0.5d0*x(3)*x(5)**2)- &
       0.8d0*x(3)*x(6)**2/(x(4)*x(5)**2)- &
       3.1d0*x(2)**0.5d0*x(6)**(1.0d0/3.0d0)/ &
       (x(1)*(x(4)*x(5))**2)
      gc(6)=1.3d0*x(2)/(x(1)**0.5d0*x(3)*x(5))+ &
       1.6d0*x(3)*x(6)/(x(4)*x(5))+ &
       (3.1d0/3.0d0)*x(2)**0.5d0/ &
       (x(1)*x(4)**2*x(5)*x(6)**(2.0d0/3.0d0))
      gc(7)=0.0d0
      return
  273 gc(1)=2.0d0*x(5)*x(7)**(4.0d0/3.0d0)/(x(3)**1.5d0*x(6))- &
       x(2)*x(3)**0.5d0*x(5)/x(1)**2
      gc(2)=0.1d0*x(5)/(x(3)**0.5d0*x(6)*x(7)**0.5d0)+ &
       x(3)**0.5d0*x(5)/x(1)- &
       1.3d0*x(3)*x(5)*x(7)/(x(2)**3*x(6))
      gc(3)=-3.0d0*x(1)*x(5)*x(7)**(4.0d0/3.0d0)/(x(3)**2.5d0* &
       x(6))-0.05d0*x(2)*x(5)/(x(3)**1.5d0*x(6)*x(7)**0.5d0)+ &
       0.5d0*x(2)*x(5)/(x(1)*x(3)**0.5d0)+ &
       0.65d0*x(5)*x(7)/(x(2)**2*x(6))
      gc(4)=0.0d0
      gc(5)=2.0d0*x(1)*x(7)**(4.0d0/3.0d0)/(x(3)**1.5d0*x(6))+ &
       0.1d0*x(2)/(x(3)**0.5d0*x(6)*x(7)**0.5d0)+ &
       x(2)*x(3)**0.5d0/x(1)+ &
       0.65d0*x(3)*x(7)/(x(2)**2*x(6))
      gc(6)=-2.0d0*x(1)*x(5)*x(7)**(4.0d0/3.0d0)/(x(3)**1.5d0* &
       x(6)**2)-0.1d0*x(2)*x(5)/(x(3)**0.5d0*x(6)**2*x(7)**0.5d0)- &
       0.65d0*x(3)*x(5)*x(7)/(x(2)**2*x(6)**2)
      gc(7)=(8.0d0/3.0d0)*x(1)*x(5)*x(7)**(1.0d0/3.0d0)/ &
       (x(3)**1.5d0*x(6))-0.05d0*x(2)*x(5)/(x(3)**0.5d0*x(6)* &
       x(7)**1.5d0)+0.65d0*x(3)*x(5)/(x(2)**2*x(6))
      return
  274 gc(1)=-0.4d0*x(2)*x(5)**0.5d0*x(7)**(1.0d0/3.0d0)/ &
       (x(1)**3*x(4))+ &
       0.15d0*x(2)**2*x(3)*x(4)**(1.0d0/3.0d0)* &
       x(7)**(1.0d0/4.0d0)/(x(1)**0.5d0*x(5)**(2.0d0/3.0d0))- &
       1.2d0*x(3)*x(5)*x(7)**(3.0d0/4.0d0)/(x(1)**4*x(2)**2)
      gc(2)=0.2d0*x(5)**0.5d0*x(7)**(1.0d0/3.0d0)/ &
       (x(1)**2*x(4))+ &
       0.6d0*x(1)**0.5d0*x(2)*x(3)*x(4)**(1.0d0/3.0d0)* &
       x(7)**0.25d0/x(5)**(2.0d0/3.0d0)- &
       0.8d0*x(3)*x(5)*x(7)**0.75d0/(x(1)*x(2))**3
      gc(3)=0.3d0*x(1)**0.5d0*x(2)**2*x(4)**(1.0d0/3.0d0)* &
       x(7)**(1.0d0/4.0d0)/x(5)**(2.0d0/3.0d0)+ &
       0.4d0*x(5)*x(7)**(3.0d0/4.0d0)/(x(1)**3*x(2)**2)- &
       x(4)*x(7)**0.5d0/x(3)**3
      gc(4)=-0.2d0*x(2)*x(5)**0.5d0*x(7)**(1.0d0/3.0d0)/ &
       (x(1)*x(4))**2+ &
       (0.3d0/3.0d0)*x(1)**0.5d0*x(2)**2*x(3)* &
       x(7)**(1.0d0/4.0d0)/(x(4)**(2.0d0/3.0d0)*x(5)** &
       (2.0d0/3.0d0))+0.5d0*x(7)**0.5d0/x(3)**2
      gc(5)=0.1d0*x(2)*x(7)**(1.0d0/3.0d0)/ &
       (x(1)**2*x(4)*x(5)**0.5d0)- &
       (0.6d0/3.0d0)*x(1)**0.5d0*x(2)**2*x(3)*x(4)**(1.0d0/ &
       3.0d0)*x(7)**(1.0d0/4.0d0)/x(5)**(5.0d0/3.0d0)+ &
       0.4d0*x(3)*x(7)**(3.0d0/4.0d0)/(x(1)**3*x(2)**2)
      gc(6)=0.0d0
      gc(7)=(0.2d0/3.0d0)*x(2)*x(5)**0.5d0/ &
       (x(1)**2*x(4)*x(7)**(2.0d0/3.0d0))+ &
       0.075d0*x(1)**0.5d0*x(2)**2*x(3)*x(4)**(1.0d0/3.0d0)/ &
       (x(5)**(2.0d0/3.0d0)*x(7)**(3.0d0/4.0d0))+ &
       0.3d0*x(3)*x(5)/(x(1)**3*x(2)**2*x(7)**(1.0d0/4.0d0))+ &
       0.25d0*x(4)/(x(3)**2*x(7)**0.5d0)
      return
  275 gc(1)=10.0d0*x(4)**2/(x(2)*x(6)**3*x(7)**0.25d0)- &
       15.0d0*x(3)*x(4)/((x(1)*x(2))**2*x(5)*x(7)**0.5d0)- &
       40.0d0*x(2)*x(6)/(x(1)**3*x(4)*x(5)**2)+ &
       50.0d0*x(1)*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6)**2)
      gc(2)=-10.0d0*x(1)*x(4)**2/(x(2)**2*x(6)**3*x(7)**0.25d0)- &
       30.0d0*x(3)*x(4)/(x(1)*x(2)**3*x(5)*x(7)**0.5d0)+ &
       20.0d0*x(6)/(x(1)**2*x(4)*x(5)**2)+ &
       50.0d0*x(1)**2*x(2)*x(5)**0.5d0*x(7)/(x(3)*x(6)**2)
      gc(3)=15.0d0*x(4)/(x(1)*x(2)**2*x(5)*x(7)**0.5d0)- &
       25.0d0*x(1)**2*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6))**2
      gc(4)=20.0d0*x(1)*x(4)/(x(2)*x(6)**3*x(7)**0.25d0)+ &
       15.0d0*x(3)/(x(1)*x(2)**2*x(5)*x(7)**0.5d0)- &
       20.0d0*x(2)*x(6)/(x(1)*x(4)*x(5))**2
      gc(5)=-15.0d0*x(3)*x(4)/(x(1)*(x(2)*x(5))**2*x(7)**0.5d0)- &
       40.0d0*x(2)*x(6)/(x(1)**2*x(4)*x(5)**3)+ &
       12.5d0*x(1)**2*x(2)**2*x(7)/(x(3)*x(5)**0.5d0*x(6)**2)
      gc(6)=-30.0d0*x(1)*x(4)**2/(x(2)*x(6)**4*x(7)**0.25d0)+ &
       20.0d0*x(2)/(x(1)**2*x(4)*x(5)**2)- &
       50.0d0*x(1)**2*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6)**3)
      gc(7)=-2.5d0*x(1)*x(4)**2/(x(2)*x(6)**3*x(7)**1.25d0)- &
       7.5d0*x(3)*x(4)/(x(1)*x(2)**2*x(5)*x(7)**1.5d0)+ &
       25.0d0*x(1)**2*x(2)**2*x(5)**0.5d0/(x(3)*x(6)**2)
      return
  280 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=0.0d0
      gc(5)=0.0d0
      gc(6)=0.0d0
      gc(7)=0.0d0
      gc(8)=0.0d0
      go to (281,282,283,284,285),kc
  281 gc(1)=0.1d0
      gc(5)=0.0588d0*x(7)
      gc(7)=0.0588d0*x(5)
      return
  282 gc(1)=0.1d0
      gc(2)=0.1d0
      gc(6)=0.0588d0*x(8)
      gc(8)=0.0588d0*x(6)
      return
  283 gc(3)=4.0d0/x(5)-1.42d0/(x(3)**1.71d0*x(5))- &
       0.0588d0*1.3d0*x(7)/x(3)**2.3d0
      gc(5)=-4.0d0*x(3)/x(5)**2-2.0d0/(x(3)**0.71d0*x(5)**2)
      gc(7)=0.0588d0/x(3)**1.3d0
      return
  284 gc(4)=4.0d0/x(6)-1.42d0/(x(4)**1.71d0*x(6))- &
       0.0588d0*1.3d0*x(8)/x(4)**2.3d0
      gc(6)=-4.0d0*x(4)/x(6)**2-2.0d0/(x(4)**0.71d0*x(6)**2)
      gc(8)=0.0588d0/x(4)**1.3d0
      return
  285 gc(1)=0.268d0*(x(7)/x(1))**0.33d0/x(7)-1.0d0
      gc(2)=0.268d0*(x(8)/x(2))**0.33d0/x(8)-1.0d0
      gc(7)=-0.268d0*(x(1)/x(7))**0.67d0/x(7)
      gc(8)=-0.268d0*(x(2)/x(8))**0.67d0/x(8)
      return
  290 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=0.0d0
      gc(5)=0.0d0
      gc(6)=0.0d0
      gc(7)=0.0d0
      gc(8)=0.0d0
      go to (291,292,293,294,295,296),kc
  291 gc(4)=0.0025d0
      gc(6)=0.0025d0
      return
  292 gc(4)=-0.0025d0
      gc(5)= 0.0025d0
      gc(7)= 0.0025d0
      return
  293 gc(5)=-0.01d0
      gc(8)= 0.01d0
      return
  294 gc(1)=x(6)-100.0d0
      gc(4)=-833.33252d0
      gc(6)=x(1)
      return
  295 gc(2)=x(7)-x(4)
      gc(4)=-x(2)+1250.0d0
      gc(5)=-1250.0d0
      gc(7)=x(2)
      return
  296 gc(3)=x(8)-x(5)
      gc(5)=-x(3)+2500.0d0
      gc(8)=x(3)
      return
  300 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=0.0d0
      gc(5)=0.0d0
      gc(6)=0.0d0
      gc(7)=0.0d0
      gc(8)=0.0d0
      gc(9)=0.0d0
      go to (301,302,303,304,305,306,307,308,309,910,911,912,913),kc
  301 gc(3)=2.0d0*x(3)
      gc(4)=2.0d0*x(4)
      return
  302 gc(5)=2.0d0*x(5)
      gc(6)=2.0d0*x(6)
      return
  303 gc(1)= 2.0d0*(x(1)-x(5))
      gc(2)= 2.0d0*(x(2)-x(6))
      gc(5)=-2.0d0*(x(1)-x(5))
      gc(6)=-2.0d0*(x(2)-x(6))
      return
  304 gc(1)= 2.0d0*(x(1)-x(7))
      gc(2)= 2.0d0*(x(2)-x(8))
      gc(7)=-2.0d0*(x(1)-x(7))
      gc(8)=-2.0d0*(x(2)-x(8))
      return
  305 gc(3)= 2.0d0*(x(3)-x(5))
      gc(4)= 2.0d0*(x(4)-x(6))
      gc(5)=-2.0d0*(x(3)-x(5))
      gc(6)=-2.0d0*(x(4)-x(6))
      return
  306 gc(3)= 2.0d0*(x(3)-x(7))
      gc(4)= 2.0d0*(x(4)-x(8))
      gc(7)=-2.0d0*(x(3)-x(7))
      gc(8)=-2.0d0*(x(4)-x(8))
      return
  307 gc(7)= 2.0d0*x(7)
      gc(8)= 2.0d0*(x(8)-x(9))
      gc(9)=-2.0d0*(x(8)-x(9))
      return
  308 gc(1)= 2.0d0*x(1)
      gc(2)= 2.0d0*(x(2)-x(9))
      gc(9)=-2.0d0*(x(2)-x(9))
      return
  309 gc(9)=2.0d0*x(9)
      return
  910 gc(3)=x(9)
      gc(9)=x(3)
      return
  911 gc(5)= x(8)
      gc(6)=-x(7)
      gc(7)=-x(6)
      gc(8)= x(5)
      return
  912 gc(1)= x(4)
      gc(2)=-x(3)
      gc(3)=-x(2)
      gc(4)= x(1)
      return
  913 gc(5)=x(9)
      gc(9)=x(5)
      return
  310 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=0.0d0
      gc(5)=0.0d0
      gc(6)=0.0d0
      gc(7)=0.0d0
      gc(8)=0.0d0
      gc(9)=0.0d0
      gc(10)=0.0d0
      go to (311,312,313,314,315,316,317,318),kc
  311 gc(1)=6.0d0*(x(1)-2.0d0)
      gc(2)=8.0d0*(x(2)-3.0d0)
      gc(3)=4.0d0*x(3)
      gc(4)=-7.0d0
      return
  312 gc(1)=1.0d1*x(1)
      gc(2)=8.0d0
      gc(3)=2.0d0*(x(3)-6.0d0)
      gc(4)=-2.0d0
      return
  313 gc(1)=1.0d0*(x(1)-8.0d0)
      gc(2)=4.0d0*(x(2)-4.0d0)
      gc(5)=6.0d0*x(5)
      gc(6)=-1.0d0
      return
  314 gc(1)=2.0d0*x(1)-2.0d0*x(2)
      gc(2)=4.0d0*(x(2)-2.0d0)-2.0d0*x(1)
      gc(5)=1.4d1
      gc(6)=-6.0d0
      return
  315 gc(1)=4.0d0
      gc(2)=5.0d0
      gc(7)=-3.0d0
      gc(8)=9.0d0
      return
  316 gc(1)=1.0d1
      gc(2)=-8.0d0
      gc(7)=-1.7d1
      gc(8)=2.0d0
      return
  317 gc(1)=-3.0d0
      gc(2)=6.0d0
      gc(9)=2.4d1*(x(9)-8.0d0)
      gc(10)=-7.0d0
      return
  318 gc(1)=-8.0d0
      gc(2)=2.0d0
      gc(9)=5.0d0
      gc(10)=-2.0d0
      return
  320 k=(kc-1)*5
      do 321 j=1,5
      gc(j+10)=2.0d0*y(k+j+60)
  321 continue
      gc(kc+10)=gc(kc+10)+6.0d0*y(kc+85)*x(kc+10)
      k=(kc-1)*10
      do 322 i=1,10
      gc(i)=-y(k+i)
  322 continue
      return
  330 do 331 i=1,16
      gc(i)=0.0d0
  331 continue
      if (kc.lt.6) then
      k=11+3*(kc-1)
      gc(kc)=y(k+1)+(y(k)+2.0d0*x(kc)*y(k+2))/x(kc+5)
      gc(kc+5)=-(y(k)+y(k+2)*x(kc))*x(kc)/x(kc+5)**2
      return
      else if (kc.lt.14) then
      go to (332,333,334,335,336,337,338,339),kc-5
  332 gc(1)=y(27)*x(12)/(x(7)*x(11))
      gc(6)=(y(26)+y(28)*x(12)/x(11))/x(7)
      gc(7)=-(y(26)*x(6)+(y(27)*x(1)+y(28)*x(6))*x(12)/x(11))/x(7)**2
      gc(11)=-((y(27)*x(1)+y(28)*x(6))*x(12)/x(7))/x(11)**2
      gc(12)=(y(27)*x(1)+y(28)*x(6))/(x(7)*x(11))
      return
  333 gc(1)=y(33)*x(12)/x(8)
      gc(2)=y(31)*x(13)/x(8)
      gc(7)=(y(29)+y(30)*x(12))/x(8)
      gc(8)=-(y(29)*x(7)+(y(30)*x(7)+y(33)*x(1))*x(12)+ &
       y(31)*x(2)*x(13))/x(8)**2
      gc(12)=(y(30)*x(7)+y(33)*x(1))/x(8)
      gc(13)=y(32)+y(31)*x(2)/x(8)
      return
  334 gc(2)=y(38)*x(13)
      gc(3)=y(36)*x(14)
      gc(8)=y(35)*x(13)+y(34)
      gc(9)=y(39)*x(14)+y(37)
      gc(13)=y(35)*x(8)+y(38)*x(2)
      gc(14)=y(36)*x(3)+y(39)*x(9)
      return
  335 gc(3)=-((y(40)+y(43)/x(14))*x(9)+((y(41)*x(4)+y(44)*x(8))*x(15)+ &
       y(42)*x(10))/x(14))/x(3)**2
      gc(4)=y(41)*x(15)/(x(3)*x(14))
      gc(8)=y(44)*x(15)/(x(3)*x(14))
      gc(9)=(y(40)+y(43)/x(14))/x(3)
      gc(10)=y(42)/(x(3)*x(14))
      gc(14)=-((y(41)*x(4)+y(44)*x(8))*x(15)+y(42)*x(10)+ &
       y(43)*x(9))/(x(3)*x(14)**2)
      gc(15)=(y(41)*x(4)+y(44)*x(8))/(x(3)*x(14))
      return
  336 gc(4)=-((y(45)*x(5)*x(16)+y(49)*x(10))/x(15)+y(46)*x(10))/x(4)**2
      gc(5)=y(45)*x(16)/(x(4)*x(15))
      gc(10)=(y(46)+y(49)/x(15))/x(4)
      gc(15)=-((y(45)*x(5)*x(16)+y(49)*x(10))/x(4)+y(47)+ &
       y(48)*x(16))/x(15)**2
      gc(16)=(y(45)*x(5)/x(4)+y(48))/x(15)
      return
  337 gc(4)=-(y(50)+y(52)*x(5)*x(16))/x(4)**2
      gc(5)=y(52)*x(16)/x(4)
      gc(16)=y(51)+y(52)*x(5)/x(4)
      return
  338 gc(11)=y(53)
      gc(12)=y(54)
      return
  339 gc(11)=-y(55)*x(12)/x(11)**2
      gc(12)=y(55)/x(11)
      return
      else if (kc.lt.18) then
      i=kc-13
      k=5-i
      gc(k)=y(i+55)/x(k+1)
      gc(k+1)=-y(i+55)*x(k)/x(k+1)**2
      return
      else if (kc.eq.18) then
      gc(9)=y(60)/x(10)
      gc(10)=-y(60)*x(9)/x(10)**2
      return
      else
      gc(8)=y(61)/x(9)
      gc(9)=-y(61)*x(8)/x(9)**2
      return
      endif
  340 gc(1)=0.0d0
      gc(2)=0.0d0
      gc(3)=0.0d0
      gc(4)=0.0d0
      gc(5)=0.0d0
      gc(6)=0.0d0
      gc(7)=0.0d0
      gc(8)=0.0d0
      gc(9)=0.0d0
      gc(10)=0.0d0
      gc(11)=0.0d0
      gc(12)=0.0d0
      gc(13)=0.0d0
      gc(14)=0.0d0
      gc(15)=0.0d0
      gc(16)=0.0d0
      gc(17)=0.0d0
      gc(18)=0.0d0
      gc(19)=0.0d0
      gc(20)=0.0d0
      go to (311,312,313,314,315,316,317,318,341,342,343,344,345,346, &
       347,348,349),kc
  341 gc(1)=1.0d0
      gc(2)=1.0d0
      gc(11)=4.0d0
      gc(12)=-2.1d1
      return
  342 gc(1)=2.0d0*x(1)
      gc(11)=1.5d1
      gc(12)=-8.0d0
      return
  343 gc(1)=4.0d0
      gc(2)=9.0d0
      gc(13)=1.0d1*x(13)
      gc(14)=-9.0d0
      return
  344 gc(1)=3.0d0
      gc(2)=4.0d0
      gc(13)=6.0d0*(x(13)-6.0d0)
      gc(14)=-1.4d1
      return
  345 gc(1)=2.8d1*x(1)
      gc(15)=3.5d1
      gc(16)=-7.9d1
      return
  346 gc(2)=3.0d1*x(2)
      gc(15)=1.1d1
      gc(16)=-6.1d1
      return
  347 gc(1)=1.0d1*x(1)
      gc(2)=2.0d0
      gc(17)=3.6d1*x(17)**3
      gc(18)=-1.0d0
      return
  348 gc(1)=2.0d0*x(1)
      gc(2)=-1.0d0
      gc(19)=1.9d1
      gc(20)=-2.0d1
      return
  349 gc(1)=1.4d1*x(1)
      gc(2)=1.0d1*x(2)
      gc(19)=2.0d0*x(19)
      gc(20)=-3.0d1
      return
      end subroutine tcgu07

! subroutine tffu07             all systems                 90/12/01
! portability : all systems
! 90/12/01 lu : original version
!
! purpose :
!  values of test functions for nonlinear approximation.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ri  x(n)  vector of variables.
!  ro  ff  value of the objective function.
!  ii  next  number of the test problem.
!
      subroutine tffu07(n,x,ff,next)
      integer n,next
      double precision x(n),ff
      double precision x1,x2,x3,x4,x5,x6,x7,x8,a,a1,a2,a3,a4,a5,a6,a7, &
       a8,b,b1,b3,b4,b5,b6,c,p,q
      integer i,j,k
      double precision y(128)
      common /empr07/ y
      go to(10,20,30,40,50,60,70,60,90,100,110,120,130,140,150,160,170, &
       180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330, &
       340),next
   10 ff=-x(1)-x(2)
      return
   20 ff=x(1)-x(2)
      return
   30 ff=0.5d0*x(1)**2+x(2)**2-x(1)*x(2)-7.0d0*x(1)-7.0d0*x(2)
      return
   40 ff=(x(1)-2.0d0)**2+x(2)**2
      return
   50 ff=0.0d0
      do 51 i=1,44
      x1=8.0d0-y(i)
      x2=exp(x1*x(2))
      x3=x(1)-0.49d0
      x4=y(i+44)-x(1)+x2*x3
      ff=ff+x4*x4
   51 continue
      return
   60 ff=100.0d0*(x(2)-x(1)**2)**2+(1.0d0-x(1))**2
      return
   70 ff=(x(1)-10.0d0)**3+(x(2)-20.0d0)**3
      return
   90 ff=-75.196d0+3.8112d0*x(1)+0.0020567d0*x(1)**3- &
       1.0345d-5*x(1)**4+6.8306d0*x(2)-0.030234d0*x(1)*x(2)+ &
       1.28134d-3*x(2)*x(1)**2+2.266d-7*x(2)*x(1)**4- &
       0.25645d0*x(2)**2+0.0034604d0*x(2)**3-1.3514d-5*x(2)**4+ &
       28.106d0/(x(2)+1.0d0)+5.2375d-6*(x(1)*x(2))**2+ &
        6.3d-8*x(1)*(x(1)*x(2))**2-7.0d-10*(x(1)*x(2))**3- &
       3.405d-4*x(1)*x(2)**2+1.6638d-6*x(1)*x(2)**3+ &
      2.8673d0*exp(0.0005d0*x(1)*x(2))-3.5256d-5*x(2)*x(1)**3
      return
  100 ff=-x(1)*x(2)*x(3)
      return
  110 ff=x(1)**2+x(2)**2+x(3)**2
      return
  120 ff=9.0d0*x(1)**2+x(2)**2+9.0d0*x(3)**2
      return
  130 ff=5.0d0*x(1)+5.0d4/x(1)+20.0d0*x(2)+72.0d3/x(2)+ &
       10.0d0*x(3)+14.4d4/x(3)
      return
  140 ff=-x(1)
      return
  150 x2=1.6d0*x(1)
  151 x3=1.22d0*x2-x(1)
      x6=(x(2)+x3)/x(1)
      x1=0.01d0*x(1)*(112.0d0+13.167d0*x6-0.6667d0*x6**2)
      if (abs(x1-x2).gt.0.001d0) then
      x2=x1
      go to 151
      endif
      x4=93.0d0
  152 x5=86.35d0+1.098d0*x6-0.038d0*x6**2+0.325d0*(x4-89.0d0)
      x8=3.0d0*x5-133.0d0
      x7=35.82d0-0.222d0*x8
      x1=98000.0d0*x(3)/(x2*x7+1000.0d0*x(3))
      if (abs(x1-x4).gt.0.001d0) then
      x4=x1
      go to 152
      endif
      ff=-(0.063d0*x2*x5-5.04d0*x(1)-3.36d0*x3-0.035d0*x(2)- &
       10.0d0*x(3))
      return
  160 ff=0.0d0
      b3=x(3)
      b4=1.0d0-x(3)
      c=1.0d0/x(4)
      b=b3+b4*x(4)
      a=b*c
      a1=log(a)+1.0d0
      b1=log(b)+1.0d0
      x1=1.0d0/(12.0d0*x(1)+1.0d0)
      x2=1.0d0/(12.0d0*x(2)+1.0d0)
      x3=sqrt(0.1591545d0*x(1))*12.0d0*x1
      x4=sqrt(0.1591545d0*x(2))*12.0d0*x2
      x5=x3*b4
      x6=x4*b3
      x7=x5*x(1)
      x8=x6*x(2)
      do 161 i=1,19
      a5=x(1)*(y(i+38)+a1-a*y(i))
      a6=exp(a5-y(i+38))
      b5=x(2)*(y(i+38)+b1-b*y(i))
      b6=exp(b5-y(i+38))
      p=x7*a6+x8*b6
      q=p-y(i+19)
      ff=ff+q*q
  161 continue
      return
  170 ff=1.0d0+x(1)+x(2)+x(3)+x(4)
      return
  180 x1=x(1)*x(1)
      x2=x(2)*x(2)
      x3=x(3)*x(3)
      x4=x(4)*x(4)
      x5=x(1)+x(1)
      x6=x(2)+x(2)
      x7=x(3)+x(3)
      x8=x(4)+x(4)
      ff=x1+x2+x3+x3+x4-5.0d0*(x(1)+x(2))-2.1d1*x(3)+7.0d0*x(4)
      return
  190 ff=10.0d0*x(1)*x(4)-6.0d0*x(3)*x(2)**2+x(2)*x(1)**3+ &
       9.0d0*sin(x(5)-x(3))+x(2)**3*x(4)**2*x(5)**4
      return
  200 ff=5.3578547d0*x(3)**2+0.8356891d0*x(1)*x(5)+ &
       37.293239d0*x(1)-40792.141d0
      return
  210 ff=24345.0d0+8720288.849d0*x(1)-150512.5253d0*x(1)*x(2)+ &
       156.6950325d0*x(1)*x(3)-476470.3222d0*x(1)*x(4)- &
       729482.8271d0*x(1)*x(5)
      return
  220 a1=-5.843d-7
      a2=1.17d-4
      a3=2.358d-5
      a4=1.502d-6
      a5=0.0321d0
      a6=0.004324d0
      a7=1.0d-4
      a8=37.48d0
      x1=146.312d3
      y(1)=x(2)+x(3)+41.6d0
      y(18)=0.024d0*x(4)-4.62d0
      y(2)=12.5d0/y(18)+12.0d0
      y(19)=0.0003535d0*x(1)**2+0.5311d0*x(1)+0.08705d0*y(2)*x(1)
      y(20)=0.052d0*x(1)+78.0d0+0.002377d0*y(2)*x(1)
      y(3)=y(19)/y(20)
      y(4)=19.0d0*y(3)
      y(21)=0.04782d0*(x(1)-y(3))+0.1956d0*(x(1)-y(3))**2/x(2)+ &
       0.6376d0*y(4)+1.594d0*y(3)
      y(22)=100.0d0*x(2)
      y(23)=x(1)-y(3)-y(4)
      y(24)=0.95d0-y(21)/y(22)
      y(5)=y(23)*y(24)
      y(6)=x(1)-y(5)-y(4)-y(3)
      y(25)=(y(5)+y(4))*0.995d0
      y(7)=y(25)/y(1)
      y(8)=y(25)/3798.0d0
      y(26)=y(7)-0.0663d0*y(7)/y(8)-0.3153d0
      y(9)=96.82d0/y(26)+0.321d0*y(1)
      y(10)=1.29d0*y(5)+1.258d0*y(4)+2.29d0*y(3)+1.71d0*y(6)
      y(11)=1.71d0*x(1)-0.452d0*y(4)+0.58d0*y(3)
      y(27)=12.3d0/752.3d0
      y(28)=1.75d0*y(2)*0.995d0*x(1)
      y(29)=0.995d0*y(10)+1998.0d0
      y(12)=y(27)*x(1)+y(28)/y(29)
      y(13)=y(29)-1.75d0*y(2)
      y(14)=3623.0d0+64.4d0*x(2)+58.4d0*x(3)+x1/(y(9)+ &
       x(5))
      y(30)=0.995d0*y(10)+60.8d0*x(2)+48.0d0*x(4)-0.1121d0*y(14)- &
       5095.0d0
      y(15)=y(13)/y(30)
      y(16)=148.0d3-331.0d3*y(15)+40.0d0*y(13)-61.0d0*y(15)*y(13)
      y(31)=2324.0d0*y(10)-2.874d7*y(2)
      y(17)=1.413d7-1328.0d0*y(10)-531.0d0*y(11)+y(31)/y(29)
      y(32)=y(13)/y(15)-y(13)/0.52d0
      y(33)=1.104d0-0.72d0*y(15)
      y(34)=y(9)+x(5)
      ff=a1*y(17)+a2*y(14)+a3*y(13)+a4*y(16)+a5*y(12)+a6*y(5)+ &
       a7*y(32)/y(33)+a8*y(2)/y(29)+0.1365d0
      return
  230 ff=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
      return
  240 ff=0.0204d0*x(1)*x(4)*(x(1)+x(2)+x(3))+0.0187d0*x(2)*x(3)* &
       (x(1)+1.57d0*x(2)+x(4))+0.0607d0*x(1)*x(4)*x(5)**2*(x(1)+ &
       x(2)+x(3))+0.0437d0*x(2)*x(3)*x(6)**2*(x(1)+1.57d0*x(2)+x(4))
      return
  250 ff=4.3d0*x(1)+31.8d0*x(2)+63.3d0*x(3)+15.8d0*x(4)+ &
       68.5d0*x(5)+4.7d0*x(6)
      return
  260 ff=(x(1)-1.0d1)**2+5.0d0*(x(2)-1.2d1)**2+x(3)**4+3.0d0* &
      (x(4)-1.1d1)**2+1.0d1*x(5)**6+7.0d0*x(6)**2+x(7)**4-4.0d0* &
      x(6)*x(7)-1.0d1*x(6)-8.0d0*x(7)
      return
  270 ff=10.0d0*x(1)*x(4)**2/(x(2)*x(6)**3*x(7)**0.25d0)+ &
       15.0d0*x(3)*x(4)/(x(1)*x(2)**2*x(5)*x(7)**0.5d0)+ &
       20.0d0*x(2)*x(6)/(x(1)**2*x(4)*x(5)**2)+ &
       25.0d0*x(1)**2*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6)**2)
      return
  280 ff=0.4d0*(x(1)/x(7))**0.67d0+0.4d0*(x(2)/x(8))**0.67d0+ &
       10.0d0-x(1)-x(2)
      return
  290 ff=x(1)+x(2)+x(3)
      return
  300 ff=-0.5d0*(x(1)*x(4)-x(2)*x(3)+x(3)*x(9)-x(5)*x(9)+x(5)*x(8)- &
       x(6)*x(7))
      return
  310 ff=x(1)**2+x(2)**2+x(1)*x(2)-1.4d1*x(1)-1.6d1*x(2)+ &
      (x(3)-1.0d1)**2+4.0d0*(x(4)-5.0d0)**2+(x(5)-3.0d0)**2+ &
      2.0d0*(x(6)-1.0d0)**2+5.0d0*x(7)**2+7.0d0*(x(8)- &
      1.1d1)**2+2.0d0*(x(9)-1.0d1)**2+(x(10)-7.0d0)**2+4.5d1
      return
  320 ff=0.0d0
      do 321 i=1,10
      ff=ff-y(i+50)*x(i)
  321 continue
      k=0
      do 323 i=1,5
      x1=2.0d0*y(i+85)*x(i+10)**2
      x2=0.0d0
      do 322 j=1,5
      x2=x2+y(k+j+60)*x(j+10)
  322 continue
      ff=ff+(x1+x2)*x(i+10)
      k=k+5
  323 continue
      return
  330 ff=0.0d0
      do 331 i=1,5
      ff=ff+y(i)*x(i+11)+y(i+5)*x(i)*x(i+11)
  331 continue
      return
  340 ff=x(1)**2+x(2)**2+x(1)*x(2)-1.4d1*x(1)-1.6d1*x(2)+(x(3)- &
      1.0d1)**2+4.0d0*(x(4)-5.0d0)**2+(x(5)-3.0d0)**2+2.0d0* &
      (x(6)-1.0d0)**2+5.0d0*x(7)**2+7.0d0*(x(8)-1.1d1)**2+ &
      2.0d0*(x(9)-1.0d1)**2+(x(10)-7.0d0)**2+(x(11)-9.0d0)**2+ &
      1.0d1*(x(12)-1.0d0)**2+5.0d0*(x(13)-7.0d0)**2+4.0d0* &
      (x(14)-1.4d1)**2+2.7d1*(x(15)-1.0d0)**2+x(16)**4+(x(17)- &
      2.0d0)**2+1.3d1*(x(18)-2.0d0)**2+(x(19)-3.d0)**2+x(20)**2+ &
      9.5d1
      return
      end subroutine tffu07

! subroutine tfgu07             all systems                 90/12/01
! portability : all systems
! 90/12/01 lu : original version
!
! purpose :
!  gradients of test functions for nonlinear approximation.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ri  x(n)  vector of variables.
!  ro  gf(n)  gradient of the objective function.
!  ii  next  number of the test problem.
!
      subroutine tfgu07(n,x,gf,next)
      integer n,next
      double precision x(n),gf(n)
      double precision x1,x2,x3,x4,x5,x6,x7,x8,a,a1,a2,a3,a4,a5,a6,a7, &
       a8,b,b1,b2,b3,b4,b5,b6,b7,c,p,q
      integer i,j,k
      double precision y(128)
      common /empr07/ y
      go to(10,20,30,40,50,60,70,60,90,100,110,120,130,140,150,160,170, &
       180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330, &
       340),next
   10 gf(1)=-1.0d0
      gf(2)=-1.0d0
      return
   20 gf(1)= 1.0d0
      gf(2)=-1.0d0
      return
   30 gf(1)=x(1)-x(2)-7.0d0
      gf(2)=2.0d0*x(2)-x(1)-7.0d0
      return
   40 gf(1)=2.0d0*(x(1)-2.0d0)
      gf(2)=2.0d0*x(2)
      return
   50 gf(1)=0.0d0
      gf(2)=0.0d0
      do 51 i=1,44
      x1=8.0d0-y(i)
      x2=exp(x1*x(2))
      x3=x(1)-0.49d0
      x4=y(i+44)-x(1)+x2*x3
      gf(1)=gf(1)+x4*(x2-1.0d0)
      gf(2)=gf(2)+x1*x2*x3*x4
   51 continue
      gf(1)=2.0d0*gf(1)
      gf(2)=2.0d0*gf(2)
      return
   60 gf(1)=-400.0d0*(x(2)-x(1)**2)*x(1)-2.0d0*(1.0d0-x(1))
      gf(2)=200.0d0*(x(2)-x(1)**2)
      return
   70 gf(1)=3.0d0*(x(1)-10.0d0)**2
      gf(2)=3.0d0*(x(2)-20.0d0)**2
      return
   90 gf(1)=3.8112d0+0.0061701d0*x(1)**2-4.1380d-5*x(1)**3- &
       0.030234d0*x(2)+2.56268d-3*x(1)*x(2)+9.064d-7*x(2)*x(1)**3+ &
       10.4750d-6*x(1)*x(2)**2+18.9d-8*(x(1)*x(2))**2- &
       21.0d-10*(x(1)*x(2))**2*x(2)-3.405d-4*x(2)**2+ &
       1.6638d-6*x(2)**3-10.5768d-5*x(2)*x(1)**2+ &
       2.8673d0*0.0005d0*exp(0.0005d0*x(1)*x(2))*x(2)
      gf(2)=6.8306d0-0.030234d0*x(1)+1.28134d-3*x(1)**2+ &
       2.266d-7*x(1)**4-0.51290d0*x(2)+0.0103812d0*x(2)**2- &
       5.4056d-5*x(2)**3-28.106d0/(x(2)+1.0d0)**2+ &
       10.4750d-6*x(2)*x(1)**2+12.6d-8*x(2)*x(1)**3- &
       21.0d-10*x(1)*(x(1)*x(2))**2-6.810d-4*x(1)*x(2)+ &
       4.9914d-6*x(1)*x(2)**2-3.5256d-5*x(1)**3+ &
       2.8673d0*0.0005d0*exp(0.0005d0*x(1)*x(2))*x(1)
      return
  100 gf(1)=-x(2)*x(3)
      gf(2)=-x(1)*x(3)
      gf(3)=-x(1)*x(2)
      return
  110 gf(1)=2.0d0*x(1)
      gf(2)=2.0d0*x(2)
      gf(3)=2.0d0*x(3)
      return
  120 gf(1)=18.0d0*x(1)
      gf(2)= 2.0d0*x(2)
      gf(3)=18.0d0*x(3)
      return
  130 gf(1)=5.0d0-5.0d4/x(1)**2
      gf(2)=2.0d1-72.0d3/x(2)**2
      gf(3)=1.0d1-14.4d4/x(3)**2
      return
  140 gf(1)=-1.0d0
      gf(2)= 0.0d0
      gf(3)= 0.0d0
      return
  150 x2=1.6d0*x(1)
      y(4)=1.6d0
      y(5)=0.0d0
      y(6)=0.0d0
  151 x3=1.22d0*x2-x(1)
      y(7)=1.22d0*y(4)-1.0d0
      y(8)=1.22d0*y(5)
      y(9)=1.22d0*y(6)
      x6=(x(2)+x3)/x(1)
      y(16)=y(7)/x(1)-(x(2)+x3)/x(1)**2
      y(17)=(1.0d0+y(8))/x(1)
      y(18)=y(9)/x(1)
      x1=0.01d0*x(1)*(112.0d0+13.167d0*x6-0.6667d0*x6**2)
      y(1)=0.01d0*(112.0d0+13.167d0*x6-0.6667d0*x6**2)+ &
       0.01d0*x(1)*(13.167*y(16)-1.3334d0*x6*y(16))
      y(2)=0.01d0*x(1)*(13.167d0*y(17)-1.3334d0*x6*y(17))
      y(3)=0.01d0*x(1)*(13.167d0*y(18)-1.3334d0*x6*y(18))
      if (abs(x1-x2).gt.0.001d0) then
      x2=x1
      y(4)=y(1)
      y(5)=y(2)
      y(6)=y(3)
      go to 151
      endif
      x4=93.0d0
      y(10)=0.0d0
      y(11)=0.0d0
      y(12)=0.0d0
  152 x5=86.35d0+1.098d0*x6-0.038d0*x6**2+0.325d0*(x4-89.0d0)
      y(13)=1.098d0*y(16)-0.076d0*x6*y(16)+0.325d0*y(10)
      y(14)=1.098d0*y(17)-0.076d0*x6*y(17)+0.325d0*y(11)
      y(15)=1.098d0*y(18)-0.076d0*x6*y(18)+0.325d0*y(12)
      x8=3.0d0*x5-133.0d0
      y(22)=3.0d0*y(13)
      y(23)=3.0d0*y(14)
      y(24)=3.0d0*y(15)
      x7=35.82d0-0.222d0*x8
      y(19)=-0.222d0*y(22)
      y(20)=-0.222d0*y(23)
      y(21)=-0.222d0*y(24)
      y(89)=x2*x7+1000.0d0*x(3)
      y(90)=y(89)**2
      x1=98000.0d0*x(3)/y(89)
      y(1)=-98000.0d0*x(3)*(y(4)*x7+x2*y(19))/y(90)
      y(2)=-98000.0d0*x(3)*(y(5)*x7+x2*y(20))/y(90)
      y(3)=-98000.0d0*x(3)*(y(6)*x7+x2*y(21)+1000.0d0)/y(90)+ &
       98000.0d0/y(89)
      if (abs(x1-x4).gt.0.001d0) then
      x4=x1
      y(10)=y(1)
      y(11)=y(2)
      y(12)=y(3)
      go to 152
      endif
      gf(1)=-(0.063d0*(y(4)*x5+x2*y(13))-3.36d0*y(7)-5.04d0)
      gf(2)=-(0.063d0*(y(5)*x5+x2*y(14))-3.36d0*y(8)-0.35d-1)
      gf(3)=-(0.063d0*(y(6)*x5+x2*y(15))-3.36d0*y(9)-10.0d0)
      return
  160 gf(1)=0.0d0
      gf(2)=0.0d0
      gf(3)=0.0d0
      gf(4)=0.0d0
      b3=x(3)
      b4=1.0d0-x(3)
      c=1.0d0/x(4)
      b=b3+b4*x(4)
      a=b*c
      a1=log(a)+1.0d0
      b1=log(b)+1.0d0
      x1=1.0d0/(12.0d0*x(1)+1.0d0)
      x2=1.0d0/(12.0d0*x(2)+1.0d0)
      x3=sqrt(0.1591545d0*x(1))*12.0d0*x1
      x4=sqrt(0.1591545d0*x(2))*12.0d0*x2
      x5=x3*b4
      x6=x4*b3
      x7=x5*x(1)
      x8=x6*x(2)
      a2=1.0d0/a
      b2=1.0d0/b
      b3=1.0d0-x(4)
      a3=c*b3
      a4=-c*c*x(3)
      x1=x1+0.5d0
      x2=x2+0.5d0
      x3=x3*x(1)
      x4=x4*x(2)
      a5=x7*x(1)
      b5=x8*x(2)
      a3=a5*a3
      b3=b5*b3
      a4=a5*a4
      b4=b5*b4
      do 161 i=1,19
      a5=x(1)*(y(i+38)+a1-a*y(i))
      a6=exp(a5-y(i+38))
      b5=x(2)*(y(i+38)+b1-b*y(i))
      b6=exp(b5-y(i+38))
      p=x7*a6+x8*b6
      q=p-y(i+19)
      a7=a6*(a2-y(i))
      b7=b6*(b2-y(i))
      gf(1)=gf(1)+q*a6*x5*(x1+a5)
      gf(2)=gf(2)+q*b6*x6*(x2+b5)
      gf(3)=gf(3)+q*(a7*a3+b7*b3-a6*x3+b6*x4)
      gf(4)=gf(4)+q*(a7*a4+b7*b4)
  161 continue
      gf(1)=2.0d0*gf(1)
      gf(2)=2.0d0*gf(2)
      gf(3)=2.0d0*gf(3)
      gf(4)=2.0d0*gf(4)
      return
  170 gf(1)=1.0d0
      gf(2)=1.0d0
      gf(3)=1.0d0
      gf(4)=1.0d0
      return
  180 x5=x(1)+x(1)
      x6=x(2)+x(2)
      x7=x(3)+x(3)
      x8=x(4)+x(4)
      gf(1)=x5-5.0d0
      gf(2)=x6-5.0d0
      gf(3)=x7+x7-2.1d1
      gf(4)=x8+7.0d0
      return
  190 gf(1)=10.0d0*x(4)+3.0d0*x(2)*x(1)**2
      gf(2)=-12.0d0*x(2)*x(3)+x(1)**3+ &
       3.0d0*(x(2)*x(4))**2*x(5)**4
      gf(3)=-6.0d0*x(2)**2-9.0d0*cos(x(5)-x(3))
      gf(4)=10.0d0*x(1)+2.0d0*x(4)*x(2)**3*x(5)**4
      gf(5)=9.0d0*cos(x(5)-x(3))+4.0d0*x(2)**3*x(4)**2*x(5)**3
      return
  200 gf(1)=0.8356891d0*x(5)+37.293239d0
      gf(2)=0.0d0
      gf(3)=10.7157094d0*x(3)
      gf(4)=0.0d0
      gf(5)=0.8356891d0*x(1)
      return
  210 gf(1)= 8720288.849d0-150512.5253d0*x(2)+156.6950325d0*x(3)- &
       476470.3222d0*x(4)-729482.8271d0*x(5)
      gf(2)=-150512.5253d0*x(1)
      gf(3)= 156.6950325d0*x(1)
      gf(4)=-476470.3222d0*x(1)
      gf(5)=-729482.8271d0*x(1)
      return
  220 a1=-5.843d-7
      a2=1.17d-4
      a3=2.358d-5
      a4=1.502d-6
      a5=0.0321d0
      a6=0.004324d0
      a7=1.0d-4
      a8=37.48d0
      x1=146.312d3
      y(1)=x(2)+x(3)+41.6d0
      y(18)=0.024d0*x(4)-4.62d0
      y(2)=12.5d0/y(18)+12.0d0
      y(19)=0.0003535d0*x(1)**2+0.5311d0*x(1)+0.08705d0*y(2)*x(1)
      y(20)=0.052d0*x(1)+78.0d0+0.002377d0*y(2)*x(1)
      y(3)=y(19)/y(20)
      y(4)=19.0d0*y(3)
      y(21)=0.04782d0*(x(1)-y(3))+0.1956d0*(x(1)-y(3))**2/x(2)+ &
       0.6376d0*y(4)+1.594d0*y(3)
      y(22)=100.0d0*x(2)
      y(23)=x(1)-y(3)-y(4)
      y(24)=0.95d0-y(21)/y(22)
      y(5)=y(23)*y(24)
      y(6)=x(1)-y(5)-y(4)-y(3)
      y(25)=(y(5)+y(4))*0.995d0
      y(7)=y(25)/y(1)
      y(8)=y(25)/3798.0d0
      y(26)=y(7)-0.0663d0*y(7)/y(8)-0.3153d0
      y(9)=96.82d0/y(26)+0.321d0*y(1)
      y(10)=1.29d0*y(5)+1.258d0*y(4)+2.29d0*y(3)+1.71d0*y(6)
      y(11)=1.71d0*x(1)-0.452d0*y(4)+0.58d0*y(3)
      y(27)=12.3d0/752.3d0
      y(28)=1.75d0*y(2)*0.995d0*x(1)
      y(29)=0.995d0*y(10)+1998.0d0
      y(12)=y(27)*x(1)+y(28)/y(29)
      y(13)=y(29)-1.75d0*y(2)
      y(14)=3623.0d0+64.4d0*x(2)+58.4d0*x(3)+x1/(y(9)+ &
       x(5))
      y(30)=0.995d0*y(10)+60.8d0*x(2)+48.0d0*x(4)-0.1121d0*y(14)- &
       5095.0d0
      y(15)=y(13)/y(30)
      y(16)=148.0d3-331.0d3*y(15)+40.0d0*y(13)-61.0d0*y(15)*y(13)
      y(31)=2324.0d0*y(10)-2.874d7*y(2)
      y(17)=1.413d7-1328.0d0*y(10)-531.0d0*y(11)+y(31)/y(29)
      y(32)=y(13)/y(15)-y(13)/0.52d0
      y(33)=1.104d0-0.72d0*y(15)
      y(34)=y(9)+x(5)
      y(35)=1.0d0
      y(36)=0.024d0
      y(37)=-12.5d0*y(36)/y(18)**2
      y(38)=0.000707d0*x(1)+0.5311d0+0.08705d0*y(2)
      y(39)=0.08705d0*x(1)*y(37)
      y(40)=0.052d0+0.002377d0*y(2)
      y(41)=0.002377d0*x(1)*y(37)
      y(42)=(y(38)*y(20)-y(19)*y(40))/y(20)**2
      y(43)=(y(39)*y(20)-y(19)*y(41))/y(20)**2
      y(44)=19.0d0*y(42)
      y(45)=19.0d0*y(43)
      y(46)=0.04782d0*(1.0d0-y(42))+0.3912d0*(x(1)-y(3))*(1.0d0- &
       y(42))/x(2)+0.6376d0*y(44)+1.594d0*y(42)
      y(47)=-0.1956d0*(x(1)-y(3))**2/x(2)**2
      y(48)=-0.04782d0*y(43)-0.3912d0*(x(1)-y(3))*y(43)/x(2)+ &
       0.6376d0*y(45)+1.594d0*y(43)
      y(49)=100.0d0
      y(50)=1.0d0-y(42)-y(44)
      y(51)=-y(43)-y(45)
      y(52)=-y(46)/y(22)
      y(53)=(y(21)*y(49)-y(47)*y(22))/y(22)**2
      y(54)=-y(48)/y(22)
      y(55)=y(50)*y(24)+y(23)*y(52)
      y(56)=y(23)*y(53)
      y(57)=y(51)*y(24)+y(23)*y(54)
      y(58)=1.0d0-y(55)-y(44)-y(42)
      y(59)=-y(56)
      y(60)=-y(57)-y(45)-y(43)
      y(61)=0.995d0*(y(55)+y(44))
      y(62)=0.995d0*y(56)
      y(63)=0.995d0*(y(57)+y(45))
      y(64)=y(61)/y(1)
      y(65)=(y(62)*y(1)-y(25)*y(35))/y(1)**2
      y(66)=-y(25)*y(35)/y(1)**2
      y(67)=y(63)/y(1)
      y(68)=y(61)/3798.0d0
      y(69)=y(62)/3798.0d0
      y(70)=y(63)/3798.0d0
      y(71)=y(64)-0.0663d0*(y(64)*y(8)-y(7)*y(68))/y(8)**2
      y(72)=y(65)-0.0663d0*(y(65)*y(8)-y(7)*y(69))/y(8)**2
      y(73)=y(66)-0.0663d0*y(66)/y(8)
      y(74)=y(67)-0.0663d0*(y(67)*y(8)-y(7)*y(70))/y(8)**2
      y(75)=-96.82d0*y(71)/y(26)**2
      y(76)=-96.82d0*y(72)/y(26)**2+0.321d0*y(35)
      y(77)=-96.82d0*y(73)/y(26)**2+0.321d0*y(35)
      y(78)=-96.82d0*y(74)/y(26)**2
      y(79)=1.29d0*y(55)+1.258d0*y(44)+2.29d0*y(42)+1.71d0*y(58)
      y(80)=1.29d0*y(56)+1.71d0*y(59)
      y(81)=1.29d0*y(57)+1.258d0*y(45)+2.29d0*y(43)+1.71d0*y(60)
      y(82)=1.71d0-0.452d0*y(44)+0.58d0*y(42)
      y(83)=-0.452d0*y(45)+0.58d0*y(43)
      y(84)=1.75d0*y(2)*0.995d0
      y(85)=1.75d0*y(37)*0.995d0*x(1)
      y(86)=0.995d0*y(79)
      y(87)=0.995d0*y(80)
      y(88)=0.995d0*y(81)
      y(89)=y(27)+(y(84)*y(29)-y(28)*y(86))/y(29)**2
      y(90)=-y(28)*y(87)/y(29)**2
      y(91)=(y(85)*y(29)-y(28)*y(88))/y(29)**2
      y(92)=y(88)-1.75d0*y(37)
      y(93)=-x1*y(75)/(y(9)+x(5))**2
      y(94)=64.4d0-x1*y(76)/(y(9)+x(5))**2
      y(95)=58.4d0-x1*y(77)/(y(9)+x(5))**2
      y(96)=-x1*y(78)/(y(9)+x(5))**2
      y(97)=-x1/(y(9)+x(5))**2
      y(98)=0.995d0*y(79)-0.1121d0*y(93)
      y(99)=0.995d0*y(80)+60.8d0-0.1121d0*y(94)
      y(100)=-0.1121d0*y(95)
      y(101)=0.995d0*y(81)+48.0d0-0.1121d0*y(96)
      y(102)=-0.1121d0*y(97)
      y(103)=(y(86)*y(30)-y(13)*y(98))/y(30)**2
      y(104)=(y(87)*y(30)-y(13)*y(99))/y(30)**2
      y(105)=-y(13)*y(100)/y(30)**2
      y(106)=(y(92)*y(30)-y(13)*y(101))/y(30)**2
      y(107)=-y(13)*y(102)/y(30)**2
      y(108)=-3.31d5*y(103)+40.0d0*y(86)-61.0d0*(y(103)*y(13)+ &
       y(15)*y(86))
      y(109)=-3.31d5*y(104)+40.0d0*y(87)-61.0d0*(y(104)*y(13)+ &
       y(15)*y(87))
      y(110)=-3.31d5*y(105)-61.0d0*y(105)*y(13)
      y(111)=-3.31d5*y(106)+40.0d0*y(92)-61.0d0*(y(106)*y(13)+ &
       y(15)*y(92))
      y(112)=-3.31d5*y(107)-61.0d0*y(107)*y(13)
      y(113)=2.324d3*y(79)
      y(114)=2.324d3*y(80)
      y(115)=2.324d3*y(81)-2.874d7*y(37)
      y(116)=-1.328d3*y(79)-531.0d0*y(82)+(y(113)*y(29)-y(31)* &
       y(86))/y(29)**2
      y(117)=-1.328d3*y(80)+(y(114)*y(29)-y(31)*y(87))/y(29)**2
      y(118)=-1.328d3*y(81)-531.0d0*y(83)+(y(115)*y(29)-y(31)* &
       y(88))/y(29)**2
      y(119)=(y(86)*y(15)-y(13)*y(103))/y(15)**2-y(86)/0.52d0
      y(120)=(y(87)*y(15)-y(13)*y(104))/y(15)**2-y(87)/0.52d0
      y(121)=-y(13)*y(105)/y(15)**2
      y(122)=(y(92)*y(15)-y(13)*y(106))/y(15)**2-y(92)/0.52d0
      y(123)=-y(13)*y(107)/y(15)**2
      y(124)=-0.72d0*y(103)
      y(125)=-0.72d0*y(104)
      y(126)=-0.72d0*y(105)
      y(127)=-0.72d0*y(106)
      y(128)=-0.72d0*y(107)
      gf(1)=a1*y(116)+a2*y(93)+a3*y(86)+a4*y(108)+a5*y(89)+a6*y(55)+ &
       a7*(y(119)*y(33)-y(32)*y(124))/y(33)**2-a8*y(2)*y(86)/y(29)**2
      gf(2)=a1*y(117)+a2*y(94)+a3*y(87)+a4*y(109)+a5*y(90)+a6*y(56)+ &
       a7*(y(120)*y(33)-y(32)*y(125))/y(33)**2-a8*y(2)*y(87)/y(29)**2
      gf(3)=a2*y(95)+a4*y(110)+a7*(y(121)*y(33)-y(32)*y(126))/y(33)**2
      gf(4)=a1*y(118)+a2*y(96)+a3*y(92)+a4*y(111)+a5*y(91)+a6*y(57)+ &
       a7*(y(122)*y(33)-y(32)*y(127))/y(33)**2+a8*(y(37)*y(29)-y(2)* &
       y(88))/y(29)**2
      gf(5)=a2*y(97)+a4*y(112)+a7*(y(123)*y(33)-y(32)*y(128))/y(33)**2
      return
  230 gf(1)= 2.0d0*(x(1)-x(4))
      gf(2)= 2.0d0*(x(2)-x(5))
      gf(3)= 2.0d0*(x(3)-x(6))
      gf(4)=-2.0d0*(x(1)-x(4))
      gf(5)=-2.0d0*(x(2)-x(5))
      gf(6)=-2.0d0*(x(3)-x(6))
      return
  240 gf(1)=0.0408d0*x(1)*x(4)+0.0204d0*x(4)*(x(2)+x(3))+ &
       0.0187d0*x(2)*x(3)+0.1214d0*x(1)*x(4)*x(5)**2+ &
       0.0607d0*x(4)*x(5)**2*(x(2)+x(3))+0.0437d0*x(2)*x(3)*x(6)**2
      gf(2)=0.0204d0*x(1)*x(4)+0.0374d0*1.57d0*x(2)*x(3)+ &
       0.0187d0*x(3)*(x(1)+x(4))+0.0607d0*x(1)*x(4)*x(5)**2+ &
       0.0874d0*1.57d0*x(2)*x(3)*x(6)**2+ &
       0.0437d0*x(3)*x(6)**2*(x(1)+x(4))
      gf(3)=0.0204d0*x(1)*x(4)+ &
       0.0187d0*x(2)*(x(1)+1.57d0*x(2)+x(4))+ &
       0.0607d0*x(1)*x(4)*x(5)**2+ &
       0.0437d0*x(2)*x(6)**2*(x(1)+1.57d0*x(2)+x(4))
      gf(4)=0.0204d0*x(1)*(x(1)+x(2)+x(3))+ &
       0.0187d0*x(2)*x(3)+ &
       0.0607d0*x(1)*x(5)**2*(x(1)+x(2)+x(3))+ &
       0.0437d0*x(2)*x(3)*x(6)**2
      gf(5)=0.1214d0*x(1)*x(4)*x(5)*(x(1)+x(2)+x(3))
      gf(6)=0.0874d0*x(2)*x(3)*x(6)*(x(1)+1.57d0*x(2)+x(4))
      return
  250 gf(1)= 4.3d0
      gf(2)=31.8d0
      gf(3)=63.3d0
      gf(4)=15.8d0
      gf(5)=68.5d0
      gf(6)= 4.7d0
      return
  260 gf(1)=2.0d0*(x(1)-1.0d1)
      gf(2)=1.0d1*(x(2)-1.2d1)
      gf(3)=4.0d0*x(3)**3
      gf(4)=6.0d0*(x(4)-1.1d1)
      gf(5)=6.0d1*x(5)**5
      gf(6)=1.4d1*x(6)-4.0d0*x(7)-1.0d1
      gf(7)=4.0d0*x(7)**3-4.0d0*x(6)-8.0d0
      return
  270 gf(1)=10.0d0*x(4)**2/(x(2)*x(6)**3*x(7)**0.25d0)- &
       15.0d0*x(3)*x(4)/((x(1)*x(2))**2*x(5)*x(7)**0.5d0)- &
       40.0d0*x(2)*x(6)/(x(1)**3*x(4)*x(5)**2)+ &
       50.0d0*x(1)*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6)**2)
      gf(2)=-10.0d0*x(1)*x(4)**2/(x(2)**2*x(6)**3*x(7)**0.25d0)- &
       30.0d0*x(3)*x(4)/(x(1)*x(2)**3*x(5)*x(7)**0.5d0)+ &
       20.0d0*x(6)/(x(1)**2*x(4)*x(5)**2)+ &
       50.0d0*x(1)**2*x(2)*x(5)**0.5d0*x(7)/(x(3)*x(6)**2)
      gf(3)=15.0d0*x(4)/(x(1)*x(2)**2*x(5)*x(7)**0.5d0)- &
       25.0d0*x(1)**2*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6))**2
      gf(4)=20.0d0*x(1)*x(4)/(x(2)*x(6)**3*x(7)**0.25d0)+ &
       15.0d0*x(3)/(x(1)*x(2)**2*x(5)*x(7)**0.5d0)- &
       20.0d0*x(2)*x(6)/(x(1)*x(4)*x(5))**2
      gf(5)=-15.0d0*x(3)*x(4)/(x(1)*(x(2)*x(5))**2*x(7)**0.5d0)- &
       40.0d0*x(2)*x(6)/(x(1)**2*x(4)*x(5)**3)+ &
       12.5d0*x(1)**2*x(2)**2*x(7)/(x(3)*x(5)**0.5d0*x(6)**2)
      gf(6)=-30.0d0*x(1)*x(4)**2/(x(2)*x(6)**4*x(7)**0.25d0)+ &
       20.0d0*x(2)/(x(1)**2*x(4)*x(5)**2)- &
       50.0d0*x(1)**2*x(2)**2*x(5)**0.5d0*x(7)/(x(3)*x(6)**3)
      gf(7)=-2.5d0*x(1)*x(4)**2/(x(2)*x(6)**3*x(7)**1.25d0)- &
       7.5d0*x(3)*x(4)/(x(1)*x(2)**2*x(5)*x(7)**1.5d0)+ &
       25.0d0*x(1)**2*x(2)**2*x(5)**0.5d0/(x(3)*x(6)**2)
      return
  280 gf(1)=0.268d0*(x(7)/x(1))**0.33d0/x(7)-1.0d0
      gf(2)=0.268d0*(x(8)/x(2))**0.33d0/x(8)-1.0d0
      gf(3)=0.0d0
      gf(4)=0.0d0
      gf(5)=0.0d0
      gf(6)=0.0d0
      gf(7)=-0.268d0*(x(1)/x(7))**0.67d0/x(7)
      gf(8)=-0.268d0*(x(2)/x(8))**0.67d0/x(8)
      return
  290 gf(1)=1.0d0
      gf(2)=1.0d0
      gf(3)=1.0d0
      gf(4)=0.0d0
      gf(5)=0.0d0
      gf(6)=0.0d0
      gf(7)=0.0d0
      gf(8)=0.0d0
      return
  300 gf(1)=-0.5d0*x(4)
      gf(2)= 0.5d0*x(3)
      gf(3)=-0.5d0*(x(9)-x(2))
      gf(4)=-0.5d0*x(1)
      gf(5)=-0.5d0*(x(8)-x(9))
      gf(6)= 0.5d0*x(7)
      gf(7)= 0.5d0*x(6)
      gf(8)=-0.5d0*x(5)
      gf(9)=-0.5d0*(x(3)-x(5))
      return
  310 gf(1)=2.0d0*x(1)+x(2)-1.4d1
      gf(2)=2.0d0*x(2)+x(1)-1.6d1
      gf(3)=2.0d0*(x(3)-1.0d1)
      gf(4)=8.0d0*(x(4)-5.0d0)
      gf(5)=2.0d0*(x(5)-3.0d0)
      gf(6)=4.0d0*(x(6)-1.0d0)
      gf(7)=1.0d1*x(7)
      gf(8)=1.4d1*(x(8)-1.1d1)
      gf(9)=4.0d0*(x(9)-1.0d1)
      gf(10)=2.0d0*(x(10)-7.0d0)
      return
  320 do 321 i=1,10
      gf(i)=-y(i+50)
  321 continue
      k=0
      do 323 i=1,5
      x1=2.0d0*y(i+85)*x(i+10)**2
      x2=0.0d0
      do 322 j=1,5
      x2=x2+y(k+j+60)*x(j+10)
  322 continue
      gf(i+10)=3.0d0*x1+2.0d0*x2
      k=k+5
  323 continue
      return
  330 do 331 i=1,5
      gf(i)=y(i+5)*x(i+11)
      gf(i+5)=0.0d0
      gf(i+11)=y(i)+y(i+5)*x(i)
  331 continue
      gf(11)=0.0d0
      return
  340 gf(1)=2.0d0*x(1)+x(2)-1.4d1
      gf(2)=2.0d0*x(2)+x(1)-1.6d1
      gf(3)=2.0d0*(x(3)-1.0d1)
      gf(4)=8.0d0*(x(4)-5.0d0)
      gf(5)=2.0d0*(x(5)-3.0d0)
      gf(6)=4.0d0*(x(6)-1.0d0)
      gf(7)=1.0d1*x(7)
      gf(8)=1.4d1*(x(8)-1.1d1)
      gf(9)=4.0d0*(x(9)-1.0d1)
      gf(10)=2.0d0*(x(10)-7.0d0)
      gf(11)=2.0d0*(x(11)-9.0d0)
      gf(12)=2.0d1*(x(12)-1.0d0)
      gf(13)=1.0d1*(x(13)-7.0d0)
      gf(14)=8.0d0*(x(14)-1.4d1)
      gf(15)=5.4d1*(x(15)-1.0d0)
      gf(16)=4.0d0*x(16)**3
      gf(17)=2.0d0*(x(17)-2.0d0)
      gf(18)=2.6d1*(x(18)-2.0d0)
      gf(19)=2.0d0*(x(19)-3.0d0)
      gf(20)=2.0d0*x(20)
      return
      end subroutine tfgu07

! subroutine tiud15                all systems                00/12/01
! portability : all systems
! 92/12/01 ra : original version
!
! purpose :
!  initial values of the variables and structure of the sparse hessian
!  matrix for unconstrained minimization.
!  sparse version with changed tests 7-10, 12.
!  changed for the tests.
!
! parameters :
!  ii  n  number of variables.
!  ii  nb  number of elements of the sparse matrix.
!  ro  x(n)  vector of variables.
!  ii  next  number of the test problem.
!  io  ierr  error indicator.
!
      subroutine tiud15(n,nb,x,xmax,next,ierr)
      integer n,nb,next,ierr
      double precision x(n),xmax
      integer i
      double precision y(20)
      common /empr15/ y
      double precision eta9
      parameter (eta9=1.0d60)
      xmax=1.0d3
      ierr=0
      goto(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160, &
       170,180,190,200,210,220,260,270),next
   10 if (n.lt.2) go to 999
      n=n-mod(n,2)
      do 11 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-1.2d0
        else
          x(i)=1.0d0
        endif
   11 continue
      nb=2*(n-1)
      return
   20 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 21 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-2.0d0
          if(i.le.4) x(i)=-3.0d0
        else
          x(i)=0.0d0
          if(i.le.4) x(i)=-1.0d0
        endif
   21 continue
      nb=3*(n-2)
      return
   30 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 31 i=1,n
        if(mod(i,4).eq.1) then
          x(i)=3.0d0
        elseif(mod(i,4).eq.2) then
          x(i)=-1.0d0
        elseif(mod(i,4).eq.3) then
          x(i)=0.0d0
        else
          x(i)=1.0d0
        endif
   31 continue
      nb=2*(n-2)
      return
   40 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 41 i=1,n
        x(i)=2.0d0
   41 continue
      x(1)=1.0d0
      nb=5*(n-2)/2
      xmax=1.0d1
      return
   50 if (n.lt.3) go to 999
      do 51 i=1,n
        x(i)=-1.0d0
   51 continue
      nb=n
      return
   60 if (n.lt.6) go to 999
      do 61 i=1,n
        x(i)=-1.0d0
   61 continue
      nb=n
      return
   70 if (n.lt.2) go to 999
      do 71 i=1,n-1
        x(i)=0.5d0
   71 continue
      x(n)=-2.0d0
      nb=2*(n-1)
      return
   80 if (n.lt.4) go to 999
      n=n-mod(n,4)
      do 81 i=1,n
        x(i)=sin(dble(i))**2
   81 continue
      nb=5*n
      return
   90 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 91 i=1,n
        x(i)=5.0d0
   91 continue
      nb=3*(n-2)
      return
  100 if (n.lt.2) go to 999
      do 101 i=1,n
        x(i)=0.2d0
  101 continue
      nb=2*(n-1)
      return
  110 continue
      if (n.lt.2) go to 999
      n=n-mod(n,2)
      do 111 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-0.8d0
        else
          x(i)=-0.8d0
        endif
  111 continue
      nb=2*(n-1)
      return
  120 continue
      if (n.lt.5) go to 999
      if (mod(n-5,3).ne.0) n=n-mod(n-5,3)
      do 121 i=1,n
        x(i)=-1.0d0
  121 continue
      nb=6*((n-5)/3+1)
      return
  130 continue
      if (n.lt.5) go to 999
      if (mod(n-5,3).ne.0) n=n-mod(n-5,3)
      do 131 i=1,n
        x(i)=-1.0d0
  131 continue
      nb=7*((n-5)/3+1)
      return
  140 continue
      if (n.lt.4) go to 999
      do 141 i=1,n
        if (mod(i,4).eq.1) then
          x(i)=-0.8d0
        else if (mod(i,4).eq.2) then
          x(i)= 1.2d0
        else if (mod(i,4).eq.3) then
          x(i)=-1.2d0
        else
          x(i)= 0.8d0
        endif
  141 continue
      y(1)=14.4d0
      y(2)=6.8d0
      y(3)=4.2d0
      y(4)=3.2d0
  142 if (mod(n-4,2).ne.0) n=n-mod(n-4,2)
      nb=4*((n-4)/2+1)
      return
  150 continue
      if (n.lt.4) go to 999
      do 151 i=1,n
        if (mod(i,4).eq.1) then
          x(i)=-0.8d0
        else if (mod(i,4).eq.2) then
          x(i)= 1.2d0
        else if (mod(i,4).eq.3) then
          x(i)=-1.2d0
        else
          x(i)= 0.8d0
        endif
  151 continue
      y(1)=35.8d0
      y(2)=11.2d0
      y(3)=6.2d0
      y(4)=4.4d0
      go to 142
  160 continue
      if (n.lt.4) go to 999
      do 161 i=1,n
        if (mod(i,4).eq.1) then
          x(i)=-0.8d0
        else if (mod(i,4).eq.2) then
          x(i)= 1.2d0
        else if (mod(i,4).eq.3) then
          x(i)=-1.2d0
        else
          x(i)= 0.8d0
        endif
  161 continue
      y(1)=30.6d0
      y(2)=72.2d0
      y(3)=124.4d0
      y(4)=187.4d0
      go to 142
  170 if (n.lt.4) go to 999
      n=n-mod(n,2)
      nb=n
      do 171 i=1,n
        if (mod(i,8).eq.1) x(i)=1.0d-1
        if (mod(i,8).eq.2.or.mod(i,8).eq.0) x(i)=2.0d-1
        if (mod(i,8).eq.3.or.mod(i,8).eq.7) x(i)=3.0d-1
        if (mod(i,8).eq.4.or.mod(i,8).eq.6) x(i)=4.0d-1
        if (mod(i,8).eq.5) x(i)=5.0d-1
  171 continue
      return
  180 if (n.lt.3) go to 999
      nb=n
      do 181 i=1,n
        x(i)=1.2d1
  181 continue
      return
  190 if (n.lt.7) go to 999
      nb=n
      do 191 i=1,n
        x(i)=-1.0d0
  191 continue
      return
  200 if (n.lt.3) go to 999
      nb=n
      do 201 i=1,n
        x(i)=dble(i)/dble(n+1)
        x(i)=x(i)*(x(i)-1.0d0)
  201 continue
      return
  210 continue
      if (n.lt.5) go to 999
      if (mod(n-5,3).ne.0) n=n-mod(n-5,3)
      do 211 i=1,n
        x(i)=-1.0d0
  211 continue
      nb=7*((n-5)/3+1)
      return
  220 if (n.lt.3) go to 999
      do 221 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-1.2d0
        else
          x(i)= 1.0d0
        endif
  221 continue
      nb=2*(n-1)
      return
  260 if (n.lt.4) go to 999
      n=n-mod(n,4)
      do 261 i=1,n
        x(i)=1.0d0
  261 continue
      nb=5*n
      return
  270 if (n.lt.7) go to 999
      do 271 i=1,n
        x(i)=5.0d0
  271 continue
      y(1)=sin(1.0d0)
      nb=13*(n-6)
      return
  999 ierr=1
      return
      end subroutine tiud15

! subroutine tafu15             all systems                92/12/01
! portability : all systems
! 92/12/01 ra : original version
!
! purpose :
!  values of partial functions in the sum of squares.
!
! parameters :
!  ii  n  number of variables.
!  ii  ka  index of the given partial function.
!  ri  x(n)  vector of variables.
!  ro  fa  value of the ka-th partial function at the point x.
!  ii  next  number of the selected test problem.
!
      subroutine tafu15(n,ka,x,fa,next)
      integer n,ka,next
      double precision x(*),fa
      double precision a,c,d,p,alfa,u,v
      integer i,j,k,l,m,ia,ib,ic
      double precision y(20)
      common /empr15/ y
      go to(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160, &
       170,180,190,200,210,220),next
   10 i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      fa=1.0d1*(x(i)**2-x(i+1))
      else
      fa=x(i)-1.0d0
      end if
      return
   20 a=sqrt(10.0d0)
      d=sqrt(90.0d0)
      i=2*((ka+5)/6)
      if (mod(ka,6).eq.1) then
      fa=1.0d1*(x(i-1)**2-x(i))
      else if (mod(ka,6).eq.2) then
      fa=x(i-1)-1.0d0
      else if (mod(ka,6).eq.3) then
      fa=d*(x(i+1)**2-x(i+2))
      else if (mod(ka,6).eq.4) then
      fa=x(i+1)-1.0d0
      else if (mod(ka,6).eq.5) then
      fa=a*(x(i)+x(i+2)-2.0d0)
      else
      fa=(x(i)-x(i+2))/a
      end if
      return
   30 a=sqrt(1.0d1)
      c=sqrt(5.0d0)
      i=2*((ka+3)/4)
      if (mod(ka,4).eq.1) then
      fa=x(i-1)+1.0d1*x(i)
      else if (mod(ka,4).eq.2) then
      fa=c*(x(i+1)-x(i+2))
      else if (mod(ka,4).eq.3) then
      fa=(x(i)-2.0d0*x(i+1))**2
      else
      fa=a*(x(i-1)-x(i+2))**2
      end if
      return
   40 i=2*((ka+4)/5)
      if (mod(ka,5).eq.1) then
      fa=(exp(x(i-1))-x(i))**2
      else if (mod(ka,5).eq.2) then
      fa=1.0d1*(x(i)-x(i+1))**3
      else if (mod(ka,5).eq.3) then
      p=x(i+1)-x(i+2)
      fa=(sin(p)/cos(p))**2
      else if (mod(ka,5).eq.4) then
      fa=x(i-1)**4
      else
      fa=x(i+2)-1.0d0
      end if
      return
   50 i=ka
      fa=(3.0d0-2.0d0*x(i))*x(i)+1.0d0
      if (i.gt.1) fa=fa-x(i-1)
      if (i.lt.n) fa=fa-x(i+1)
      return
   60 i=ka
      fa=(2.0d0+5.0d0*x(i)**2)*x(i)+1.0d0
      do 61 j=max(1,i-5),min(n,i+1)
      fa=fa+x(j)*(1.0d0+x(j))
   61 continue
      return
   70 i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      fa=x(i)+x(i+1)*((5.0d0-x(i+1))*x(i+1)-2.0d0)-1.3d1
      else
      fa=x(i)+x(i+1)*((1.0d0+x(i+1))*x(i+1)-1.4d1)-2.9d1
      end if
      return
   80 i=mod(ka,n/2)+1
      j=i+n/2
      m=5*n
      if (ka.le.m/2) then
      ia=1
      else
      ia=2
      end if
      ib=5-ka/(m/4)
      ic=mod(ka,5)+1
      fa=(x(i)**ia-x(j)**ib)**ic
      return
   90 i=2*((ka+5)/6)-1
      if (mod(ka,6).eq.1) then
      fa=x(i)+3.0d0*x(i+1)*(x(i+2)-1.0d0)+x(i+3)**2-1.0d0
      else if (mod(ka,6).eq.2) then
      fa=(x(i)+x(i+1))**2+(x(i+2)-1.0d0)**2-x(i+3)-3.0d0
      else if (mod(ka,6).eq.3) then
      fa=x(i)*x(i+1)-x(i+2)*x(i+3)
      else if (mod(ka,6).eq.4) then
      fa=2.0d0*x(i)*x(i+2)+x(i+1)*x(i+3)-3.0d0
      else if (mod(ka,6).eq.5) then
      fa=(x(i)+x(i+1)+x(i+2)+x(i+3))**2+(x(i)-1.0d0)**2
      else
      fa=x(i)*x(i+1)*x(i+2)*x(i+3)+(x(i+3)-1.0d0)**2-1.0d0
      end if
      return
  100 i=(ka+1)/2
      j=mod(ka,2)
      if (j.eq.0) then
      fa=6.0d0-exp(2.0d0*x(i))-exp(2.0d0*x(i+1))
      else if (i.eq.1) then
      fa=4.0d0-exp(x(i))-exp(x(i+1))
      else if (i.eq.n) then
      fa=8.0d0-exp(3.0d0*x(i-1))-exp(3.0d0*x(i))
      else
      fa=8.0d0-exp(3.0d0*x(i-1))-exp(3.0d0*x(i))+ &
       4.0d0-exp(x(i))-exp(x(i+1))
      end if
      return
  110 i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      fa=1.0d1*(2.0d0*x(i)/(1.0d0+x(i)**2)-x(i+1))
      else
      fa=x(i)-1.0d0
      end if
      return
  120 i=3*((ka+5)/6)-2
      if (mod(ka,6).eq.1) then
      fa=1.0d1*(x(i)**2-x(i+1))
      else if (mod(ka,6).eq.2) then
      fa=x(i+2)-1.0d0
      else if (mod(ka,6).eq.3) then
      fa=(x(i+3)-1.0d0)**2
      else if (mod(ka,6).eq.4) then
      fa=(x(i+4)-1.0d0)**3
      else if (mod(ka,6).eq.5) then
      fa=x(i)**2*x(i+3)+sin(x(i+3)-x(i+4))-1.0d1
      else
      fa=x(i+1)+(x(i+2)**2*x(i+3))**2-2.0d1
      end if
      return
  130 i=3*((ka+6)/7)-2
      if (mod(ka,7).eq.1) then
      fa=1.0d1*(x(i)**2-x(i+1))
      else if (mod(ka,7).eq.2) then
      fa=1.0d1*(x(i+1)**2-x(i+2))
      else if (mod(ka,7).eq.3) then
      fa=(x(i+2)-x(i+3))**2
      else if (mod(ka,7).eq.4) then
      fa=(x(i+3)-x(i+4))**2
      else if (mod(ka,7).eq.5) then
      fa=x(i)+x(i+1)**2+x(i+2)-3.0d1
      else if (mod(ka,7).eq.6) then
      fa=x(i+1)-x(i+2)**2+x(i+3)-1.0d1
      else
      fa=x(i)*x(i+4)-1.0d1
      end if
      return
  140 i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      fa=-y(l)
      do 142 k=1,3
      a=dble(k*k)/dble(l)
      do 141 j=1,4
      if (x(i+j).eq.0) x(i+j)=1.0d-16
      a=a*sign(1.0d0,x(i+j))*abs(x(i+j))**(dble(j)/dble(k*l))
  141 continue
      fa=fa+a
  142 continue
      return
  150 i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      fa=-y(l)
      do 152 k=1,3
      a=0.0d0
      do 151 j=1,4
      a=a+x(i+j)*(dble(j)/dble(k*l))
  151 continue
      fa=fa+exp(a)*dble(k*k)/dble(l)
  152 continue
      return
  160 i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      fa=-y(l)
      do 161 j=1,4
      fa=fa+dble((1-2*mod(j,2))*l*j*j)*sin(x(i+j))+ &
       dble(l*l*j)*cos(x(i+j))
  161 continue
      return
  170 alfa=0.5d0
      if (ka.eq.1) then
      fa=alfa-(1.0d0-alfa)*x(3)-x(1)*(1.0d0+4.0d0*x(2))
      else if (ka.eq.2) then
      fa=-(2.0d0-alfa)*x(4)-x(2)*(1.0d0+4.0d0*x(1))
      else if (ka.eq.n-1) then
      fa=alfa*x(n-3)-x(n-1)*(1.0d0+4.0d0*x(n))
      else if (ka.eq.n) then
      fa=alfa*x(n-2)-(2.0d0-alfa)-x(n)*(1.0d0+4.0d0*x(n-1))
      else if (mod(ka,2).eq.1) then
      fa=alfa*x(ka-2)-(1.0d0-alfa)*x(ka+2)- &
       x(ka)*(1.0d0+4.0d0*x(ka+1))
      else
      fa=alfa*x(ka-2)-(2.0d0-alfa)*x(ka+2)- &
       x(ka)*(1.0d0+4.0d0*x(ka-1))
      end if
      return
  180 if (ka.lt.2) then
      fa=4.0d0*(x(ka)-x(ka+1)**2)
      else if (ka.lt.n) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)
      else
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))
      end if
      return
  190 if (ka.eq.1) then
      fa=-2.0d0*x(ka)**2+3.0d0*x(ka)-2.0d0*x(ka+1)+ &
       3.0d0*x(n-4)-x(n-3)-x(n-2)+0.5d0*x(n-1)-x(n)+1.0d0
      else if (ka.le.n-1) then
      fa=-2.0d0*x(ka)**2+3.0d0*x(ka)-x(ka-1)-2.0d0*x(ka+1)+ &
       3.0d0*x(n-4)-x(n-3)-x(n-2)+0.5d0*x(n-1)-x(n)+1.0d0
      else
      fa=-2.0d0*x(n)**2+3.0d0*x(n)-x(n-1)+ &
       3.0d0*x(n-4)-x(n-3)-x(n-2)+0.5d0*x(n-1)-x(n)+1.0d0
      end if
      return
  200 u=1.0d0/dble(n+1)
      v=dble(ka)*u
      fa=2.0d0*x(ka)+0.5d0*u*u*(x(ka)+v+1.0d0)**3+1.0d0
      if (ka.gt.1) fa=fa-x(ka-1)
      if (ka.lt.n) fa=fa-x(ka+1)
      return
  210 i=3*((ka+6)/7)-2
      if (mod(ka,7).eq.1) then
      fa=1.0d1*(x(i)**2-x(i+1))
      else if (mod(ka,7).eq.2) then
      fa=x(i+1)+x(i+2)-2.0d0
      else if (mod(ka,7).eq.3) then
      fa=x(i+3)-1.0d0
      else if (mod(ka,7).eq.4) then
      fa=x(i+4)-1.0d0
      else if (mod(ka,7).eq.5) then
      fa=x(i)+3.0d0*x(i+1)
      else if (mod(ka,7).eq.6) then
      fa=x(i+2)+x(i+3)-2.0d0*x(i+4)
      else
      fa=1.0d1*(x(i+1)**2-x(i+4))
      end if
      return
  220 i=ka/2
      if (ka.eq.1) then
      fa=x(ka)-1.0d0
      else if (mod(ka,2).eq.0) then
      fa=1.0d1*(x(i)**2-x(i+1))
      else
      fa=2.0d0*exp(-(x(i)-x(i+1))**2)+ &
       exp(-2.0d0*(x(i+1)-x(i+2))**2)
      end if
      return
      end subroutine tafu15

! subroutine tagu15                all systems                92/12/01
! portability : all systems
! 92/12/01 ra : original version
!
! purpose :
!  gradients of partial functions in the sum of squares.
!
! parameters :
!  ii  n  number of variables.
!  ii  ka  index of the given partial function.
!  ri  x(n)  vector of variables.
!  ro  ga(n)  gradient of the ka-th partial function at the point x.
!  ii  next  number of the selected test problem.
!
      subroutine tagu15(n,ka,x,ga,next)
      integer n,ka,next
      double precision x(*),ga(*)
      double precision a,b,c,d,q,r,e,alfa,u,v
      integer i,j,k,l,m,ia,ib,ic
      double precision y(20)
      common /empr15/ y
      go to(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160, &
       170,180,190,200,210,220,260,270),next
   10 i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      ga(i)=2.0d1*x(i)
      ga(i+1)=-1.0d1
      else
      ga(i)=1.0d0
      end if
      return
   20 i=2*((ka+5)/6)
      a=sqrt(9.0d1)
      b=sqrt(1.0d1)
      if (mod(ka,6).eq.1) then
      ga(i-1)=2.0d1 *x(i-1)
      ga(i)=-1.0d1
      else if (mod(ka,6).eq.2) then
      ga(i-1)=1.0d0
      else if (mod(ka,6).eq.3) then
      ga(i+1)=2.0d0*a*x(i+1)
      ga(i+2)=-a
      else if (mod(ka,6).eq.4) then
      ga(i+1)=1.0d0
      else if (mod(ka,6).eq.5) then
      ga(i)=b
      ga(i+2)=b
      else
      ga(i)=1.0d0/b
      ga(i+2)=-1.0d0/b
      end if
      return
   30 i=2*((ka+3)/4)
      a=sqrt(5.0d0)
      b=sqrt(1.0d1)
      if (mod(ka,4).eq.1) then
      ga(i-1)=1.0d0
      ga(i)=1.0d1
      else if (mod(ka,4).eq.2) then
      ga(i+1)=a
      ga(i+2)=-a
      else if (mod(ka,4).eq.3) then
      c=x(i)-2.0d0*x(i+1)
      ga(i)=2.0d0*c
      ga(i+1)=-4.0d0*c
      else
      c=x(i-1)-x(i+2)
      d=2.0d0*b*c
      ga(i-1)=d
      ga(i+2)=-d
      end if
      return
   40 i=2*((ka+4)/5)
      if (mod(ka,5).eq.1) then
      a=exp(x(i-1))
      b=a-x(i)
      c=2.0d0*b
      ga(i-1)=c*a
      ga(i)=-c
      else if (mod(ka,5).eq.2) then
      a=x(i)-x(i+1)
      b=3.0d1*a**2
      ga(i)=b
      ga(i+1)=-b
      else if (mod(ka,5).eq.3) then
      c=x(i+1)-x(i+2)
      q=sin(c)/cos(c)
      r=cos(c)
      d=2.0d0*q/r**2
      ga(i+1)=d
      ga(i+2)=-d
      else if (mod(ka,5).eq.4) then
      ga(i-1)=4.0d0*x(i-1)**3
      else
      ga(i+2)=1.0d0
      end if
      return
   50 i=ka
      ga(i)=3.0d0-4.0d0*x(i)
      if (i.gt.1) ga(i-1)=-1.0d0
      if (i.lt.n) ga(i+1)=-1.0d0
      return
   60 i=ka
      do 61 j=max(1,i-5),min(n,i+1)
      ga(j)=1.0d0+2.0d0*x(j)
   61 continue
      ga(i)=ga(i)+2.0d0+1.5d1*x(i)**2
      return
   70 i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      ga(i)=1.0d0
      ga(i+1)=1.0d1*x(i+1)-3.0d0*x(i+1)**2-2.0d0
      else
      ga(i)=1.0d0
      ga(i+1)=2.0d0*x(i+1)+3.0d0*x(i+1)**2-1.4d1
      end if
      return
   80 i=mod(ka,n/2)+1
      j=i+n/2
      m=5*n
      if (ka.le.m/2) then
      ia=1
      else
      ia=2
      end if
      ib=5-ka/(m/4)
      ic=mod(ka,5)+1
      a=dble(ia)
      b=dble(ib)
      c=dble(ic)
      d=x(i)**ia-x(j)**ib
      if (d.eq.0.0d0) then
      ga(i)=0.0d0
      ga(j)=0.0d0
      else
      e=c*d**(ic-1)
      if (x(i).eq.0.0d0.and.ia.le.1) then
      ga(i)=0.0d0
      else
      ga(i)= e*a*x(i)**(ia-1)
      end if
      if (x(j).eq.0.0d0.and.ib.le.1) then
      ga(j)=0.0d0
      else
      ga(j)=-e*b*x(j)**(ib-1)
      end if
      end if
      return
   90 i=2*((ka+5)/6)-1
      if (mod(ka,6).eq.1) then
      ga(i)=1.0d0
      ga(i+1)=3.0d0*(x(i+2)-1.0d0)
      ga(i+2)=3.0d0*x(i+1)
      ga(i+3)=2.0d0*x(i+3)
      else if (mod(ka,6).eq.2) then
      ga(i)=2.0d0*(x(i)+x(i+1))
      ga(i+1)=2.0d0*(x(i)+x(i+1))
      ga(i+2)=2.0d0*(x(i+2)-1.0d0)
      ga(i+3)=-1.0d0
      else if (mod(ka,6).eq.3) then
      ga(i)=x(i+1)
      ga(i+1)=x(i)
      ga(i+2)=-x(i+3)
      ga(i+3)=-x(i+2)
      else if (mod(ka,6).eq.4) then
      ga(i)=2.0d0*x(i+2)
      ga(i+1)=x(i+3)
      ga(i+2)=2.0d0*x(i)
      ga(i+3)=x(i+1)
      else if (mod(ka,6).eq.5) then
      ga(i)=2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))+2.0d0*(x(i)-1.0d0)
      ga(i+1)=2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))
      ga(i+2)=2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))
      ga(i+3)=2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))
      else
      ga(i)=x(i+1)*x(i+2)*x(i+3)
      ga(i+1)=x(i)*x(i+2)*x(i+3)
      ga(i+2)=x(i)*x(i+1)*x(i+3)
      ga(i+3)=x(i)*x(i+1)*x(i+2)+2.0d0*(x(i+3)-1.0d0)
      end if
      return
  100 if (n.ge.2) then
      i=(ka+1)/2
      j=mod(ka,2)
      if (j.eq.0) then
      ga(i)=-2.0d0*exp(2.0d0*x(i))
      ga(i+1)=-2.0d0*exp(2.0d0*x(i+1))
      else if (i.eq.1) then
      ga(i)=-exp(x(i))
      ga(i+1)=-exp(x(i+1))
      else if (i.eq.n) then
      ga(i-1)=-3.0d0*exp(3.0d0*x(i-1))
      ga(i)=-3.0d0*exp(3.0d0*x(i))
      else
      ga(i-1)=-3.0d0*exp(3.0d0*x(i-1))
      ga(i)=-3.0d0*exp(3.0d0*x(i))-exp(x(i))
      ga(i+1)=-exp(x(i+1))
      end if
      end if
      return
  110 i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      ga(i)=2.0d1*(1.0d0-x(i)**2)/(1.0d0+x(i)**2)**2
      ga(i+1)=-1.0d1
      else
      ga(i)=1.0d0
      end if
      return
  120 i=3*((ka+5)/6)-2
      if (mod(ka,6).eq.1) then
      ga(i)=2.0d1*x(i)
      ga(i+1)=-1.0d1
      else if (mod(ka,6).eq.2) then
      ga(i+2)=1.0d0
      else if (mod(ka,6).eq.3) then
      ga(i+3)=2.0d0*(x(i+3)-1.0d0)
      else if (mod(ka,6).eq.4) then
      ga(i+4)=3.0d0*(x(i+4)-1.0d0)**2
      else if (mod(ka,6).eq.5) then
      ga(i)=2.0d0*x(i)*x(i+3)
      ga(i+3)=x(i)**2+cos(x(i+3)-x(i+4))
      ga(i+4)=-cos(x(i+3)-x(i+4))
      else
      ga(i+1)=1.0d0
      ga(i+2)=4.0d0*x(i+2)*(x(i+2)*x(i+3))**2
      ga(i+3)=2.0d0*x(i+2)**4*x(i+3)
      end if
      return
  130 i=3*((ka+6)/7)-2
      if (mod(ka,7).eq.1) then
      ga(i)=2.0d1*x(i)
      ga(i+1)=-1.0d1
      else if (mod(ka,7).eq.2) then
      ga(i+1)=2.0d1*x(i+1)
      ga(i+2)=-1.0d1
      else if (mod(ka,7).eq.3) then
      ga(i+2)= 2.0d0*(x(i+2)-x(i+3))
      ga(i+3)=-2.0d0*(x(i+2)-x(i+3))
      else if (mod(ka,7).eq.4) then
      ga(i+3)= 2.0d0*(x(i+3)-x(i+4))
      ga(i+4)=-2.0d0*(x(i+3)-x(i+4))
      else if (mod(ka,7).eq.5) then
      ga(i)=1.0d0
      ga(i+1)=2.0d0*x(i+1)
      ga(i+2)=1.0d0
      else if (mod(ka,7).eq.6) then
      ga(i+1)=1.0d0
      ga(i+2)=-2.0d0*x(i+2)
      ga(i+3)=1.0d0
      else
      ga(i)=x(i+4)
      ga(i+4)=x(i)
      end if
      return
  140 i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      do 141 j=1,4
      ga(i+j)=0.0d0
  141 continue
      do 144 k=1,3
      a=dble(k*k)/dble(l)
      do 142 j=1,4
      a=a*sign(1.0d0,x(i+j))*abs(x(i+j))**(dble(j)/dble(k*l))
  142 continue
      do 143 j=1,4
      ga(i+j)=ga(i+j)+(dble(j)/dble(k*l))*a/x(i+j)
  143 continue
  144 continue
      return
  150 i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      do 151 j=1,4
      ga(i+j)=0.0d0
  151 continue
      do 154 k=1,3
      a=0.0d0
      do 152 j=1,4
      a=a+x(i+j)*(dble(j)/dble(k*l))
  152 continue
      a=exp(a)*dble(k*k)/dble(l)
      do 153 j=1,4
      ga(i+j)=ga(i+j)+a*(dble(j)/dble(k*l))
  153 continue
  154 continue
      return
  160 i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      do 161 j=1,4
      ga(i+j)=dble((1-2*mod(j,2))*l*j*j)*cos(x(i+j))- &
       dble(l*l*j)*sin(x(i+j))
  161 continue
      return
  170 alfa=0.5d0
      if (ka.eq.1) then
      ga(1)=-1.0d0-4.0d0*x(2)
      ga(2)=-4.0d0*x(1)
      ga(3)=alfa-1.0d0
      else if (ka.eq.2) then
      ga(1)=-4.0d0*x(2)
      ga(2)=-1.0d0-4.0d0*x(1)
      ga(4)=-2.0d0+alfa
      else if (ka.eq.n-1) then
      ga(n-3)=alfa
      ga(n-1)=-1.0d0-4.0d0*x(n)
      ga(n)=-4.0d0*x(n-1)
      else if (ka.eq.n) then
      ga(n-2)=alfa
      ga(n-1)=-4.0d0*x(n)
      ga(n)=-1.0d0-4.0d0*x(n-1)
      else if (mod(ka,2).eq.1) then
      ga(ka-2)=alfa
      ga(ka)=-1.0d0-4.0d0*x(ka+1)
      ga(ka+1)=-4.0d0*x(ka)
      ga(ka+2)=-1.0d0+alfa
      else
      ga(ka-2)=alfa
      ga(ka-1)=-4.0d0*x(ka)
      ga(ka)=-1.0d0-4.0d0*x(ka-1)
      ga(ka+2)=-2.0d0+alfa
      end if
      return
  180 if (ka.lt.2) then
      ga(ka)=4.0d0
      ga(ka+1)=-8.0d0*x(ka+1)
      else if (ka.lt.n) then
      ga(ka-1)=-8.0d0*x(ka)
      ga(ka)=24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0
      ga(ka+1)=-8.0d0*x(ka+1)
      else
      ga(ka-1)=-8.0d0*x(ka)
      ga(ka)=24.0d0*x(ka)**2-8.0d0*x(ka-1)+2.0d0
      end if
      return
  190 if (ka.eq.1) then
      ga(n-4)=3.0d0
      ga(n-3)=-1.0d0
      ga(n-2)=-1.0d0
      ga(n-1)=0.50d0
      ga(n)=-1.0d0
      ga(ka)=-4.0d0*x(ka)+3.0d0
      ga(ka+1)=-2.0d0
      else if (ka.le.n-1) then
      ga(ka-1)=0.0d0
      ga(ka)=0.0d0
      ga(ka+1)=0.0d0
      ga(n-4)=3.0d0
      ga(n-3)=-1.0d0
      ga(n-2)=-1.0d0
      ga(n-1)=0.50d0
      ga(n)=-1.0d0
      ga(ka-1)=ga(ka-1)-1.0d0
      ga(ka)=ga(ka)-4.0d0*x(ka)+3.0d0
      ga(ka+1)=ga(ka+1)-2.0d0
      else
      ga(n-4)=3.0d0
      ga(n-3)=-1.0d0
      ga(n-2)=-1.0d0
      ga(n-1)=0.50d0
      ga(n)=-4.0d0*x(n)+2.0d0
      end if
      return
  200 u=1.0d0/dble(n+1)
      v=dble(ka)*u
      ga(ka)=2.0d0+1.5d0*u**2*(x(ka)+v+1.0d0)**2
      if (ka.gt.1) ga(ka-1)=-1.0d0
      if (ka.lt.n) ga(ka+1)=-1.0d0
      return
  210 i=3*((ka+6)/7)-2
      if (mod(ka,7).eq.1) then
      ga(i)=2.0d1*x(i)
      ga(i+1)=-1.0d1
      else if (mod(ka,7).eq.2) then
      ga(i+1)=1.0d0
      ga(i+2)=1.0d0
      else if (mod(ka,7).eq.3) then
      ga(i+3)=1.0d0
      else if (mod(ka,7).eq.4) then
      ga(i+4)=1.0d0
      else if (mod(ka,7).eq.5) then
      ga(i)=1.0d0
      ga(i+1)=3.0d0
      else if (mod(ka,7).eq.6) then
      ga(i+2)=1.0d0
      ga(i+3)=1.0d0
      ga(i+4)=-2.0d0
      else
      ga(i+1)=2.0d1*x(i+1)
      ga(i+4)=-1.0d1
      end if
      return
  220 i=ka/2
      if (ka.eq.1) then
      ga(ka)=1.0d0
      else if (mod(ka,2).eq.0) then
      ga(i)=2.0d1*x(i)
      ga(i+1)=-1.0d1
      else
      a=2.0d0*exp(-(x(i)-x(i+1))**2)
      b=exp(-2.0d0*(x(i+1)-x(i+2))**2)
      ga(i)= -2.0d0*(x(i)-x(i+1))*a
      ga(i+1)=2.0d0*(x(i)-x(i+1))*a-4.0d0*(x(i+1)-x(i+2))*b
      ga(i+2)=4.0d0*(x(i+1)-x(i+2))*b
      end if
      return
  260 i=mod(ka,n/2)+1
      j=i+n/2
      m=5*n
      ia=ka/(m/4)+1
      a=dble(ia)
      b=dble(ka/(m/5)+1)
      ga(i)=a*x(i)**(ia-1)*exp(b*x(j))
      ga(j)=b*x(i)**ia*exp(b*x(j))+1.0d0
      return
  270 ia=min(max(mod(ka,13)-2,1),7)
      ib=(ka+12)/13
      i=ia+ib-1
      if (ia.eq.7) then
      j=ib
      else
      j=ia+ib
      end if
      c=3.0d0*dble(ia)/1.0d1
      a=cos(c)
      b=exp(sin(c*x(j)))
      ga(i)=-cos(x(i))*(1.0d0+a)+5.0d0*b
      ga(j)=(1.0d0+a)+5.0d0*(x(i)-2.0d0)*c*cos(c*x(j))*b
      do 271 l=0,6
      if (ib+l.ne.i.and.ib+l.ne.j) ga(ib+l)=0.5d0*cos(x(ib+l))
  271 continue
      return
      end subroutine tagu15

! subroutine tiud16             all systems                 99/12/01
! portability : all systems
! 99/12/01 ga : original version
!
! purpose :
!  initial values of the variables for nonlinear equations.
!  unconstrained and dense version.
!
! parameters :
!  io  n  number of variables.
!  io  na  number of equations.
!  ro  x(n)  vector of variables.
!  ro  fmin  lower bound for value of the objective function.
!  ro  xmax  maximum stepsize.
!  io  next  number of the test problem.
!  io  ierr  error indicator.
!
      subroutine tiud16(n,na,x,fmin,xmax,next,ierr)
      integer n,na,next,ierr
      double precision x(n),fmin,xmax
      double precision y(1000),f,s,t
      integer i,n1,alf,bet
      na=n
      fmin=0.0d0
      xmax=1.0d3
      ierr=0
      go to (201,202,203,204,205,206,207,208,209,210,211,212,213,214, &
       217,218,219,220,221,222,223,224,225,226,227,228, &
       229,230,231,231),next
201   n1=n-1
      do 3030 i=1,n1
        x(i)=-1.2d0
3030  continue
      x(n)=-1.0d0
      return
202   do 3040 i=1,n
        x(i)=2.0d0
3040  continue
      return
203   do 3050 i=1,n
        x(i)=0.5d0
3050  continue
      return
204   n1=n-1
      do 3060 i=1,n1,2
        x(i)=-1.2d0
        x(i+1)=1.0d0
3060  continue
      return
205   do 3070 i=1,n
        x(i)=1.5d0
3070  continue
      return
206   do 3080 i=1,n
        x(i)=0.0d0
3080  continue
      return
207   do 3090 i=1,n
        x(i)=-1.0d0
3090  continue
      return
208   do 4000 i=1,n
        x(i)=-1.0d0
4000  continue
      return
209   do 4010 i=1,n
        x(i)=0.5d0
4010  continue
      return
210   do 4020 i=1,n
        x(i)=0.0d0
4020  continue
      return
211   do 4030 i=1,n
        x(i)=0.0d0
4030  continue
      return
212   do 4040 i=1,n
        x(i)=0.5d0
4040  continue
      return
213   do 4050 i=1,n
        x(i)=1.0d0
4050  continue
      return
214   do 4060 i=1,n
        x(i)=-1.0d0
4060  continue
      return
215   do 4070 i=1,n
        x(i)=0.0d0
4070  continue
      return
216   do 4080 i=1,n
         x(i)=dble(i)/dble(n+1)
4080  continue
      return
217   do 4090 i=1,n
        x(i)=10.0d0
4090  continue
      return
218   do 5000 i=1,n
        t=dble(i)/dble(n+1)
        x(i)=t*(t-1.0d0)
5000  continue
      return
219   do 5010 i=1,n
        t=dble(i)/dble(n+1)
        x(i)=t*(t-1.d0)
5010  continue
      return
220   do 5020 i=1,n
        x(i)=1.0d0/dble(n)
5020  continue
      return
221   do 5030 i=1,n
        x(i)=1.0d0-dble(i)/dble(n)
5030  continue
      return
222   do 5040 i=1,n
        x(i)=-1.0d0
5040  continue
      return
223   alf=5
      bet=14
      do 5051 i=1,n
        x(i)=0.0d0
5051  continue
      do 5052 i=1,n
        call eafu16(n,i,x,f,next)
        y(i)=-f
5052  continue
      s=dble(bet*n)/dble(bet**2*n**2-(alf+1)**2*(n-1)**2)
      do 5050 i=1,n
        x(i)=s*y(i)
5050  continue
      return
224   do 5060 i=1,n
        x(i)=1.0d0
5060  continue
      return
225   do 5070 i=1,n
        x(i)=0.0d0
5070  continue
      return
226   do 5080 i=1,n
        x(i)=1.0d0
5080   continue
      return
227   do 5090 i=1,n
        x(i)=1.0d0
5090   continue
      return
228   do 6000 i=1,n
        x(i)=1.0d0
6000   continue
      return
229   do 6010 i=1,n
        x(i)=1.0d0
6010  continue
      return
230   t=dble(2)/dble(n+2)
      n1=n/2
      do 6020 i=1,n1
        s=dble(i)*t
        x(i)=s*(1.0d0-s)
        x(n1+i)=x(i)
6020  continue
      return
231   n1=int(sqrt(dble(n)))
      n=n1*n1
      do 6030 i=1,n
        x(i)=1.0d0
6030  continue
      return
      end subroutine tiud16

! subroutine eafu16             all systems                 99/12/01
! portability : all systems
! 99/12/01 ga : original version
!
! purpose :
!  values of test functions for nonlinear equations.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ii  i  index of the approximated function.
!  ri  x(n)  vector of variables.
!  ro  f  value of the approximated function at the
!          selected point.
!  ii  next  number of the test problem.
!
      subroutine eafu16(n,i,x,f,next)
      integer n,i,next
      double precision x(n),f
      double precision s,t,s1,s2,s3,t1,c,h,al1,al2,be1,be2
      double precision al,be,a,b,ga,ca,cb,ff,u,la,h2
      integer j,n1,i1,i2,j1,j2,k,alf,bet,gam,l,nd
      go to (201,202,203,204,205,206,207,208,209,210,211,212,213,214, &
       217,218,219,220,221,222,223,224,225,226,227,228, &
       229,230,231,232),next
201   if(i.eq.1) then
        f=1.0d0-x(1)
      else
        f=10.0d0*dble(i-1)*(x(i)-x(i-1))**2
      endif
      return
202   if(i.eq.n) then
        f=x(i)-0.1d0*x(1)**2
      else
        f=x(i)-0.1d0*x(i+1)**2
      endif
      return
203   t=-1.0d0
      if(i.lt.n) then
        t=t+x(i)
        do 3010 j=1,n
           t=t+x(j)
3010    continue
        f=t-dble(n)
      else
        s=1.0d0
        do 3020 j=1,n
           s=s*x(j)
3020    continue
        f=t+s
      endif
      return
204   if(i/2*2.lt.i) then
        f=1.d0-x(i)
      else
        f=10.0d0*(x(i)-x(i-1)**2)
      endif
      return
205   s=0.0d0
      do 3030 j=1,n
        s=s+x(j)**3
3030  continue
      f=x(i)-1.0d0/dble(2*n)*(s+dble(i))
      return
206   s=1.0d0/dble(n+1)
      if(n.eq.1) then
        f=-2.0d0*x(i)-s**2*exp(x(i))
      else if(i.eq.1) then
        f=-2.0d0*x(i)+x(i+1)-s**2*exp(x(i))
      else if(i.eq.n) then
        f=x(i-1)-2.0d0*x(i)-s**2*exp(x(i))
      else
        f=x(i-1)-2.0d0*x(i)+x(i+1)-s**2*exp(x(i))
      endif
      return
207   s=0.1d0
      if(n.eq.1) then
        f=(3.0d0-s*x(i))*x(i)+1.0d0
      else if(i.eq.1) then
        f=(3.0d0-s*x(i))*x(i)+1.0d0-2.0d0*x(i+1)
      else if(i.eq.n) then
        f=(3.0d0-s*x(i))*x(i)+1.0d0-x(i-1)
      else
        f=(3.0d0-s*x(i))*x(i)+1.0d0-x(i-1)-2.0d0*x(i+1)
      endif
      return
208   s1=1.0d0
      s2=1.0d0
      s3=1.0d0
      j1=3
      j2=3
      if(i-j1.gt.1) then
        i1=i-j1
      else
        i1=1
      endif
      if(i+j2.lt.n) then
        i2=i+j2
      else
        i2=n
      endif
      s=0.0d0
      do 3040 j=i1,i2
        if(j.ne.i) s=s+x(j)+x(j)**2
3040  continue
      f=(s1+s2*x(i)**2)*x(i)+1.d0-s3*s
      return
209   if(i.eq.1) then
        f=x(1)**2-1.0d0
      else
        f=x(i-1)**2+log(x(i))-1.0d0
      endif
      return
210   s=0.0d0
      do 3050 j=1,n
        s=s+x(j)
3050  continue
      f=exp(cos(dble(i)*s))
      return
211   s=0.0d0
      do 3060 j=1,n
        s=s+x(j)**3
3060  continue
      f=1.0d0/dble(2*n)*(s+dble(i))
      return
212   if(i.eq.1) then
        f=x(1)
      else
        f=cos(x(i-1))+x(i)-1.0d0
      endif
      return
213   s=(1.0d0/dble(n+1))**2
      if(n.eq.1) then
        f=2.0d0*x(1)-1.d0+s*(x(1)+sin(x(1)))
      else if(i.eq.1) then
        f=2.0d0*x(1)-x(i+1)+s*(x(1)+sin(x(1)))
      else if(i.eq.n) then
        f=-x(i-1)+2.0d0*x(i)-1.d0+s*(x(i)+sin(x(i)))
      else
        f=-x(i-1)+2.0d0*x(i)-x(i+1)+s*(x(i)+sin(x(i)))
      endif
      return
214   if(i-5.gt.1) then
        i1=i-5
      else
        i1=1
      endif
      if(i+1.lt.n) then
        i2=i+1
      else
        i2=n
      endif
      s=0.d0
      do 3070 j=i1,i2
        if(j.ne.i) s=s+x(j)*(1.0d0+x(j))
3070  continue
      f=x(i)*(2.0d0+5.0d0*x(i)**2)+1.0d0-s
      return
215   s3=0.0d0
      do 4000 k=1,29
        s1=0.0d0
        s=dble(k)/dble(29)
        do 3080 j=1,n
           s1=s1+s**(j-1)*x(j)
3080    continue
        t1=0.0d0
        do 3090 j=2,n
           t1=t1+(j-1)*s**(j-2)*x(j)
3090    continue
        s3=s3+s**(i-2)*((i-1.0d0-2.0d0*s*s1)*(t1-s1**2-1.0d0))
4000  continue
      if (i.eq.1) then
         f=s3+x(1)*(1.0d0-2*(x(2)-x(1)*x(1)-1.0d0))
      else if(i.eq.2) then
        f=s3+x(2)-x(1)*x(1)-1.0d0
      else
        f=s3
      endif
      return
216   s=0.d0
      do 4020 j=1,n
        s1=2.d0*x(j)-1.d0
        t1=s1
        if(i.gt.1) s2=2.d0*s1**2-1.d0
        if(i.eq.1) then
           s=s+s1
        else if(i.eq.2) then
           s=s+s2
        else
           do 4010 k=3,i
              s3=2*t1*s2-s1
              s1=s2
              s2=s3
4010       continue
           s=s+s3
        endif
4020  continue
      if( i/2*2.eq.i) then
         f=s/dble(n)+1.0d0/dble(i**2-1)
      else
        f=s/dble(n)
      endif
      return
217   if(n.eq.1) then
        f=3.d0*x(1)*(20.d0-2.d0*x(1))+100.d0
      else if(i.eq.1) then
        f=3.d0*x(i)*(x(i+1)-2.d0*x(i))+(x(i+1))**2/4.0d0
      else if(i.eq.n) then
        f=3.d0*x(i)*(20.0d0-2.d0*x(i)+x(i-1))+ &
         (20.0d0-x(i-1))**2/4.0d0
      else
        f=3.d0*x(i)*(x(i+1)-2.d0*x(i)+x(i-1))+ &
         (x(i+1)-x(i-1))**2/4.0d0
      endif
      return
218   s=1.0d0/dble(n+1)
      if(n.eq.1) then
        f=2.d0*x(1)+s**2/2.0d0*(x(1)+s+1.0d0)**3
      else if(i.eq.1) then
       f=2.d0*x(1)-x(i+1)+s**2/2.0d0*(x(1)+dble(i)*s+1.0d0)**3
      else if(i.eq.n) then
        f=2.d0*x(i)-x(i-1)+s**2/2.0d0*(x(i)+dble(i)*s+1.0d0)**3
      else
        f=2.d0*x(i)-x(i-1)-x(i+1)+s**2/2.0d0*(x(i)+&
          dble(i)*s+1.0d0)**3
      endif
      return
219   s1=1.0d0/(dble(n)+1.0d0)
      s2=0.0d0
      s3=0.0d0
      do 4030 j=1,i
        s2=s2+dble(j)*s1*(x(j)+dble(j)*s1+1.0d0)**3
4030  continue
      do 4040 j=i+1,n
        s3=s3+(1.0d0-dble(j)*s1)*(x(j)+dble(j)*s1+1.0d0)**3
4040  continue
      f=x(i)+s1/2.0d0*((1.0d0-dble(i)*s1)*s2+dble(i)*s1*s3)
      return
220   s=0.0d0
      do 4050 j=1,n
        s=s+cos(x(j))
4050  continue
      f=dble(n)-s+dble(i)*(1.0d0-cos(x(i)))-sin(x(i))
      return
221   s=0.0d0
      do 4060 j=1,n
        s=s+dble(j)*(x(j)-1.0d0)
4060  continue
      f=x(i)-1.0d0+dble(i)*s*(1.0d0+2.0d0*s*s)
      return
222   if(i.eq.1.and.n.gt.1) then
        f=(3.d0-2.d0*x(i))*x(i)-2.d0*x(i+1)+1.d0
      else if(i.eq.n.and.n.gt.1) then
        f=(3.d0-2.d0*x(i))*x(i)-x(i-1)+1.d0
      else if(n.gt.1) then
        f=(3.d0-2.d0*x(i))*x(i)-x(i-1)-2.d0*x(i+1)+1.d0
      else
        f=(3.d0-2.d0*x(i))*x(i)+1.d0
      endif
      return
223   alf=5
      bet=14
      gam=3
      f=dble(bet*n)*x(i)+(dble(i)-dble(n)/2.0d0)**gam
      do 4070 j=1,n
        if(j.ne.i) then
        t=sqrt(x(j)**2+dble(i)/dble(j))
        s1=log(t)
        f=f+t*(sin(s1)**alf+cos(s1)**alf)
      endif
4070  continue
      return
224   c=0.5d0
      h=1.0d0/dble(n)
      f=(1.0d0-c*h/4.0d0)
      do 4080 j=1,n
        s=c*h*dble(i)/dble(2*(i+j))
        if(j.eq.n) s=s/2.0d0
        f=f-s*x(j)
4080  continue
      f=-1.0d0+x(i)*f
      return
225   if(i.eq.1) then
        f=3.0d0*x(i)**3+2.0d0*x(i+1)-5.0d0+&
        sin(x(i)-x(i+1))*sin(x(i)+x(i+1))
      else if(i.ne.n) then
        f=3.0d0*x(i)**3+2.0d0*x(i+1)-5.0d0+&
        sin(x(i)-x(i+1))*sin(x(i)+x(i+1))+&
        4.0d0*x(i)-3.0d0-x(i-1)*exp(x(i-1)-x(i))
      else
        f= 4.0d0*x(i)-3.0d0-x(i-1)*exp(x(i-1)-x(i))
      endif
      return
226   t=10.0d0
      h=t/dble(n+1)**2
      s=t*x(i)
      if(i.eq.1) then
        f=2.0d0*x(i)-x(i+1)+h*(exp(s)-exp(-s))/2.d0
      endif
      if(1.lt.i.and.i.lt.n) then
        f=-x(i-1)+2.0d0*x(i)-x(i+1)+h*(exp(s)-exp(-s))/2.d0
      endif
      if(i.eq.n) then
        f=-x(i-1)+2.0d0*x(i)-1.0d0+h*(exp(s)-exp(-s))/2.d0
      endif
      return
227   s=0.5d0
      h=1.0d0/dble(n+1)
      t=h**2/s
      h=2.0d0*h
      if(i.eq.1) then
        f=2.0d0*x(i)-x(i+1)-t*(x(i)**2+x(i+1)/h)
      endif
      if(1.lt.i.and.i.lt.n) then
        f=-x(i-1)+2.0d0*x(i)-x(i+1)-t*(x(i)**2+(x(i+1)-x(i-1))/h)
      endif
      if(i.eq.n) then
        f=-x(i-1)+2.0d0*x(i)-0.5d0-t*(x(i)**2+(0.5d0-x(i-1))/h)
      endif
      return
228   s=0.5d0
      h=1.0d0/dble(n+1)
      t=h**2/s
      t1=2.0d0*h
      al=0.0d0
      be=0.5d0
      s1=0.0d0
      do 4090 j=1,i
        if(j.eq.1) then
           s1=s1+dble(j)*(x(j)**2+(x(j+1)-al)/t1)
        endif
        if(1.lt.j.and.j.lt.n) then
           s1=s1+dble(j)*(x(j)**2+(x(j+1)-x(j-1))/t1)
        endif
        if(j.eq.n) then
           s1=s1+dble(j)*(x(j)**2+(be-x(j-1))/t1)
        endif
4090  continue
      s1=(1.0d0-dble(i)*h)*s1
      if(i.eq.n) go to 5010
      s2=0.0d0
      do 5000 j=i+1,n
        if(j.lt.n) then
           s2=s2+(1.0d0-dble(j)*h)*(x(j)**2+(x(j+1)-x(j-1))/t1)
        else
           s2=s2+(1.0d0-dble(j)*h)*(x(j)**2+(be-x(j-1))/t1)
        endif
5000  continue
      s1=s1+dble(i)*s2
5010  f=x(i)-0.5d0*dble(i)*h-t*s1
      return
229   a=-9.0d-3
      b=1.0d-3
      al=0.0d0
      be=25.0d0
      ga=20.0d0
      ca=0.3d0
      cb=0.3d0
      h=(b-a)/dble(n+1)
      t=a+dble(i)*h
      h=h**2
      s=dble(i)/dble(n+1)
      u=al*(1.0d0-s)+be*s+x(i)
      ff=cb*exp(ga*(u-be))-ca*exp(ga*(al-u))
      if(t.le.0) then
        ff=ff+ca
      else
        ff=ff-cb
      endif
      if(n.eq.1) then
        f=-al+2.0d0*x(i)-be+h*ff
      elseif(i.eq.1) then
        f=-al+2.0d0*x(i)-x(i+1)+h*ff
      elseif(i.lt.n) then
        f=-x(i-1)+2.0d0*x(i)-x(i+1)+h*ff
      else
        f=-x(i-1)+2.0d0*x(i)-be+h*ff
      endif
      return
230   al1=0.0d0
      al2=0.0d0
      be1=0.0d0
      be2=0.0d0
      n1=n/2
      h=1.0d0/dble(n1+1)
      t=dble(i)*h
      if(i.eq.1) then
        s1=2.0d0*x(i)-x(i+1)
        b=al1
      else if(i.eq.n1+1) then
        s1=2.0d0*x(i)-x(i+1)
        b=al2
      else if(i.eq.n1) then
        s1=-x(i-1)+2.0d0*x(i)
        b=be1
      else if(i.eq.n) then
        s1=-x(i-1)+2.0d0*x(i)
        b=be2
      else
        s1=-x(i-1)+2.0d0*x(i)-x(i+1)
        b=0.0d0
      endif
      if(i.le.n1) then
        s2=x(i)**2+x(i)+0.1d0*x(n1+i)**2-1.2d0
      else
        s2=0.2d0*x(i-n1)**2+x(i)**2+2.0d0*x(i)-0.6d0
      endif
      f=s1+h**2*s2-b
      return
231   nd=int(sqrt(dble(n)))
      l=mod(i,nd)
      if(l.eq.0) then
         k=i/nd
         l=nd
      else
         k=int(i/nd)+1
      endif
      la=1.0d0
      h=1.0d0/dble(nd+1)
      h2=la*h*h
      if(l.eq.1.and.k.eq.1) then
         f=4.0d0*x(1)-x(2)-x(nd+1)+h2*exp(x(1))
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         f=4.0d0*x(l)-x(l-1)-x(l+1)-x(l+nd)+h2*exp(x(l))
      endif
      if(l.eq.nd.and.k.eq.1) then
         f=4.0d0*x(nd)-x(nd-1)-x(nd+nd)+h2*exp(x(nd))
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         f=4.0d0*x(i)-x(i+1)-x(i-nd)-x(i+nd)+h2*exp(x(i))
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         f=4.0d0*x(i)-x(i-nd)-x(i-1)-x(i+nd)+h2*exp(x(i))
      endif
      if(l.eq.1.and.k.eq.nd) then
         f=4.0d0*x(i)-x(i+1)-x(i-nd)+h2*exp(x(i))
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         f=4.0d0*x(i)-x(i-1)-x(i+1)-x(i-nd)+h2*exp(x(i))
      endif
      if(l.eq.nd.and.k.eq.nd) then
        f=4.0d0*x(i)-x(i-1)-x(i-nd)+h2*exp(x(i))
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
         f=4.0d0*x(i)-x(i-1)-x(i+1)-x(i-nd)-x(i+nd)+h2*exp(x(i))
      endif
      return
232   nd=int(sqrt(dble(n)))
      l=mod(i,nd)
      if(l.eq.0) then
         k=i/nd
         l=nd
      else
         k=int(i/nd)+1
      endif
      h=1.0d0/dble(nd+1)
      h2=h*h
      if(l.eq.1.and.k.eq.1) then
         f=4.0d0*x(1)-x(2)-x(nd+1)+h2*x(1)**2-24.0d0/(h+1.0d0)**2
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         f=4.0d0*x(l)-x(l-1)-x(l+1)-x(l+nd)+h2*x(l)**2 &
      -12.0d0/(dble(l)*h+1.0d0)**2
      endif
      if(l.eq.nd.and.k.eq.1) then
         f=4.0d0*x(nd)-x(nd-1)-x(nd+nd)+h2*x(nd)**2 &
      -12.0d0/(dble(nd)*h+1.0d0)**2-12.0d0/(h+2.0d0)**2
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         f=4.0d0*x(i)-x(i+1)-x(i-nd)-x(i+nd)+h2*x(i)**2 &
      -12.0d0/(dble(k)*h+1.0d0)**2
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         f=4.0d0*x(i)-x(i-nd)-x(i-1)-x(i+nd)+h2*x(i)**2 &
      -12.0d0/(dble(k)*h+2.0d0)**2
      endif
      if(l.eq.1.and.k.eq.nd) then
         f=4.0d0*x(i)-x(i+1)-x(i-nd)+h2*x(i)**2 &
      -12.0d0/(dble(nd)*h+1.0d0)**2-12.0d0/(h+2.0d0)**2
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         f=4.0d0*x(i)-x(i-1)-x(i+1)-x(i-nd)+h2*x(i)**2 &
      -12.0d0/(dble(l)*h+2.0d0)**2
      endif
      if(l.eq.nd.and.k.eq.nd) then
        f=4.0d0*x(i)-x(i-1)-x(i-nd)+h2*x(i)**2 &
      -24.0d0/(dble(nd)*h+2.0d0)**2
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
         f=4.0d0*x(i)-x(i-1)-x(i+1)-x(i-nd)-x(i+nd)+h2*x(i)**2
      endif
      return
      end subroutine eafu16

! subroutine eagu16             all systems                 99/12/01
! portability : all systems
! 99/12/01 ga : original version
!
! purpose :
!  gradients of test functions for nonlinear equations.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ii  i  index of the approximated function.
!  ri  x(n)  vector of variables.
!  ro  g(n)  gradient of the approximated function at the
!          selected point.
!  ii  next  number of the test problem.
!
      subroutine eagu16(n,i,x,g,next)
      integer n,i,next
      double precision x(n),g(n)
      double precision brow(1000),t,s,s1,s2,s3,t1,t2,c,h,h2,la,a,b,al,be
      double precision ga,ca,cb,u,ff,ti,tim1,tim2,si,sim1,sim2,v
      integer j,k,i1,i2,j1,j2,l,alf,bet,n1,nd
      go to (201,202,203,204,205,206,207,208,209,210,211,212,213,214, &
       217,218,219,220,221,222,223,224,225,226,227,228, &
       229,230,231,232),next
201   if(i.eq.1) then
      g(1)=-1.0d0
      do 3010 j=2,n
        g(j)=0.0d0
3010  continue
      else
      s=20.0d0*dble(i-1)*(x(i)-x(i-1))
      do 3020 j=1,n
        if(j.eq.i) then
           g(j)=s
        else if(j.eq.i-1) then
           g(j)=-s
        else
          g(j)=0.0d0
        endif
3020  continue
      endif
      return
202   if(i.eq.n) then
        g(1)=-0.20d0*x(1)
        n1=n-1
        do 3030 j=2,n1
           g(j)=0.0d0
3030    continue
        g(n)=1.0d0
      else
        do 3050 j=1,n
           if(j.eq.i) then
              g(j)=1.0d0
           else if(j.eq.i+1) then
              g(j)=-0.2d0*x(i+1)
           else
              g(j)=0.0d0
           endif
3050    continue
      endif
      return
203   if(i.eq.n) then
        do 3070 j=1,n
           t=1.0d0
           do 3060 k=1,n
              if(k.ne.j) t=t*x(k)
3060       continue
        g(j)=t
3070    continue
      else
        do 3080 j=1,n
           if(j.eq.i) then
              g(j)=2.0d0
           else
              g(j)=1.0d0
           endif
3080    continue
      endif
      return
204   do 3090 j=1,n
        if(j.eq.i.and.i/2*2.eq.i) then
           g(j)=10.0d0
        else if(j.eq.i) then
          g(j)=-1.0d0
        else if(j.eq.i-1.and.i/2*2.eq.i) then
           g(j)=-20.0d0*x(i-1)
        else
           g(j)=0.0d0
        endif
3090  continue
      return
205   do 4000 j=1,n
        if(j.eq.i) then
           g(j)=1.0d0-3.0d0*x(i)**2/(2.0d0*dble(n))
        else
           g(j)=-3.0d0*x(j)**2/(2.0d0*dble(n))
        endif
4000  continue
      return
206   do 4010 j=1,n
        if(j.eq.i) then
           g(j)=-2.0d0-1.0d0/dble(n+1)**2*exp(x(i))
        else if(j.eq.i-1.or.j.eq.i+1) then
           g(j)=1.0d0
        else
           g(j)=0.0d0
        endif
4010  continue
      return
207   s=0.1d0
      do 4020 j=1,n
        if(j.eq.i-1) then
           g(j)=-1.0d0
        else if(j.eq.i) then
           g(j)=3.0d0-2.0d0*s*x(i)
        else if(j.eq.i+1) then
           g(j)=-2.0d0
        else
           g(j)=0.0d0
        endif
4020  continue
      return
208   s1=1.0d0
      s2=1.0d0
      s3=1.0d0
      j1=3
      j2=3
      if(i-j1.gt.1) then
        i1=i-j1
      else
        i1=1
      endif
      if(i+j2.lt.n) then
        i2=i+j2
      else
        i2=n
      endif
      do 4030 j=1,n
        if(j.eq.i) then
           g(j)=s1+3.0d0*s2*x(i)**2
        else if(j.ge.i1.and.j.le.i2) then
           g(j)=-s3*(1.0d0+2.0d0*x(j))
        else
           g(j)=0.0d0
        endif
4030  continue
      return
209   do 4040 j=1,n
        if(i.eq.1.and.j.eq.1) then
           g(j)=2.0d0*x(1)
        else if(i.eq.1) then
           g(j)=0.0d0
        else if(j.eq.i-1) then
           g(j)=2.0d0*x(i-1)
        else if(j.eq.i) then
           g(j)=1.0d0/x(i)
        else
           g(j)=0.0d0
        endif
4040  continue
      return
210   s1=0.0d0
      do 4050 j=1,n
        s1=s1+x(j)
4050  continue
      do 4060 j=1,n
        s2=dble(i)
        s3=s2*s1
        g(j)=-s2*sin(s3)*exp(cos(s3))
4060  continue
      return
211   do 4070 j=1,n
        g(j)=3.0d0*x(j)**2/(dble(n)*2.0d0)
4070  continue
      return
212   do 4080 j=1,n
        if(j.eq.i) then
           g(j)=1.0d0
        else if(j.eq.i-1) then
           g(j)=-sin(x(i-1))
        else
           g(j)=0.0d0
        endif
4080  continue
      return
213   s=(1.0d0/dble(n+1))**2
      do 4090 j=1,n
        if(j.eq.i-1.or.j.eq.i+1) then
           g(j)=-1.0d0
        else if(j.eq.i) then
           g(j)=2.0d0+s*(1.0d0+cos(x(j)))
        else
           g(j)=0.0d0
        endif
4090  continue
      return
214   if(i-5.gt.1) then
        i1=i-5
      else
        i1=1
      endif
      if(i+1.lt.n) then
        i2=i+1
      else
        i2=n
      endif
      do 5000 j=1,n
        if(j.eq.i) then
           g(j)=2.0d0+15.0d0*x(j)**2
        else if(j.ge.i1.and.j.le.i2) then
           g(j)=-1.0d0-2.0d0*x(j)
        else
           g(j)=0.0d0
        endif
5000  continue
      return
215   do 5010 j=1,n
        s3=0.0d0
        do 5011 k=1,29
           s1=0.0d0
           t1=0.d0
           s=dble(k)/29.0d0
           do 5012 l=1,n
              s1=s1+s**(l-1)*x(l)
5012       continue
           do 5013 l=2,n
              t1=t1+s**(l-2)*(l-1)*x(l)
5013       continue
         s3=s3+s**(i+j-2)*((dble(29*(i-1))/dble(k)-2.0d0*s1)* &
         (dble(29*(j - 1))/dble(k)-2.0d0*s1)- &
         2.0d0*(t1-s1*s1-1.0d0))
         if ((i.eq.1).and.(j.eq.1)) then
             g(j)=s3+6.0d0*x(1)*x(1)-2.0d0*x(2)+3.0d0
         else if (((i.eq.1).and.(j.eq.2)).or.((i.eq.2) &
          .and.(j.eq.1))) then
            g(j)=s3-2.0d0*x(1)
         else if ((i.eq.2).and.(j.eq.2)) then
           g(j)=s3+1.0d0
         else
           g(j)=s3
         endif
5011   continue
5010  continue
      return
216   do 5020 j=1,n
            si=2.d0*x(j)-1.d0
            v=2.d0*si
            tim1=2.d0
            ti=2.d0
            tim2=0.d0
            sim1=1.d0
            do 5021 l=2,i
              ti=4.d0*si+v*tim1-tim2
              tim2=tim1
              tim1=ti
              sim2=sim1
              sim1=si
              si=v*sim1-sim2
5021    continue
5022  g(j)=ti/dble(n)
5020  continue
      return
217   if(i.eq.1) then
        s1=0.0d0
      else
        s1=x(i-1)
      endif
      if(i.eq.n) then
        s2=20.0d0
      else
        s2=x(i+1)
      endif
      do 5030 j=1,n
        if(j.eq.i) then
           g(j)=3.0d0*(s2-4.d0*x(i)+s1)
        else if(j.eq.i-1) then
           g(j)=3.0d0*x(i)-(s2-s1)/2.0d0
        else if(j.eq.i+1) then
           g(j)=3.0d0*x(i)+(s2-s1)/2.0d0
        else
           g(j)=0.0d0
        endif
5030  continue
      return
218   s=1.0d0/dble(n+1)
      do 5040 j=1,n
        if(j.eq.i-1.or.j.eq.i+1) then
           g(j)=-1.0d0
        else if(j.eq.i) then
          g(j)=2.0d0+3.0d0*s*s*(1.0d0+dble(i)*s+x(i))**2/2.0d0
       else
          g(j)=0.0d0
       endif
5040  continue
      return
219   s=1.0d0/dble(n+1)
      s1=3.0d0*s/2.0d0
      t1=dble(i)*s
      do 5050 j=1,n
        t2=dble(j)*s
        if(j.lt.i) then
           g(j)=s1*(1.0d0-t1)*t2*(x(j)+t2+1.0d0)**2
        else if(j.eq.i) then
           g(j)=1.0d0+s1*(1.0d0-t1)*t2*(x(j)+t2+1.0d0)**2
        else
           g(j)=s1*(1.0d0-t2)*t1*(x(j)+t2+1.0d0)**2
        endif
5050  continue
      return
220   do 5060 j=1,n
        s=sin(x(j))
        if(j.eq.i) then
           g(j)=s+dble(i)*s-cos(x(j))
        else
           g(j)=s
        endif
5060  continue
      return
221   s=0.0d0
      do 5070 j=1,n
        s=s+dble(j)*(x(j)-1.0d0)
5070  continue
      do 5080 j=1,n
        if(j.ne.i) then
           g(j)=dble(i*j)*(1.0d0+6.0d0*s*s)
        else
           g(j)=1.0d0+dble(i*i)*(1.0d0+6.0d0*s*s)
        endif
5080  continue
      return
222   do 5090 j=1,n
        if(j.eq.i-1) then
           g(j)=-1.0d0
        else if(j.eq.i)then
           g(j)=3.d0-4.d0*x(i)
        else if(j.eq.i+1) then
           g(j)=-2.0d0
        else
           g(j)=0.0d0
        endif
5090  continue
      return
223   alf=5
      bet=14
      do 6000 j=1,n
        if(j.ne.i) then
        t=sqrt(x(j)**2+dble(i)/dble(j))
        s1=log(t)
        g(j)=x(j)*(sin(s1)**alf+cos(s1)**alf+ &
           alf*sin(s1)**(alf-1)*cos(s1)- &
           alf*sin(s1)*cos(s1)**(alf-1))/t
        else
          g(j)=dble(bet*n)
        endif
6000  continue
      return
224   c=0.5d0
      h=1.0d0/dble(n)
      do 6008 j=1,n
        brow(j)=c*h*dble(i)/dble(2*(i+j))
6008  continue
      brow(n)=0.5d0*brow(n)
      do 6010 j=1,n
        if(i.ne.j) then
           g(j)=-brow(j)*x(i)
        else
           t=1.0d0-c*h/4.0d0
           do 6011 l=1,n
             if(l.eq.i) then
               t=t-2.0d0*brow(i)*x(i)
             else
               t=t-brow(l)*x(l)
             endif
6011       continue
           g(j)=t
        endif
6010  continue
      return
225   do 6020 j=1,n
        if(i.eq.1.and.j.eq.1) then
        g(j)=9.0d0*x(i)**2+cos(x(i)-x(i+1))*sin(x(i)+x(i+1))+ &
           sin(x(i)-x(i+1))*cos(x(i)+x(i+1))
        else if(i.eq.1.and.j.eq.2) then
           g(j)=2.0d0-cos(x(i)-x(i+1))*sin(x(i)+x(i+1))+ &
           sin(x(i)-x(i+1))*cos(x(i)+x(i+1))
        else if(i.lt.n.and.j.eq.i) then
           g(j)=9.0d0*x(i)**2+cos(x(i)-x(i+1))*sin(x(i)+x(i+1))+ &
           sin(x(i)-x(i+1))*cos(x(i)+x(i+1))+ &
           4.0d0+x(i-1)*exp(x(i-1)-x(i))
        else if(i.lt.n.and.j.eq.i+1) then
           g(j)=2.0d0-cos(x(i)-x(i+1))*sin(x(i)+x(i+1))+ &
           sin(x(i)-x(i+1))*cos(x(i)+x(i+1))
        else if(j.eq.i-1) then
           g(j)=-exp(x(i-1)-x(i))-x(i-1)*exp(x(i-1)-x(i))
        else if(i.eq.n.and.j.eq.n) then
           g(j)=4.0d0+x(i-1)*exp(x(i-1)-x(i))
        else
           g(j)=0.0d0
        endif
6020  continue
      return
226   t=10.0d0
      h=(t/dble(n+1))**2
      s=t*x(i)
      do 6030 j=1,n
        if(j.eq.i-1.or.j.eq.i+1) then
           g(j)=-1.0d0
        elseif(j.eq.i) then
           g(j)=2.0d0+h*(exp(s)+exp(-s))/2.d0
        else
           g(j)=0.0d0
        endif
6030  continue
      return
227   s=0.5d0
      h=1.0d0/dble(n+1)
      do 6040 j=1,n
        if(j.eq.i-1) then
           g(j)=-1.0d0+h/(2.0d0*s)
        elseif(j.eq.i+1) then
           g(j)=-1.0d0-h/(2.0d0*s)
        elseif(j.eq.i) then
           g(j)=2.0d0*(1.0d0-h**2*x(i)/s)
        else
           g(j)=0.0d0
        endif
6040  continue
      return
228   s=0.5d0
      h=1.0d0/dble(n+1)
      t1=2.0d0*h
      s1=h**2/s
      s2=1.0d0-dble(i)*h
      do 6050 j=1,i
        brow(j)=dble(j)*s2
6050  continue
      if(i.eq.n) go to 6055
      do 6051 j=i+1,n
        brow(j)=dble(i)*(1.0d0-dble(j)*h)
6051  continue
6055  g(1)=-s1*(brow(1)*2.0d0*x(1)-brow(2)/t1)
      g(n)=-s1*(brow(n-1)/t1+brow(n)*2.0d0*x(n))
      do 6052 j=2,n-1
        g(j)=-s1*((brow(j-1)-brow(j+1))/t1+brow(j)*2.0d0*x(j))
6052  continue
      g(i)=1.0d0+g(i)
      return
229   a=-9.0d-3
      b=1.0d-3
      al=0.0d0
      be=25.0d0
      ga=20.0d0
      ca=0.3d0
      cb=0.3d0
      h=(b-a)/dble(n+1)
      s=dble(i)/dble(n+1)
      h=h**2
      u=al*(1.0d0-s)+be*s+x(i)
      ff=ga*(cb*exp(ga*(u-be))+ca*exp(ga*(al-u)))
      do 6060 j=1,n
        if(j.eq.i-1.or.j.eq.i+1) then
           g(j)=-1.0d0
        elseif(j.eq.i) then
           g(j)=2.0d0+h*ff
        else
           g(j)=0.0d0
        endif
6060  continue
      return
230   n1=n/2
      do 6070 j=1,n
        g(j)=0.0d0
6070  continue
      h=1.0d0/dble(n1+1)**2
      if(i.eq.1) then
        g(i)=2.0d0+h*(2.0d0*x(i)+1.0d0)
        g(i+1)=-1.0d0
        g(n1+1)=h*0.2d0*x(n1+i)
      else if(i.eq.n1+1) then
        g(1)=h*0.4d0*x(1)
        g(i)=2.0d0+h*(2.0d0*x(i)+2.0d0)
        g(i+1)=-1.0d0
      else if(i.eq.n1) then
        g(i-1)=-1.0d0
        g(i)=2.0d0+h*(2.0d0*x(i)+1.0d0)
        g(n1+i)=h*0.2d0*x(n1+i)
      else if(i.eq.n) then
        g(n1)=h*0.4d0*x(n1)
        g(i-1)=-1.0d0
        g(i)=2.0d0+h*(2.0d0*x(i)+2.0d0)
      else if(i.lt.n1) then
        g(i-1)=-1.0d0
        g(i)=2.0d0+h*(2.0d0*x(i)+1.0d0)
        g(i+1)=-1.0d0
        g(n1+i)=h*0.2d0*x(n1+i)
      else
         g(i-n1)=h*0.4d0*x(i-n1)
         g(i-1)=-1.0d0
         g(i)=2.0d0+h*(2.0d0*x(i)+2.0d0)
         g(i+1)=-1.0d0
      endif
      return
231   nd=int(sqrt(dble(n)))
      l=mod(i,nd)
      if(l.eq.0) then
         k=i/nd
         l=nd
      else
         k=int(i/nd)+1
      endif
      la=1.0d0
      h=1.0d0/dble(nd+1)
      h2=la*h*h
      do 2010 j=1,n
        g(j)=0.0d0
2010  continue
      if(l.eq.1.and.k.eq.1) then
         g(1)=4.0d0+h2*exp(x(1))
         g(2)=-1.0d0
         g(nd+1)=-1.0d0
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         g(l)=4.0d0+h2*exp(x(l))
         g(l-1)=-1.0d0
         g(l+1)=-1.0d0
         g(l+nd)=-1.0d0
      endif
      if(l.eq.nd.and.k.eq.1) then
         g(nd)=4.0d0+h2*exp(x(nd))
         g(nd-1)=-1.0d0
         g(nd+nd)=-1.0d0
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         g(i)=4.0d0+h2*exp(x(i))
         g(i-nd)=-1.0d0
         g(i+1)=-1.0d0
         g(i+nd)=-1.0d0
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         g(i)=4.0d0+h2*exp(x(i))
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
         g(i+nd)=-1.0d0
      endif
      if(l.eq.1.and.k.eq.nd) then
         g(i)=4.0d0+h2*exp(x(i))
         g(i-nd)=-1.0d0
         g(i+1)=-1.0d0
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         g(i)=4.0d0+h2*exp(x(i))
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
         g(i+1)=-1.0d0
      endif
      if(l.eq.nd.and.k.eq.nd) then
         g(i)=4.0d0+h2*exp(x(i))
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
         g(i)=4.0d0+h2*exp(x(i))
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
         g(i+1)=-1.0d0
         g(i+nd)=-1.0d0
      endif
      return
232   nd=int(sqrt(dble(n)))
      l=mod(i,nd)
      if(l.eq.0) then
         k=i/nd
         l=nd
      else
         k=int(i/nd)+1
      endif
      h=1.0d0/dble(nd+1)
      h2=h*h
      do 2020 j=1,n
        g(j)=0.0d0
2020  continue
      if(l.eq.1.and.k.eq.1) then
         g(1)=4.0d0+h2*x(1)*2.0d0
         g(2)=-1.0d0
         g(nd+1)=-1.0d0
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         g(l)=4.0d0+h2*x(l)*2.0d0
         g(l-1)=-1.0d0
         g(l+1)=-1.0d0
         g(l+nd)=-1.0d0
      endif
      if(l.eq.nd.and.k.eq.1) then
         g(nd)=4.0d0+h2*x(nd)*2.0d0
         g(nd-1)=-1.0d0
         g(nd+nd)=-1.0d0
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         g(i)=4.0d0+h2*x(i)*2.0d0
         g(i-nd)=-1.0d0
         g(i+1)=-1.0d0
         g(i+nd)=-1.0d0
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         g(i)=4.0d0+h2*x(i)*2.0d0
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
         g(i+nd)=-1.0d0
      endif
      if(l.eq.1.and.k.eq.nd) then
         g(i)=4.0d0+h2*x(i)*2.0d0
         g(i-nd)=-1.0d0
         g(i+1)=-1.0d0
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         g(i)=4.0d0+h2*x(i)*2.0d0
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
         g(i+1)=-1.0d0
      endif
      if(l.eq.nd.and.k.eq.nd) then
         g(i)=4.0d0+h2*x(i)*2.0d0
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
         g(i)=4.0d0+h2*x(i)*2.0d0
         g(i-nd)=-1.0d0
         g(i-1)=-1.0d0
         g(i+1)=-1.0d0
         g(i+nd)=-1.0d0
      endif
      return
      end subroutine eagu16

! subroutine tiud28                all systems                92/12/01
! portability : all systems
! 92/12/01 lu : original version
!
! purpose :
!  initial values of the variables for unconstrained minimization.
!  unconstrained and dense version.
!
! parameters :
!  ii  n  number of variables.
!  ii  na  number of approximated functions.
!  ro  x(n)  vector of variables.
!  ii  next  number of the test problem.
!  io  ierr  error indicator.
!
      subroutine tiud28(n,x,fmin,xmax,next,ierr)
      integer n,na,next,ierr
      double precision x(n),fmin,xmax
      double precision p,q,alf,bet,gam,f,s,s1,t,z(1000)
      integer i,j,k,m,n1
      double precision y(20),par
      common /empr28/ y,par,na,m
      double precision eta9
      parameter (eta9=1.0d60)
      fmin=0.0d0
      xmax=1.0d3
      ierr=0
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160, &
        170,180,190,200,210,220,230,250,310,320,330,350,370,390,400, &
        450,460,470,480,490,500,510,520,530,540,550,560,570,580,590, &
        600,610,620,630,720,740,750,760,780,790,810,830,840,860,870, &
        880,900,910,920,930,940,950,960,970,980,990,800,240,410,420, &
        650,660,670,680,690,340,360,380,430,440,270,280,290,300,710, &
        820),next
   10 if (n.lt.2) go to 999
      n=n-mod(n,2)
      do 11 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-1.2d0
        else
          x(i)=1.0d0
        endif
   11 continue
      return
   20 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 21 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-2.0d0
          if(i.le.4) x(i)=-3.0d0
        else
          x(i)=0.0d0
          if(i.le.4) x(i)=-1.0d0
        endif
   21 continue
      return
   30 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 31 i=1,n
        if(mod(i,4).eq.1) then
          x(i)=3.0d0
        elseif(mod(i,4).eq.2) then
          x(i)=-1.0d0
        elseif(mod(i,4).eq.3) then
          x(i)=0.0d0
        else
          x(i)=1.0d0
        endif
   31 continue
      return
   40 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 41 i=1,n
        x(i)=2.0d0
   41 continue
      x(1)=1.0d0
      return
   50 if (n.lt.3) go to 999
      do 51 i=1,n
        x(i)=-1.0d0
   51 continue
      return
   60 if (n.lt.7) go to 999
      do 61 i=1,n
        x(i)=-1.0d0
   61 continue
      return
   70 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 71 i=1,n
        x(i)=-1.0d0
   71 continue
      return
   80 if (n.lt.6) go to 999
      do 81 i=1,n
        x(i)=1.0d0/dble(n)
   81 continue
      return
   90 if (n.lt.6) go to 999
      do 91 i=1,n
        x(i)=1.0d0/dble(n)
   91 continue
      fmin=-eta9
      return
  100 if (n.lt.6) go to 999
      do 101 i=1,n
        x(i)=1.0d0
  101 continue
      fmin=-eta9
      return
  110 if (n.lt.5) go to 999
      n=n-mod(n,5)
      do 111 i=0,n-5,5
        x(i+1)=-1.0d0
        x(i+2)=-1.0d0
        x(i+3)=2.0d0
        x(i+4)=-1.0d0
        x(i+5)=-1.0d0
  111 continue
      x(1)=-2.0d0
      x(2)=2.0d0
      xmax=1.0d0
      return
  120 if (n.lt.2) go to 999
      n=n-mod(n,2)
      do 121 i=2,n,2
        x(i-1)=0.0d0
        x(i)=-1.0d0
  121 continue
      xmax=1.0d1
      return
  130 if (n.lt.2) go to 999
      n=n-mod(n,2)
      do 131 i=2,n,2
        x(i-1)=-1.0d0
        x(i)=1.0d0
  131 continue
      xmax=1.0d1
      return
  140 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 141 i=1,n
        q=p*dble(i)
        x(i)=q*(q-1.0d0)
  141 continue
      return
  150 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 151 i=1,n
        q=p*dble(i)
        x(i)=0.1d0*q*dble(n+1-i)
  151 continue
      fmin=-eta9
      return
  160 if (n.lt.3) go to 999
      do 161 i=1,n
        x(i)=1.0d0
  161 continue
      fmin=-eta9
      return
  170 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 171 i=1,n
        x(i)=dble(i)*dble(n+1-i)*p**2
  171 continue
      fmin=-eta9
      return
  180 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 181 i=1,n
        x(i)=dble(i)*dble(n+1-i)*p**2
  181 continue
      fmin=-eta9
      return
  190 if (n.lt.3) go to 999
      p=exp(2.0d0)/dble(n+1)
      do 191 i=1,n
        x(i)=(p*dble(i)+1.0d0)/3.0d0
  191 continue
      fmin=-eta9
      return
  200 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 201 i=1,n
        x(i)=p*dble(i)
  201 continue
      fmin=-eta9
      return
  210 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 211 i=1,n
        x(i)=p*dble(i)+1.0d0
  211 continue
      fmin=-eta9
      return
  220 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 221 i=1,n
        x(i)=dble(i)*dble(n+1-i)*p**2
  221 continue
      fmin=-eta9
      return
  230 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 231 i=1,n
        x(i)=dble(i)*p
  231 continue
      return
  250 if (n.lt.3) go to 999
      p=1.0d0/dble(n+1)
      do 251 i=1,n
        x(i)=dble(i)*p
  251 continue
      xmax=1.0d0
      return
  310 if(n.ge.2) then
      n=n-mod(n,2)
      na=n
      do 311 i=1,n
      if(mod(i,2).eq.1) then
      x(i)=-1.2d0
      else
      x(i)=1.0d0
      endif
  311 continue
      else
      ierr=1
      endif
      return
  320 if(n.ge.4) then
      n=n-mod(n,4)
      na=n
      do 321 i=1,n
      if (mod(i,4).eq.1) then
      x(i)=3.0d0
      else if(mod(i,4).eq.2) then
      x(i)=-1.0d0
      elseif(mod(i,4).eq.3) then
      x(i)=0.0d0
      else
      x(i)=1.0d0
      endif
  321 continue
      else
      ierr=1
      endif
      return
  330 if(n.ge.2) then
      na=n+1
      do 331 i=1,n
      x(i)=dble(i)
  331 continue
      else
      ierr=1
      endif
      return
  350 if(n.ge.2) then
      na=n+2
      do 351 i=1,n
      x(i)=1.0d0-dble(i)/dble(n)
  351 continue
      xmax=1.0d2
      else
      ierr=1
      endif
      return
  370 if(n.ge.2) then
      na=n
      do 371 i=1,n
      x(i)=0.5d0
  371 continue
      else
      ierr=1
      endif
      return
  390 if(n.ge.2) then
      na=n
      do 391 i=1,n
      x(i)=dble(i)/dble(n+1)
      x(i)=x(i)*(x(i)-1.0d0)
  391 continue
      else
      ierr=1
      endif
      return
  400 if(n.ge.2) then
      na=n
      do 401 i=1,n
      x(i)=-1.0d0
  401 continue
      else
      ierr=1
      endif
  450 if (n.lt.3) go to 999
      do 451 i=1,n
        x(i)=-1.0d0
  451 continue
      na=n
      return
  460 if (n.lt.6) go to 999
      do 461 i=1,n
        x(i)=-1.0d0
  461 continue
      na=n
      return
  470 if (n.lt.2) go to 999
      do 471 i=1,n-1
        x(i)=0.5d0
  471 continue
      x(n)=-2.0d0
      na=2*(n-1)
      return
  480 if (n.lt.4) go to 999
      n=n-mod(n,4)
      do 481 i=1,n
        x(i)=sin(dble(i))**2
  481 continue
      na=5*n
      return
  490 if (n.lt.4) go to 999
      n=n-mod(n,2)
      do 491 i=1,n
        x(i)=5.0d0
  491 continue
      na=3*(n-2)
      return
  500 if (n.lt.2) go to 999
      do 501 i=1,n
        x(i)=0.2d0
  501 continue
      na=2*(n-1)
      return
  510 continue
      if (n.lt.2) go to 999
      n=n-mod(n,2)
      do 511 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-0.8d0
        else
          x(i)=-0.8d0
        endif
  511 continue
      na=2*(n-1)
      return
  520 continue
      if (n.lt.5) go to 999
      if (mod(n-5,3).ne.0) n=n-mod(n-5,3)
      do 521 i=1,n
        x(i)=-1.0d0
  521 continue
      na=6*((n-5)/3+1)
      return
  530 continue
      if (n.lt.5) go to 999
      if (mod(n-5,3).ne.0) n=n-mod(n-5,3)
      do 531 i=1,n
        x(i)=-1.0d0
  531 continue
      na=7*((n-5)/3+1)
      return
  540 continue
      if (n.lt.4) go to 999
      do 541 i=1,n
        if (mod(i,4).eq.1) then
          x(i)=-0.8d0
        else if (mod(i,4).eq.2) then
          x(i)= 1.2d0
        else if (mod(i,4).eq.3) then
          x(i)=-1.2d0
        else
          x(i)= 0.8d0
        endif
  541 continue
      y(1)=14.4d0
      y(2)=6.8d0
      y(3)=4.2d0
      y(4)=3.2d0
  542 if (mod(n-4,2).ne.0) n=n-mod(n-4,2)
      na=4*((n-4)/2+1)
      return
  550 continue
      if (n.lt.4) go to 999
      do 551 i=1,n
        if (mod(i,4).eq.1) then
          x(i)=-0.8d0
        else if (mod(i,4).eq.2) then
          x(i)= 1.2d0
        else if (mod(i,4).eq.3) then
          x(i)=-1.2d0
        else
          x(i)= 0.8d0
        endif
  551 continue
      y(1)=35.8d0
      y(2)=11.2d0
      y(3)=6.2d0
      y(4)=4.4d0
      go to 542
  560 continue
      if (n.lt.4) go to 999
      do 561 i=1,n
        if (mod(i,4).eq.1) then
          x(i)=-0.8d0
        else if (mod(i,4).eq.2) then
          x(i)= 1.2d0
        else if (mod(i,4).eq.3) then
          x(i)=-1.2d0
        else
          x(i)= 0.8d0
        endif
  561 continue
      y(1)=30.6d0
      y(2)=72.2d0
      y(3)=124.4d0
      y(4)=187.4d0
      go to 542
  570 if (n.lt.4) go to 999
      n=n-mod(n,2)
      na=n
      do 571 i=1,n
        if (mod(i,8).eq.1) x(i)=1.0d-1
        if (mod(i,8).eq.2.or.mod(i,8).eq.0) x(i)=2.0d-1
        if (mod(i,8).eq.3.or.mod(i,8).eq.7) x(i)=3.0d-1
        if (mod(i,8).eq.4.or.mod(i,8).eq.6) x(i)=4.0d-1
        if (mod(i,8).eq.5) x(i)=5.0d-1
  571 continue
      return
  580 if (n.lt.3) go to 999
      na=n
      do 581 i=1,n
        x(i)=1.2d1
  581 continue
      xmax=1.0d1
      return
  590 if (n.lt.7) go to 999
      na=n
      do 591 i=1,n
        x(i)=-1.0d0
  591 continue
      return
  600 if (n.lt.3) go to 999
      na=n
      do 601 i=1,n
        x(i)=dble(i)/dble(n+1)
        x(i)=x(i)*(x(i)-1.0d0)
  601 continue
      return
  610 continue
      if (n.lt.5) go to 999
      if (mod(n-5,3).ne.0) n=n-mod(n-5,3)
      do 611 i=1,n
        x(i)=-1.0d0
  611 continue
      na=7*((n-5)/3+1)
      return
  620 if (n.lt.3) go to 999
      do 621 i=1,n
        if(mod(i,2).eq.1) then
          x(i)=-1.2d0
        else
          x(i)= 1.0d0
        endif
  621 continue
      na=2*(n-1)
      return
  630 if (n.lt.14) go to 999
      n=n/2
      do 631 i=1,n
        x(i)=5.0d0
  631 continue
      y(1)=sin(1.0d0)
      na=13*(n-6)
      return
  720 if(n.ge.5) then
      n=n-mod(n,2)
      na=n
      do 721 i=1,n
      if (mod(i,8).eq.1) x(i)=1.0d-1
      if (mod(i,8).eq.2.or.mod(i,8).eq.0) x(i)=2.0d-1
      if (mod(i,8).eq.3.or.mod(i,8).eq.7) x(i)=3.0d-1
      if (mod(i,8).eq.4.or.mod(i,8).eq.6) x(i)=4.0d-1
      if (mod(i,8).eq.5) x(i)=5.0d-1
  721 continue
      else
      ierr=1
      endif
      return
  740 if(n.ge.3) then
      do 741 i=1,n
      x(i)=0.0d0
  741 continue
      na=n
      else
      ierr=1
      endif
      return
  750 if(n.ge.3) then
      if (mod(n,2).ne.1) n=n-1
      na=n
      do 751 i=1,n
      x(i)=1.0d0
  751 continue
      else
      ierr=1
      endif
      return
  760 if(n.ge.3) then
      do 761 i=1,n
      x(i)=-1.0d0
  761 continue
      na=n
      else
      ierr=1
      endif
      return
  780 if(n.ge.5) then
      na=n
      do 781 i=1,n
      x(i)=-2.0d0
  781 continue
      else
      ierr=1
      endif
      return
  790 if(n.ge.7) then
      na=n
      do 791 i=1,n
      x(i)=-3.0d0
  791 continue
      xmax=1.0d1
      else
      ierr=1
      endif
      return
  810 if(n.ge.2) then
      n=n-mod(n,2)
      na=n
      do 811 i=1,n
      if(mod(i,2).eq.1) then
      x(i)=9.0d1
      else
      x(i)=6.0d1
      endif
  811 continue
      else
      ierr=1
      endif
      return
  830 if(n.ge.4) then
      n=n-mod(n,4)
      na=n
      do 831 i=1,n
      if (mod(i,4).eq.1) then
      x(i)=1.0d0
      else if(mod(i,4).eq.2) then
      x(i)=2.0d0
      elseif(mod(i,4).eq.3) then
      x(i)=2.0d0
      else
      x(i)=2.0d0
      endif
  831 continue
      xmax=1.0d1
      else
      ierr=1
      endif
      return
  840 if(n.ge.3) then
      do 841 i=1,n
      x(i)=-1.0d0
  841 continue
      na=n
      else
      ierr=1
      endif
      return
  860 if(n.ge.2) then
      n=n-mod(n,2)
      na=n
      do 861 i=1,n
      if (mod(i,2).eq.1) then
      x(i)=0.0d0
      else
      x(i)=1.0d0
      endif
  861 continue
      else
      ierr=1
      endif
      return
  870 if (n.ge.4) then
      n=n-mod(n,4)
      na=n
      do 871 i=2,n,2
      x(i-1)=-3.0d0
      x(i)=-1.0d0
  871 continue
      else
      ierr=1
      endif
      return
  880 if (n.ge.3) then
      do 881 i=1,n
      x(i)=1.5d0
  881 continue
      xmax=1.0d0
      na=n
      else
        ierr=1
      endif
      return
  900 if(n.ge.3) then
      do 901 i=1,n
      x(i)=1.0d1
  901 continue
      na=n
      else
      ierr=1
      endif
      return
  910 if (n.ge.3) then
      do 911 i=1,n
      x(i)=1.0d0
  911 continue
      par=1.0d1
      na=n
      else
        ierr=1
      endif
      return
  920 if (n.ge.5) then
      par=5.0d2/dble(n+2)
      na=n
      do 921 i=1,n
      x(i)=((dble(i)+0.5d0)/dble(n+2)-0.5d0)**2
  921 continue
      else
        ierr=1
      endif
      return
  930 if (n.ge.10) then
      n=n-mod(n,2)
      m=n/2
      par=5.0d2
      na=n
      do 931 i=1,m
      x(i)=(dble(i)/dble(m+1)-0.5d0)**2
  931 continue
      do 932 i=m+1,n
      k=i-m
      x(i)=dble(k)/dble(m+1)-0.5d0
  932 continue
      else
        ierr=1
      endif
      return
  940 if (n.ge.16) then
      m=int(sqrt(dble(n)))
      par=6.8d0/dble(m+1)**2
      n=m*m
      k=0
      do 942 j=1,m
      do 941 i=1,m
      k=k+1
      x(k)=0.0d0
  941 continue
  942 continue
      na=n
      else
      ierr=1
      endif
      return
  950 if (n.ge.16) then
      m=int(sqrt(dble(n)))
      par=1.0d0/dble(m+1)**2
      n=m*m
      k=0
      do 952 j=1,m
      do 951 i=1,m
      k=k+1
      x(k)=-1.0d0
  951 continue
  952 continue
      na=n
      else
      ierr=1
      endif
      return
  960 if (n.ge.16) then
      m=int(sqrt(dble(n)))
      par=1.0d0/dble(m+1)**2
      n=m*m
      k=0
      do 962 j=1,m
      do 961 i=1,m
      k=k+1
      x(k)=0.0d0
  961 continue
  962 continue
      na=n
      else
      ierr=1
      endif
      return
  970 if (n.ge.16) then
      m=int(sqrt(dble(n)))
      par=5.0d1/dble(m+1)
      n=m*m
      k=0
      do 972 j=1,m
      do 971 i=1,m
      k=k+1
      x(k)=1.0d0-dble(i)*dble(j)/dble(m+1)**2
  971 continue
  972 continue
      na=n
      else
      ierr=1
      endif
      return
  980 if (n.ge.16) then
      m=int(sqrt(dble(n)))
      par=1.0d0/dble(m+1)
      n=m*m
      k=0
      do 982 j=1,m
      do 981 i=1,m
      k=k+1
      x(k)=0.0d0
  981 continue
  982 continue
      na=n
      else
      ierr=1
      endif
      return
  990 if (n.ge.16) then
      m=int(sqrt(dble(n)))
      n=m*m
      par=500.0d0/dble(m+2)**4
      k=0
      do 992 j=1,m
      do 991 i=1,m
      k=k+1
      x(k)=0.0d0
  991 continue
  992 continue
      na=n
      else
        ierr=1
      endif
      return
  800 if (n.ge.16) then
      m=int(sqrt(dble(n)))
      n=m*m
      par=500.0d0
      k=0
      do 802 j=1,m
      do 801 i=1,m
      k=k+1
      x(k)=0.0d0
  801 continue
  802 continue
      na=n
      else
        ierr=1
      endif
      return
  240 if(n.ge.2) then
      na=n
      do 241 i=1,n
      if (mod(i,2).eq.1) then
      x(i)=1.0d0
      else
      x(i)=3.0d0
      endif
  241 continue
      else
      ierr=1
      endif
      return
  410 n1=n-1
      do 411 i=1,n1
        x(i)=-1.2d0
  411 continue
      x(n)=-1.0d0
      return
  420 do 421 i=1,n
        x(i)=2.0d0
  421 continue
      return
  650 do 651 i=1,n
        x(i)=1.5d0
  651 continue
      return
  660 do 661 i=1,n
        x(i)=0.0d0
  661 continue
      return
  670 do 671 i=1,n
        x(i)=-1.0d0
  671 continue
      return
  680 do 681 i=1,n
        x(i)=-1.0d0
  681 continue
      return
  690 do 691 i=1,n
        x(i)=0.5d0
  691 continue
      return
  340 do 341 i=1,n
        x(i)=0.5d0
  341 continue
      return
  360 do 361 i=1,n
        x(i)=1.0d0
  361 continue
      return
  380 do 381 i=1,n
        x(i)=-1.0d0
  381 continue
      return
  430 alf=5
      bet=14
      gam=3
      do 431 i=1,n
        x(i)=0.0d0
  431 continue
      do 433 i=1,n
      f=dble(bet*n)*x(i)+(dble(i)-dble(n)/2.0d0)**gam
      do 432 j=1,n
        if(j.ne.i) then
        t=sqrt(x(j)**2+dble(i)/dble(j))
        s1=log(t)
        f=f+t*(sin(s1)**alf+cos(s1)**alf)
      endif
  432 continue
        z(i)=-f
  433 continue
      s=dble(bet*n)/dble(bet**2*n**2-(alf+1)**2*(n-1)**2)
      do 434 i=1,n
        x(i)=s*z(i)
  434 continue
      return
  440 do 441 i=1,n
        x(i)=1.0d0
  441 continue
      return
  270 do 271 i=1,n
        x(i)=1.0d0
  271 continue
      return
  280 do 281 i=1,n
        x(i)=1.0d0
  281 continue
      return
  290 do 291 i=1,n
        x(i)=1.0d0
  291 continue
      return
  300 t=dble(2)/dble(n+2)
      n1=n/2
      do 301 i=1,n1
        s=dble(i)*t
        x(i)=s*(1.0d0-s)
        x(n1+i)=x(i)
  301 continue
      return
  710 n1=int(sqrt(dble(n)))
      n=n1*n1
      do 711 i=1,n
        x(i)=1.0d0
  711 continue
      return
  820 n1=int(sqrt(dble(n)))
      n=n1*n1
      do 821 i=1,n
        x(i)=1.0d0
  821 continue
      return
  999 ierr=1
      return
      end subroutine tiud28

! subroutine tfbu28                all systems                92/12/01
! portability : all systems
! 92/12/01 lu : original version
!
! purpose :
!  values of model functions for unconstrained minimization.
!  universal version.
!
! parameters :
!  ii  n  number of variables.
!  ii  na  number of approximated functions.
!  ri  x(n)  vector of variables.
!  ro  f  value of the model function.
!  ii  next  number of the test problem.
!
      subroutine tfbu28(n,x,f,g,next)
      use matrix_routines
      integer n,na,next
      double precision x(n),g(n),f
      double precision fa
      double precision a,b,c,d,e,p,q,r,u,v,w,alfa,sx(1000)
      double precision a1,a2,a3,a4,ex,d1s,d2s,h,pi
      data pi /3.14159265358979323d0/
      double precision ga1(2),ga2(2),ga3(6),ga4(6)
      double precision al,al1,al2,alf,be,be1,be2,bet,ca,cb,ff,fg,ga, &
       gam,h2,s,s1,s2,s3,t,t1
      integer i,j,k,l,m,ia,ib,ic,i1,i2,j1,j2,ka,la,n1,nd
      double precision y(20),par
      common /empr28/ y,par,na,m
      f=0.0d0
      call mxvset(n,0.0d0,g)
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160, &
        170,180,190,200,210,220,230,250,310,320,330,350,370,390,400, &
        450,460,470,480,490,500,510,520,530,540,550,560,570,580,590, &
        600,610,620,630,720,740,750,760,780,790,810,830,840,860,870, &
        880,900,910,920,930,940,950,960,970,980,990,800,240,410,420, &
        650,660,670,680,690,340,360,380,430,440,270,280,290,300,710, &
        820),next
   10 do 11 j=2,n
      a=x(j-1)**2-x(j)
      b=x(j-1)-1.0d0
      f=f+1.0d2*a**2+b**2
      g(j-1)=g(j-1)+4.0d2*x(j-1)*a+2.0d0*b
      g(j)=g(j)-2.0d2*a
   11 continue
      return
   20 do 21 j=2,n-2,2
      a=x(j-1)**2-x(j)
      b=x(j-1)-1.0d0
      c=x(j+1)**2-x(j+2)
      d=x(j+1)-1.0d0
      u=x(j)+x(j+2)-2.0d0
      v=x(j)-x(j+2)
      f=f+1.0d2*a**2+b**2+9.0d1*c**2+d**2+ &
        1.0d1*u**2+0.1d0*v**2
      g(j-1)=g(j-1)+4.0d2*x(j-1)*a+2.0d0*b
      g(j)=g(j)-2.0d2*a+2.0d1*u+0.2d0*v
      g(j+1)=g(j+1)+3.6d2*x(j+1)*c+2.0d0*d
      g(j+2)=g(j+2)-1.8d2*c+2.0d1*u-0.2d0*v
   21 continue
      return
   30 do 31 j=2,n-2,2
      a=x(j-1)+1.0d1*x(j)
      b=x(j+1)-x(j+2)
      c=x(j)-2.0d0*x(j+1)
      d=x(j-1)-x(j+2)
      f=f+a**2+5.0d0*b**2+c**4+1.0d1*d**4
      g(j-1)=g(j-1)+2.0d0*a+4.0d1*d**3
      g(j)=g(j)+2.0d1*a+4.0d0*c**3
      g(j+1)=g(j+1)-8.0d0*c**3+1.0d1*b
      g(j+2)=g(j+2)-4.0d1*d**3-1.0d1*b
   31 continue
      return
   40 do 41 j=2,n-2,2
      a=exp(x(j-1))
      b=a-x(j)
      d=x(j)-x(j+1)
      p=x(j+1)-x(j+2)
      c=cos(p)
      q=sin(p)/cos(p)
      u=x(j-1)
      v=x(j+2)-1.0d0
      f=f+b**4+1.0d2*d**6+q**4+u**8+v**2
      b=4.0d0*b**3
      d=6.0d2*d**5
      q=4.0d0*q**3/c**2
      g(j-1)=g(j-1)+a*b+8.0d0*u**7
      g(j)=g(j)+d-b
      g(j+1)=g(j+1)+q-d
      g(j+2)=g(j+2)+2.0d0*v-q
   41 continue
      return
   50 p=7.0d0/3.0d0
      do 51 j=1,n
      a=(3.0d0-2.0d0*x(j))*x(j)+1.0d0
      if (j.gt.1) a=a-x(j-1)
      if (j.lt.n) a=a-x(j+1)
      f=f+abs(a)**p
      b=p*abs(a)**(p-1.0d0)*sign(1.0d0,a)
      g(j)=g(j)+b*(3.0d0-4.0d0*x(j))
      if (j.gt.1) g(j-1)=g(j-1)-b
      if (j.lt.n) g(j+1)=g(j+1)-b
   51 continue
      return
   60 p=7.0d0/3.0d0
      do 63 j=1,n
      a=(2.0d0+5.0d0*x(j)**2)*x(j)+1.0d0
      do 61 i=max(1,j-5),min(n,j+1)
      a=a+x(i)*(1.0d0+x(i))
   61 continue
      b=p*abs(a)**(p-1.0d0)*sign(1.0d0,a)
      f=f+abs(a)**p
      g(j)=g(j)+b*(2.0d0+1.5d1*x(j)**2)
      do 62 i=max(1,j-5),min(n,j+1)
      g(i)=g(i)+b*(1.0d0+2.0d0*x(i))
   62 continue
   63 continue
      return
   70 p=7.0d0/3.0d0
      k=n/2
      do 71 j=1,n
      a=(3.0d0-2.0d0*x(j))*x(j)+1.0d0
      if (j.gt.1) a=a-x(j-1)
      if (j.lt.n) a=a-x(j+1)
      b=p*abs(a)**(p-1.0d0)*sign(1.0d0,a)
      f=f+abs(a)**p
      g(j)=g(j)+b*(3.0d0-4.0d0*x(j))
      if (j.gt.1) g(j-1)=g(j-1)-b
      if (j.lt.n) g(j+1)=g(j+1)-b
      if (j.le.k) then
      f=f+abs(x(j)+x(j+k))**p
      a=x(j)+x(j+k)
      b=p*abs(a)**(p-1.0d0)*sign(1.0d0,a)
      g(j)=g(j)+b
      g(j+k)=g(j+k)+b
      endif
   71 continue
      return
   80 k=n/2
      do 83 j=1,n
      p=0.0d0
      do 81 i=j-2,j+2
      if (i.lt.1.or.i.gt.n) go to 81
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      p=p+a*sin(x(i))+b*cos(x(i))
   81 continue
      if (j.gt.k) then
      i=j-k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      p=p+a*sin(x(i))+b*cos(x(i))
      else
      i=j+k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      p=p+a*sin(x(i))+b*cos(x(i))
      endif
      f=f+(dble(n+j)-p)**2/dble(n)
      p=2.0d0*(dble(n+j)-p)/dble(n)
      do 82 i=j-2,j+2
      if (i.lt.1.or.i.gt.n) go to 82
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      g(i)=g(i)-p*(a*cos(x(i))-b*sin(x(i)))
   82 continue
      if (j.gt.k) then
      i=j-k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      g(i)=g(i)-p*(a*cos(x(i))-b*sin(x(i)))
      else
      i=j+k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      g(i)=g(i)-p*(a*cos(x(i))-b*sin(x(i)))
      endif
   83 continue
      return
   90 k=n/2
      q=1.0d0/dble(n)
      do 92 j=1,n
      p=0.0d0
      do 91 i=j-2,j+2
      if (i.lt.1.or.i.gt.n) go to 91
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      p=p+a*sin(x(i))+b*cos(x(i))
      g(i)=g(i)+q*(a*cos(x(i))-b*sin(x(i)))
   91 continue
      if (j.gt.k) then
      i=j-k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      p=p+a*sin(x(i))+b*cos(x(i))
      g(i)=g(i)+q*(a*cos(x(i))-b*sin(x(i)))
      else
      i=j+k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=dble(i+j)/1.0d1
      p=p+a*sin(x(i))+b*cos(x(i))
      g(i)=g(i)+q*(a*cos(x(i))-b*sin(x(i)))
      endif
      f=f+(p+dble(j)*(1.0d0-cos(x(j))))*q
      g(j)=g(j)+q*dble(j)*sin(x(j))
   92 continue
      return
  100 k=n/2
      do 102 j=1,n
      p=0.0d0
      q=1.0d0+dble(j)/1.0d1
      do 101 i=j-2,j+2
      if (i.lt.1.or.i.gt.n) go to 101
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=1.0d0+dble(i)/1.0d1
      c=dble(i+j)/1.0d1
      p=p+a*sin(q*x(j)+b*x(i)+c)
      r=a*cos(q*x(j)+b*x(i)+c)/dble(n)
      g(j)=g(j)+r*q
      g(i)=g(i)+r*b
  101 continue
      if (j.gt.k) then
      i=j-k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=1.0d0+dble(i)/1.0d1
      c=dble(i+j)/1.0d1
      p=p+a*sin(q*x(j)+b*x(i)+c)
      r=a*cos(q*x(j)+b*x(i)+c)/dble(n)
      g(j)=g(j)+r*q
      g(i)=g(i)+r*b
      else
      i=j+k
      a=5.0d0*(1.0d0+mod(i,5)+mod(j,5))
      b=1.0d0+dble(i)/1.0d1
      c=dble(i+j)/1.0d1
      p=p+a*sin(q*x(j)+b*x(i)+c)
      r=a*cos(q*x(j)+b*x(i)+c)/dble(n)
      g(j)=g(j)+r*q
      g(i)=g(i)+r*b
      endif
      f=f+p
  102 continue
      f=f/dble(n)
      return
  110 p=-0.2008d-2
      q=-0.1900d-2
      r=-0.0261d-2
      do 112 i=0,n-5,5
      a=1.0d0
      b=0.0d0
      do 111 j=1,5
      a=a*x(i+j)
      b=b+x(i+j)**2
  111 continue
      w=exp(a)
      a=a*w
      b=b-1.0d1-p
      c=x(i+2)*x(i+3)-5.0d0*x(i+4)*x(i+5)-q
      d=x(i+1)**3+x(i+2)**3+1.0d0-r
      f=f+w+1.0d1*(b**2+c**2+d**2)
      g(i+1)=g(i+1)+a/x(i+1)+2.0d1*(2.0d0*b*x(i+1)+ &
        3.0d0*d*x(i+1)**2)
      g(i+2)=g(i+2)+a/x(i+2)+2.0d1*(2.0d0*b*x(i+2)+ &
        c*x(i+3)+3.0d0*d*x(i+2)**2)
      g(i+3)=g(i+3)+a/x(i+3)+2.0d1*(2.0d0*b*x(i+3)+c*x(i+2))
      g(i+4)=g(i+4)+a/x(i+4)+2.0d1*(2.0d0*b*x(i+4)- &
        5.0d0*c*x(i+5))
      g(i+5)=g(i+5)+a/x(i+5)+2.0d1*(2.0d0*b*x(i+5)-5.0d0*c*x(i+4))
  112 continue
      return
  120 c=0.0d0
      do 121 j=2,n,2
      a=x(j-1)-3.0d0
      b=x(j-1)-x(j)
      c=c+a
      f=f+1.0d-4*a**2-b+exp(2.0d1*b)
      g(j-1)=g(j-1)+2.0d-4*a-1.0d0+2.0d1*exp(2.0d1*b)
      g(j)=g(j)+1.0d0-2.0d1*exp(2.0d1*b)
  121 continue
      f=f+c**2
      do 122 j=2,n,2
      g(j-1)=g(j-1)+2.0d0*c
  122 continue
      return
  130 do 131 j=2,n,2
      a=x(j)**2
      if (a.eq.0.0d0) a=1.0d-60
      b=x(j-1)**2
      if (b.eq.0.0d0) b=1.0d-60
      c=a+1.0d0
      d=b+1.0d0
      f=f+b**c+a**d
      p=0.0d0
      if (a.gt.p) p=log(a)
      q=0.0d0
      if (b.gt.q) q=log(b)
      g(j-1)=g(j-1)+2.0d0*x(j-1)*(c*b**a+p*a**d)
      g(j)=g(j)+2.0d0*x(j)*(d*a**b+q*b**c)
  131 continue
      return
  140 p=1.0d0/dble(n+1)
      q=0.5d0*p**2
      do 141 j=1,n
      a=2.0d0*x(j)+q*(x(j)+dble(j)*p+1.0d0)**3
      if(j.gt.1) a=a-x(j-1)
      if(j.lt.n) a=a-x(j+1)
      f=f+a**2
      g(j)=g(j)+a*(4.0d0+6.0d0*q*(x(j)+dble(j)*p+1.0d0)**2.0d0)
      if(j.gt.1) g(j-1)=g(j-1)-2.0d0*a
      if(j.lt.n) g(j+1)=g(j+1)-2.0d0*a
  141 continue
      return
  150 p=1.0d0/dble(n+1)
      q=2.0d0/p
      r=2.0d0*p
      do 151 j=2,n
      a=x(j-1)-x(j)
      f=f+q*x(j-1)*a
      g(j-1)=g(j-1)+q*(2.0d0*x(j-1)-x(j))
      g(j)=g(j)-q*x(j-1)
      if (abs(a).le.1.0d-6) then
      f=f+r*exp(x(j))*(1.0d0+a/2.0d0*(1.0d0+a/3.0d0*(1.0d0+ &
        a/4.0d0)))
      g(j-1)=g(j-1)+r*exp(x(j))*(1.0d0/2.0d0+a*(1.0d0/3.0d0+ &
        a/8.0d0))
      g(j)=g(j)+r*exp(x(j))*(1.0d0/2.0d0+a*(1.0d0/6.0d0+ &
        a/24.0d0))
      else
      b=exp(x(j-1))-exp(x(j))
      f=f+r*b/a
      g(j-1)=g(j-1)+r*(exp(x(j-1))*a-b)/a**2
      g(j)=g(j)-r*(exp(x(j))*a-b)/a**2
      endif
  151 continue
      f=f+q*x(n)**2+r*(exp(x(1))-1.0d0)/x(1) &
                   +r*(exp(x(n))-1.0d0)/x(n)
      g(1)=g(1)+r*(exp(x(1))*(x(1)-1.0d0)+1.0d0)/x(1)**2
      g(n)=g(n)+2.0d0*q*x(n) &
               +r*(exp(x(n))*(x(n)-1.0d0)+1.0d0)/x(n)**2
      return
  160 do 161 j=1,n
      a=dble(j)*(1.0d0-cos(x(j)))
      if(j.gt.1) a=a+dble(j)*sin(x(j-1))
      if(j.lt.n) a=a-dble(j)*sin(x(j+1))
      f=f+a
      a=dble(j)*sin(x(j))
      g(j)=g(j)+a
      if(j.gt.1) g(j-1)=g(j-1)+dble(j)*cos(x(j-1))
      if(j.lt.n) g(j+1)=g(j+1)-dble(j)*cos(x(j+1))
  161 continue
      return
  170 p=1.0d0/dble(n+1)
      do 171 j=1,n
      if (j.eq.1) then
      f=f+0.25d0*x(j)**2/p+1.25d-1*x(j+1)**2/p+ &
       p*(exp(x(j))-1.0d0)
      g(j)=g(j)+0.5d0*x(j)/p+p*exp(x(j))
      g(j+1)=g(j+1)+0.25d0*x(j+1)/p
      else if (j.eq.n) then
      f=f+0.25d0*x(j)**2/p+1.25d-1*x(j-1)**2/p+ &
       p*(exp(x(j))-1.0d0)
      g(j)=g(j)+0.5d0*x(j)/p+p*exp(x(j))
      g(j-1)=g(j-1)+0.25d0*x(j-1)/p
      else
      f=f+1.25d-1*(x(j+1)-x(j-1))**2/p+p*(exp(x(j))-1.0d0)
      a=0.25d0*(x(j+1)-x(j-1))/p
      g(j)=g(j)+p*exp(x(j))
      g(j-1)=g(j-1)-a
      g(j+1)=g(j+1)+a
      endif
  171 continue
      return
  180 p=1.0d0/dble(n+1)
      do 181 j=1,n
      q=dble(j)*p
      if (j.eq.1) then
      f=f+0.5d0*x(j)**2/p+0.25d0*x(j+1)**2/p- &
       p*(x(j)**2+2.0d0*x(j)*q)
      g(j)=g(j)+x(j)/p-2.0d0*p*(x(j)+q)
      g(j+1)=g(j+1)+0.5d0*x(j+1)/p
      else if (j.eq.n) then
      f=f+0.5d0*x(j)**2/p+0.25d0*x(j-1)**2/p- &
       p*(x(j)**2+2.0d0*x(j)*q)
      g(j)=g(j)+x(j)/p-2.0d0*p*(x(j)+q)
      g(j-1)=g(j-1)+0.5d0*x(j-1)/p
      else
      f=f+2.5d-1*(x(j+1)-x(j-1))**2/p-p*(x(j)**2+ &
       2.0d0*x(j)*q)
      a=0.5d0*(x(j+1)-x(j-1))/p
      g(j)=g(j)-2.0d0*p*(x(j)+q)
      g(j-1)=g(j-1)-a
      g(j+1)=g(j+1)+a
      endif
  181 continue
      return
  190 p=1.0d0/dble(n+1)
      do 191 j=1,n
      q=exp(2.0d0*dble(j)*p)
      if (j.eq.1) then
      r=1.0d0/3.0d0
      f=f+0.5d0*(x(j)-r)**2/p+7.0d0*r**2+ &
       2.5d-1*(x(j+1)-r)**2/p+p*(x(j)**2+2.0d0*x(j)*q)
      a=0.5d0*(x(j+1)-r)/p
      g(j)=g(j)+2.0d0*p*(x(j)+q)+(x(j)-r)/p
      g(j+1)=g(j+1)+a
      else if (j.eq.n) then
      r=exp(2.0d0)/3.0d0
      f=f+0.5d0*(x(j)-r)**2/p+7.0d0*r**2+ &
       2.5d-1*(x(j-1)-r)**2/p+p*(x(j)**2+2.0d0*x(j)*q)
      a=0.5d0*(x(j-1)-r)/p
      g(j)=g(j)+2.0d0*p*(x(j)+q)+(x(j)-r)/p
      g(j-1)=g(j-1)+a
      else
      f=f+2.5d-1*(x(j+1)-x(j-1))**2/p+p*(x(j)**2+ &
       2.0d0*x(j)*q)
      a=0.5d0*(x(j+1)-x(j-1))/p
      g(j)=g(j)+2.0d0*p*(x(j)+q)
      g(j-1)=g(j-1)-a
      g(j+1)=g(j+1)+a
      endif
  191 continue
      return
  200 p=1.0d0/dble(n+1)
      do 201 j=1,n
      a=exp(-2.0d0*x(j)**2)
      if (j.eq.1) then
      f=f+(0.5d0*x(j)**2/p-p)+ &
       (2.5d-1*x(j+1)**2/p-p)*a
      b=0.5d0*x(j+1)/p
      g(j)=g(j)+x(j)/p-4.0d0*x(j)*a*p*(b**2-1.0d0)
      g(j+1)=g(j+1)+a*b
      else if (j.eq.n) then
      f=f+(0.5d0*x(j)**2/p-p)*exp(-2.0d0)+ &
       (2.5d-1*x(j-1)**2/p-p)*a
      b=0.5d0*x(j-1)/p
      g(j)=g(j)+x(j)/p*exp(-2.0d0)-4.0d0*x(j)*a*p*(b**2-1.0d0)
      g(j-1)=g(j-1)+a*b
      else
      f=f+(2.5d-1*(x(j+1)-x(j-1))**2/p-p)*a
      b=0.5d0*(x(j+1)-x(j-1))/p
      g(j)=g(j)-4.0d0*x(j)*a*p*(b**2-1.0d0)
      g(j-1)=g(j-1)-a*b
      g(j+1)=g(j+1)+a*b
      endif
  201 continue
      return
  210 p=1.0d0/dble(n+1)
      do 211 j=1,n
      if (j.eq.1) then
      a=0.5d0*(x(j+1)-1.0d0)/p
      b=(x(j)-1.0d0)/p
      u=atan(a)
      v=atan(b)
      f=f+p*(x(j)**2+a*u-log(sqrt(1.0d0+a**2)))+ &
         0.5d0*p*(1.0d0+b*v-log(sqrt(1.0d0+b**2)))
      g(j)=g(j)+2.0d0*p*x(j)+0.5d0*v
      g(j+1)=g(j+1)+0.5d0*u
      else if (j.eq.n) then
      a=0.5d0*(2.0d0-x(j-1))/p
      b=(2.0d0-x(j))/p
      u=atan(a)
      v=atan(b)
      f=f+p*(x(j)**2+a*u-log(sqrt(1.0d0+a**2)))+ &
         0.5d0*p*(4.0d0+b*v-log(sqrt(1.0d0+b**2)))
      g(j)=g(j)+2.0d0*p*x(j)-0.5d0*v
      g(j-1)=g(j-1)-0.5d0*u
      else
      a=0.5d0*(x(j+1)-x(j-1))/p
      u=atan(a)
      f=f+p*(x(j)**2+a*u-log(sqrt(1.0d0+a**2)))
      g(j)=g(j)+2.0d0*p*x(j)
      g(j-1)=g(j-1)-0.5d0*u
      g(j+1)=g(j+1)+0.5d0*u
      endif
  211 continue
      return
  220 p=1.0d0/dble(n+1)
      do 221 j=1,n
      if (j.eq.1) then
      a= 0.5d0*x(j+1)/p
      b= x(j)/p
      f=f+p*(1.0d2*(x(j)-a**2)**2+(1.0d0-a)**2)+ &
         0.5d0*p*(1.0d2*b**4+(1.0d0-b)**2)
      g(j)=g(j)+2.0d2*p*(x(j)-a**2)+2.0d2*b**3-(1.0d0-b)
      g(j+1)=g(j+1)-2.0d2*(x(j)-a**2)*a-(1.0d0-a)
      else if (j.eq.n) then
      a=-0.5d0*x(j-1)/p
      b=-x(j)/p
      f=f+p*(1.0d2*(x(j)-a**2)**2+(1.0d0-a)**2)+ &
         0.5d0*p*(1.0d2*b**4+(1.0d0-b)**2)
      g(j)=g(j)+2.0d2*p*(x(j)-a**2)-2.0d2*b**3+(1.0d0-b)
      g(j-1)=g(j-1)+2.0d2*(x(j)-a**2)*a+(1.0d0-a)
      else
      a=0.5d0*(x(j+1)-x(j-1))/p
      f=f+p*(1.0d2*(x(j)-a**2)**2+(1.0d0-a)**2)
      g(j)=g(j)+2.0d2*p*(x(j)-a**2)
      g(j-1)=g(j-1)+2.0d2*(x(j)-a**2)*a+(1.0d0-a)
      g(j+1)=g(j+1)-2.0d2*(x(j)-a**2)*a-(1.0d0-a)
      endif
  221 continue
      return
  230 a=1.0d0
      b=1.0d-3
      c=0.0d0
      d=0.0d0
      do 231 j=1,n
      c=c+(x(j)-1.0d0)**2
      d=d+x(j)**2
  231 continue
      f=a*c+b*(d-2.5d-1)**2
      do 232 j=1,n
      g(j)=2.0d0*a*(x(j)-1.0d0)+4.0d0*b*(d-2.5d-1)*x(j)
  232 continue
      return
  250 a=1.0d0
      b=0.0d0
      c=0.0d0
      d=0.0d0
      f=0.0d0
      u=exp(x(n))
      v=exp(x(n-1))
      do 251 j=1,n
      if (j.le.n/2) f=f+(x(j)-1.0d0)**2
      if (j.le.n-2) then
      b=b+(x(j)+2.0d0*x(j+1)+1.0d1*x(j+2)-1.0d0)**2
      c=c+(2.0d0*x(j)+x(j+1)-3.0d0)**2
      endif
      d=d+x(j)**2-dble(n)
  251 continue
      f=f+a*(1.0d0+u*b+b*c+v*c)+d**2
      do 252 j=1,n
      if (j.le.n/2) g(j)=g(j)+2.0d0*(x(j)-1.0d0)
      if (j.le.n-2) then
      p=a*(u+c)*(x(j)+2.0d0*x(j+1)+1.0d1*x(j+2)-1.0d0)
      q=a*(v+b)*(2.0d0*x(j)+x(j+1)-3.0d0)
      g(j)=g(j)+2.0d0*p+4.0d0*q
      g(j+1)=g(j+1)+4.0d0*p+2.0d0*q
      g(j+2)=g(j+2)+2.0d1*p
      endif
      g(j)=g(j)+4.0d0*d*x(j)
  252 continue
      g(n-1)=g(n-1)+a*v*c
      g(n)=g(n)+a*u*b
      return
  310 do 311 ka=1,na
      if(mod(ka,2).eq.1) then
      fa=1.0d1*(x(ka+1)-x(ka)**2)
      g(ka)=g(ka)-2.0d1*x(ka)*fa
      g(ka+1)=g(ka+1)+1.0d1*fa
      else
      fa=1.0d0-x(ka-1)
      g(ka-1)=g(ka-1)-fa
      endif
      f=f+fa**2
  311 continue
      f=0.5d0*f
      return
  320 do 321 ka=1,na
      if(mod(ka,4).eq.1) then
      fa=x(ka)+1.0d1*x(ka+1)
      g(ka)=g(ka)+fa
      g(ka+1)=g(ka+1)+1.0d1*fa
      elseif(mod(ka,4).eq.2) then
      fa=2.23606797749979d0*(x(ka+1)-x(ka+2))
      g(ka+1)=g(ka+1)+2.23606797749979d0*fa
      g(ka+2)=g(ka+2)-2.23606797749979d0*fa
      elseif(mod(ka,4).eq.3) then
      a=x(ka-1)-2.0d0*x(ka)
      fa=a**2
      g(ka-1)=g(ka-1)+2.0d0*a*fa
      g(ka)=g(ka)-4.0d0*a*fa
      else
      fa=3.16227766016838d0*(x(ka-3)-x(ka))**2
      a=2.0d0*(x(ka-3)-x(ka))
      g(ka-3)=g(ka-3)+3.16227766016838d0*a*fa
      g(ka)=g(ka)-3.16227766016838d0*a*fa
      endif
      f=f+fa**2
  321 continue
      f=0.5d0*f
      return
  330 do 333 ka=1,na
      if(ka.le.n) then
      fa=(x(ka)-1.0d0)/3.16227766016838d0**5
      g(ka)=g(ka)+1.0d0/3.16227766016838d0**5*fa
      else
      fa=-0.25d0
      do 331 j=1,n
      fa=fa+x(j)**2
  331 continue
      do 332 j=1,n
      g(j)=g(j)+2.0d0*x(j)*fa
  332 continue
      endif
      f=f+fa**2
  333 continue
      f=0.5d0*f
      return
  350 do 354 ka=1,na
      if(ka.le.n) then
      fa=x(ka)-1.0d0
      g(ka)=g(ka)+fa
      else
      fa=0.0d0
      do 351 j=1,n
      fa=fa+dble(j)*(x(j)-1.0d0)
  351 continue
      if(ka.eq.n+1) then
      do 352 j=1,n
      g(j)=g(j)+dble(j)*fa
  352 continue
      else if(ka.eq.n+2) then
      do 353 j=1,n
      g(j)=g(j)+2.0d0*dble(j)*fa**3
  353 continue
      fa=fa**2
      endif
      endif
      f=f+fa**2
  354 continue
      f=0.5d0*f
      return
  370 do 376 ka=1,na
      if(ka.lt.n) then
      a=0.0d0
      do 371 j=1,n
      a=a+x(j)
  371 continue
      fa=x(ka)+a-dble(n+1)
      do 372 j=1,n
      g(j)=g(j)+fa
  372 continue
      g(ka)=g(ka)+fa
      else
      a=1.0d0
      do 373 j=1,n
      a=a*x(j)
  373 continue
      fa=a-1.0d0
      i=0
      do 374 j=1,n
      b=x(j)
      if(b.eq.0.0d0.and.i.eq.0) i=j
  374 continue
      if(i.ne.j) a=a*b
      if(i.eq.0) then
      do 375 j=1,n
      g(j)=g(j)+a/x(j)*fa
  375 continue
      else
      g(i)=g(i)+a*fa
      endif
      endif
      f=f+fa**2
  376 continue
      f=0.5d0*f
      return
  390 do 393 ka=1,na
      u=1.0d0/dble(n+1)
      v=dble(ka)*u
      a=0.0d0
      b=0.0d0
      do 391 j=1,n
      w=dble(j)*u
      if(j.le.ka) then
      a=a+w*(x(j)+w+1.0d0)**3
      else
      b=b+(1.0d0-w)*(x(j)+w+1.0d0)**3
      endif
  391 continue
      fa=x(ka)+u*((1.0d0-v)*a+v*b)/2.0d0
      f=f+fa**2
      do 392 j=1,n
      w=dble(j)*u
      a=(x(j)+w+1.0d0)**2
      if(j.le.ka) then
      g(j)=g(j)+1.5d0*u*(1.0d0-v)*w*a*fa
      else
      g(j)=g(j)+1.5d0*u*(1.0d0-w)*v*a*fa
      endif
  392 continue
      g(ka)=g(ka)+fa
  393 continue
      f=0.5d0*f
      return
  400 do 401 ka=1,na
      fa=(3.0d0-2.0d0*x(ka))*x(ka)+1.0d0
      if(ka.gt.1) fa=fa-x(ka-1)
      if(ka.lt.n) fa=fa-2.0d0*x(ka+1)
      f=f+fa**2
      g(ka)=g(ka)+(3.0d0-4.0d0*x(ka))*fa
      if(ka.gt.1) g(ka-1)=g(ka-1)-fa
      if(ka.lt.n) g(ka+1)=g(ka+1)-2.0d0*fa
  401 continue
      f=0.5d0*f
      return
  450 do 451 ka=1,na
      i=ka
      fa=(3.0d0-2.0d0*x(i))*x(i)+1.0d0
      if (i.gt.1) fa=fa-x(i-1)
      if (i.lt.n) fa=fa-x(i+1)
      f=f+fa**2
      g(i)=g(i)+(3.0d0-4.0d0*x(i))*fa
      if (i.gt.1) g(i-1)=g(i-1)-fa
      if (i.lt.n) g(i+1)=g(i+1)-fa
  451 continue
      f=0.5d0*f
      return
  460 do 463 ka=1,na
      i=ka
      fa=(2.0d0+5.0d0*x(i)**2)*x(i)+1.0d0
      do 461 j=max(1,i-5),min(n,i+1)
      fa=fa+x(j)*(1.0d0+x(j))
  461 continue
      f=f+fa**2
      do 462 j=max(1,i-5),min(n,i+1)
      g(j)=g(j)+(1.0d0+2.0d0*x(j))*fa
  462 continue
      g(i)=g(i)+(2.0d0+1.5d1*x(i)**2)*fa
  463 continue
      f=0.5d0*f
      return
  470 do 471 ka=1,na
      i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      fa=x(i)+x(i+1)*((5.0d0-x(i+1))*x(i+1)-2.0d0)-1.3d1
      g(i)=g(i)+fa
      g(i+1)=g(i+1)+(1.0d1*x(i+1)-3.0d0*x(i+1)**2-2.0d0)*fa
      else
      fa=x(i)+x(i+1)*((1.0d0+x(i+1))*x(i+1)-1.4d1)-2.9d1
      g(i)=g(i)+fa
      g(i+1)=g(i+1)+(2.0d0*x(i+1)+3.0d0*x(i+1)**2-1.4d1)*fa
      endif
      f=f+fa**2
  471 continue
      f=0.5d0*f
      return
  480 do 481 ka=1,na
      i=mod(ka,n/2)+1
      j=i+n/2
      m=5*n
      if (ka.le.m/2) then
      ia=1
      else
      ia=2
      endif
      ib=5-ka/(m/4)
      ic=mod(ka,5)+1
      fa=(x(i)**ia-x(j)**ib)**ic
      f=f+fa**2
      a=dble(ia)
      b=dble(ib)
      c=dble(ic)
      d=x(i)**ia-x(j)**ib
      if (d.ne.0.0d0) then
      e=c*d**(ic-1)
      if (x(i).eq.0.0d0.and.ia.le.1) then
      else
      g(i)=g(i)+e*a*x(i)**(ia-1)*fa
      endif
      if (x(j).eq.0.0d0.and.ib.le.1) then
      else
      g(j)=g(j)-e*b*x(j)**(ib-1)*fa
      endif
      endif
  481 continue
      f=0.5d0*f
      return
  490 do 491 ka=1,na
      i=2*((ka+5)/6)-1
      if (mod(ka,6).eq.1) then
      fa=x(i)+3.0d0*x(i+1)*(x(i+2)-1.0d0)+x(i+3)**2-1.0d0
      g(i)=g(i)+fa
      g(i+1)=g(i+1)+3.0d0*(x(i+2)-1.0d0)*fa
      g(i+2)=g(i+2)+3.0d0*x(i+1)*fa
      g(i+3)=g(i+3)+2.0d0*x(i+3)*fa
      elseif (mod(ka,6).eq.2) then
      fa=(x(i)+x(i+1))**2+(x(i+2)-1.0d0)**2-x(i+3)-3.0d0
      g(i)=g(i)+2.0d0*(x(i)+x(i+1))*fa
      g(i+1)=g(i+1)+2.0d0*(x(i)+x(i+1))*fa
      g(i+2)=g(i+2)+2.0d0*(x(i+2)-1.0d0)*fa
      g(i+3)=g(i+3)-fa
      elseif (mod(ka,6).eq.3) then
      fa=x(i)*x(i+1)-x(i+2)*x(i+3)
      g(i)=g(i)+x(i+1)*fa
      g(i+1)=g(i+1)+x(i)*fa
      g(i+2)=g(i+2)-x(i+3)*fa
      g(i+3)=g(i+3)-x(i+2)*fa
      elseif (mod(ka,6).eq.4) then
      fa=2.0d0*x(i)*x(i+2)+x(i+1)*x(i+3)-3.0d0
      g(i)=g(i)+2.0d0*x(i+2)*fa
      g(i+1)=g(i+1)+x(i+3)*fa
      g(i+2)=g(i+2)+2.0d0*x(i)*fa
      g(i+3)=g(i+3)+x(i+1)*fa
      elseif (mod(ka,6).eq.5) then
      fa=(x(i)+x(i+1)+x(i+2)+x(i+3))**2+(x(i)-1.0d0)**2
      g(i)=g(i)+(2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))+ &
       2.0d0*(x(i)-1.0d0))*fa
      g(i+1)=g(i+1)+2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))*fa
      g(i+2)=g(i+2)+2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))*fa
      g(i+3)=g(i+3)+2.0d0*(x(i)+x(i+1)+x(i+2)+x(i+3))*fa
      else
      fa=x(i)*x(i+1)*x(i+2)*x(i+3)+(x(i+3)-1.0d0)**2-1.0d0
      g(i)=g(i)+x(i+1)*x(i+2)*x(i+3)*fa
      g(i+1)=g(i+1)+x(i)*x(i+2)*x(i+3)*fa
      g(i+2)=g(i+2)+x(i)*x(i+1)*x(i+3)*fa
      g(i+3)=g(i+3)+(x(i)*x(i+1)*x(i+2)+2.0d0*(x(i+3)-1.0d0))*fa
      endif
      f=f+fa**2
  491 continue
      f=0.5d0*f
      return
  500 do 501 ka=1,na
      i=(ka+1)/2
      j=mod(ka,2)
      if (j.eq.0) then
      fa=6.0d0-exp(2.0d0*x(i))-exp(2.0d0*x(i+1))
      g(i)=g(i)-2.0d0*exp(2.0d0*x(i))*fa
      g(i+1)=g(i+1)-2.0d0*exp(2.0d0*x(i+1))*fa
      elseif (i.eq.1) then
      fa=4.0d0-exp(x(i))-exp(x(i+1))
      g(i)=g(i)-exp(x(i))*fa
      g(i+1)=g(i+1)-exp(x(i+1))*fa
      elseif (i.eq.n) then
      fa=8.0d0-exp(3.0d0*x(i-1))-exp(3.0d0*x(i))
      g(i-1)=g(i-1)-3.0d0*exp(3.0d0*x(i-1))*fa
      g(i)=g(i)-3.0d0*exp(3.0d0*x(i))*fa
      else
      fa=8.0d0-exp(3.0d0*x(i-1))-exp(3.0d0*x(i))+ &
       4.0d0-exp(x(i))-exp(x(i+1))
      g(i-1)=g(i-1)-3.0d0*exp(3.0d0*x(i-1))*fa
      g(i)=g(i)-(3.0d0*exp(3.0d0*x(i))+exp(x(i)))*fa
      g(i+1)=g(i+1)-exp(x(i+1))*fa
      endif
      f=f+fa**2
  501 continue
      f=0.5d0*f
      return
  510 do 511 ka=1,na
      i=(ka+1)/2
      if (mod(ka,2).eq.1) then
      fa=1.0d1*(2.0d0*x(i)/(1.0d0+x(i)**2)-x(i+1))
      g(i)=g(i)+2.0d1*(1.0d0-x(i)**2)/(1.0d0+x(i)**2)**2*fa
      g(i+1)=g(i+1)-1.0d1*fa
      else
      fa=x(i)-1.0d0
      g(i)=g(i)+fa
      endif
      f=f+fa**2
  511 continue
      f=0.5d0*f
      return
  520 do 521 ka=1,na
      i=3*((ka+5)/6)-2
      if (mod(ka,6).eq.1) then
      fa=1.0d1*(x(i)**2-x(i+1))
      g(i)=g(i)+2.0d1*x(i)*fa
      g(i+1)=g(i+1)-1.0d1*fa
      elseif (mod(ka,6).eq.2) then
      fa=x(i+2)-1.0d0
      g(i+2)=g(i+2)+fa
      elseif (mod(ka,6).eq.3) then
      fa=(x(i+3)-1.0d0)**2
      g(i+3)=g(i+3)+2.0d0*(x(i+3)-1.0d0)*fa
      elseif (mod(ka,6).eq.4) then
      fa=(x(i+4)-1.0d0)**3
      g(i+4)=g(i+4)+3.0d0*(x(i+4)-1.0d0)**2*fa
      elseif (mod(ka,6).eq.5) then
      fa=x(i)**2*x(i+3)+sin(x(i+3)-x(i+4))-1.0d1
      g(i)=g(i)+2.0d0*x(i)*x(i+3)*fa
      g(i+3)=g(i+3)+(x(i)**2+cos(x(i+3)-x(i+4)))*fa
      g(i+4)=g(i+4)-cos(x(i+3)-x(i+4))*fa
      else
      fa=x(i+1)+(x(i+2)**2*x(i+3))**2-2.0d1
      g(i+1)=g(i+1)+fa
      g(i+2)=g(i+2)+4.0d0*x(i+2)*(x(i+2)*x(i+3))**2*fa
      g(i+3)=g(i+3)+2.0d0*x(i+2)**4*x(i+3)*fa
      endif
      f=f+fa**2
  521 continue
      f=0.5d0*f
      return
  530 do 531 ka=1,na
      i=3*((ka+6)/7)-2
      if (mod(ka,7).eq.1) then
      fa=1.0d1*(x(i)**2-x(i+1))
      g(i)=g(i)+2.0d1*x(i)*fa
      g(i+1)=g(i+1)-1.0d1*fa
      elseif (mod(ka,7).eq.2) then
      fa=1.0d1*(x(i+1)**2-x(i+2))
      g(i+1)=g(i+1)+2.0d1*x(i+1)*fa
      g(i+2)=g(i+2)-1.0d1*fa
      elseif (mod(ka,7).eq.3) then
      fa=(x(i+2)-x(i+3))**2
      g(i+2)=g(i+2)+2.0d0*(x(i+2)-x(i+3))*fa
      g(i+3)=g(i+3)-2.0d0*(x(i+2)-x(i+3))*fa
      elseif (mod(ka,7).eq.4) then
      fa=(x(i+3)-x(i+4))**2
      g(i+3)=g(i+3)+2.0d0*(x(i+3)-x(i+4))*fa
      g(i+4)=g(i+4)-2.0d0*(x(i+3)-x(i+4))*fa
      elseif (mod(ka,7).eq.5) then
      fa=x(i)+x(i+1)**2+x(i+2)-3.0d1
      g(i)=g(i)+fa
      g(i+1)=g(i+1)+2.0d0*x(i+1)*fa
      g(i+2)=g(i+2)+fa
      elseif (mod(ka,7).eq.6) then
      fa=x(i+1)-x(i+2)**2+x(i+3)-1.0d1
      g(i+1)=g(i+1)+fa
      g(i+2)=g(i+2)-2.0d0*x(i+2)*fa
      g(i+3)=g(i+3)+fa
      else
      fa=x(i)*x(i+4)-1.0d1
      g(i)=g(i)+x(i+4)*fa
      g(i+4)=g(i+4)+x(i)*fa
      endif
      f=f+fa**2
  531 continue
      f=0.5d0*f
      return
  540 do 546 ka=1,na
      i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      fa=-y(l)
      do 542 k=1,3
      a=dble(k*k)/dble(l)
      do 541 j=1,4
      if (x(i+j).eq.0) x(i+j)=1.0d-16
      a=a*sign(1.0d0,x(i+j))*abs(x(i+j))**(dble(j)/dble(k*l))
  541 continue
      fa=fa+a
  542 continue
      f=f+fa**2
      do 545 k=1,3
      a=dble(k*k)/dble(l)
      do 543 j=1,4
      a=a*sign(1.0d0,x(i+j))*abs(x(i+j))**(dble(j)/dble(k*l))
  543 continue
      do 544 j=1,4
      g(i+j)=g(i+j)+(dble(j)/dble(k*l))*a/x(i+j)*fa
  544 continue
  545 continue
  546 continue
      f=0.5d0*f
      return
  550 do 556 ka=1,na
      i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      fa=-y(l)
      do 552 k=1,3
      a=0.0d0
      do 551 j=1,4
      a=a+x(i+j)*(dble(j)/dble(k*l))
  551 continue
      fa=fa+exp(a)*dble(k*k)/dble(l)
  552 continue
      f=f+fa**2
      do 555 k=1,3
      a=0.0d0
      do 553 j=1,4
      a=a+x(i+j)*(dble(j)/dble(k*l))
  553 continue
      a=exp(a)*dble(k*k)/dble(l)
      do 554 j=1,4
      g(i+j)=g(i+j)+a*(dble(j)/dble(k*l))*fa
  554 continue
  555 continue
  556 continue
      f=0.5d0*f
      return
  560 do 563 ka=1,na
      i=2*((ka+3)/4)-2
      l=mod((ka-1),4)+1
      fa=-y(l)
      do 561 j=1,4
      fa=fa+dble((1-2*mod(j,2))*l*j*j)*sin(x(i+j))+ &
       dble(l*l*j)*cos(x(i+j))
  561 continue
      f=f+fa**2
      do 562 j=1,4
      g(i+j)=g(i+j)+(dble((1-2*mod(j,2))*l*j*j)*cos(x(i+j))- &
       dble(l*l*j)*sin(x(i+j)))*fa
  562 continue
  563 continue
      f=0.5d0*f
      return
  570 do 571 ka=1,na
      alfa=0.5d0
      if (ka.eq.1) then
      fa=alfa-(1.0d0-alfa)*x(3)-x(1)*(1.0d0+4.0d0*x(2))
      g(1)=g(1)-(1.0d0+4.0d0*x(2))*fa
      g(2)=g(2)-4.0d0*x(1)*fa
      g(3)=g(3)+(alfa-1.0d0)*fa
      elseif(ka.eq.2) then
      fa=-(2.0d0-alfa)*x(4)-x(2)*(1.0d0+4.0d0*x(1))
      g(1)=g(1)-4.0d0*x(2)*fa
      g(2)=g(2)-(1.0d0+4.0d0*x(1))*fa
      g(4)=g(4)+(alfa-2.0d0)*fa
      elseif(ka.eq.n-1) then
      fa=alfa*x(n-3)-x(n-1)*(1.0d0+4.0d0*x(n))
      g(n-3)=g(n-3)+alfa*fa
      g(n-1)=g(n-1)-(1.0d0+4.0d0*x(n))*fa
      g(n)=g(n)-4.0d0*x(n-1)*fa
      elseif (ka.eq.n) then
      fa=alfa*x(n-2)-(2.0d0-alfa)-x(n)*(1.0d0+4.0d0*x(n-1))
      g(n-2)=g(n-2)+alfa*fa
      g(n-1)=g(n-1)-4.0d0*x(n)*fa
      g(n)=g(n)-(1.0d0+4.0d0*x(n-1))*fa
      elseif (mod(ka,2).eq.1) then
      fa=alfa*x(ka-2)-(1.0d0-alfa)*x(ka+2)- &
       x(ka)*(1.0d0+4.0d0*x(ka+1))
      g(ka-2)=g(ka-2)+alfa*fa
      g(ka)=g(ka)-(1.0d0+4.0d0*x(ka+1))*fa
      g(ka+1)=g(ka+1)-4.0d0*x(ka)*fa
      g(ka+2)=g(ka+2)+(alfa-1.0d0)*fa
      else
      fa=alfa*x(ka-2)-(2.0d0-alfa)*x(ka+2)- &
       x(ka)*(1.0d0+4.0d0*x(ka-1))
      g(ka-2)=g(ka-2)+alfa*fa
      g(ka-1)=g(ka-1)-4.0d0*x(ka)*fa
      g(ka)=g(ka)-(1.0d0+4.0d0*x(ka-1))*fa
      g(ka+2)=g(ka+2)+(alfa-2.0d0)*fa
      endif
      f=f+fa**2
  571 continue
      f=0.5d0*f
      return
  580 do 581 ka=1,na
      if (ka.lt.2) then
      fa=4.0d0*(x(ka)-x(ka+1)**2)
      g(ka)=g(ka)+4.0d0*fa
      g(ka+1)=g(ka+1)-8.0d0*x(ka+1)*fa
      elseif (ka.lt.n) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)
      g(ka-1)=g(ka-1)-8.0d0*x(ka)*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-8.0d0*x(ka+1)*fa
      else
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))
      g(ka-1)=g(ka-1)-8.0d0*x(ka)*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+2.0d0)*fa
      endif
      f=f+fa**2
  581 continue
      f=0.5d0*f
      return
  590 do 591 ka=1,na
      if (ka.eq.1) then
      fa=-2.0d0*x(ka)**2+3.0d0*x(ka)-2.0d0*x(ka+1)+ &
       3.0d0*x(n-4)-x(n-3)-x(n-2)+0.5d0*x(n-1)-x(n)+1.0d0
      g(n-4)=g(n-4)+3.0d0*fa
      g(n-3)=g(n-3)-fa
      g(n-2)=g(n-2)-fa
      g(n-1)=g(n-1)+0.50d0*fa
      g(n)=g(n)-fa
      g(ka)=g(ka)-(4.0d0*x(ka)-3.0d0)*fa
      g(ka+1)=g(ka+1)-2.0d0*fa
      elseif (ka.le.n-1) then
      fa=-2.0d0*x(ka)**2+3.0d0*x(ka)-x(ka-1)-2.0d0*x(ka+1)+ &
       3.0d0*x(n-4)-x(n-3)-x(n-2)+0.5d0*x(n-1)-x(n)+1.0d0
      g(n-4)=g(n-4)+3.0d0*fa
      g(n-3)=g(n-3)-fa
      g(n-2)=g(n-2)-fa
      g(n-1)=g(n-1)+0.50d0*fa
      g(n)=g(n)-fa
      g(ka-1)=g(ka-1)-fa
      g(ka)=g(ka)-(4.0d0*x(ka)-3.0d0)*fa
      g(ka+1)=g(ka+1)-2.0d0*fa
      else
      fa=-2.0d0*x(n)**2+3.0d0*x(n)-x(n-1)+ &
       3.0d0*x(n-4)-x(n-3)-x(n-2)+0.5d0*x(n-1)-x(n)+1.0d0
      g(n-4)=g(n-4)+3.0d0*fa
      g(n-3)=g(n-3)-fa
      g(n-2)=g(n-2)-fa
      g(n-1)=g(n-1)+0.50d0*fa
      g(n)=g(n)-(4.0d0*x(n)-2.0d0)*fa
      endif
      f=f+fa**2
  591 continue
      f=0.5d0*f
      return
  600 do 601 ka=1,na
      u=1.0d0/dble(n+1)
      v=dble(ka)*u
      fa=2.0d0*x(ka)+0.5d0*u*u*(x(ka)+v+1.0d0)**3+1.0d0
      if(ka.gt.1) fa=fa-x(ka-1)
      if(ka.lt.n) fa=fa-x(ka+1)
      f=f+fa**2
      g(ka)=g(ka)+(2.0d0+1.5d0*u**2*(x(ka)+v+1.0d0)**2)*fa
      if(ka.gt.1) g(ka-1)=g(ka-1)-fa
      if(ka.lt.n) g(ka+1)=g(ka+1)-fa
  601 continue
      f=0.5d0*f
      return
  610 do 611 ka=1,na
      i=3*((ka+6)/7)-2
      if (mod(ka,7).eq.1) then
      fa=1.0d1*(x(i)**2-x(i+1))
      g(i)=g(i)+2.0d1*x(i)*fa
      g(i+1)=g(i+1)-1.0d1*fa
      elseif (mod(ka,7).eq.2) then
      fa=x(i+1)+x(i+2)-2.0d0
      g(i+1)=g(i+1)+fa
      g(i+2)=g(i+2)+fa
      elseif (mod(ka,7).eq.3) then
      fa=x(i+3)-1.0d0
      g(i+3)=g(i+3)+fa
      elseif (mod(ka,7).eq.4) then
      fa=x(i+4)-1.0d0
      g(i+4)=g(i+4)+fa
      elseif (mod(ka,7).eq.5) then
      fa=x(i)+3.0d0*x(i+1)
      g(i)=g(i)+fa
      g(i+1)=g(i+1)+3.0d0*fa
      elseif (mod(ka,7).eq.6) then
      fa=x(i+2)+x(i+3)-2.0d0*x(i+4)
      g(i+2)=g(i+2)+fa
      g(i+3)=g(i+3)+fa
      g(i+4)=g(i+4)-2.0d0*fa
      else
      fa=1.0d1*(x(i+1)**2-x(i+4))
      g(i+1)=g(i+1)+2.0d1*x(i+1)*fa
      g(i+4)=g(i+4)-1.0d1*fa
      endif
      f=f+fa**2
  611 continue
      f=0.5d0*f
      return
  620 do 621 ka=1,na
      i=ka/2
      if (ka.eq.1) then
      fa=x(ka)-1.0d0
      g(ka)=g(ka)+fa
      else if (mod(ka,2).eq.0) then
      fa=1.0d1*(x(i)**2-x(i+1))
      g(i)=g(i)+2.0d1*x(i)*fa
      g(i+1)=g(i+1)-1.0d1*fa
      else
      a=2.0d0*exp(-(x(i)-x(i+1))**2)
      b=exp(-2.0d0*(x(i+1)-x(i+2))**2)
      fa=a+b
      g(i)=g(i)-2.0d0*(x(i)-x(i+1))*a*fa
      g(i+1)=g(i+1)+(2.0d0*(x(i)-x(i+1))*a-4.0d0*(x(i+1)-x(i+2))*b)*fa
      g(i+2)=g(i+2)+4.0d0*(x(i+1)-x(i+2))*b*fa
      endif
      f=f+fa**2
  621 continue
      f=0.5d0*f
      return
  630 do 633 ka=1,na
      ia=min(max(mod(ka,13)-2,1),7)
      ib=(ka+12)/13
      i=ia+ib-1
      if (ia.eq.7) then
      j=ib
      else
      j=ia+ib
      endif
      c=3.0d0*dble(ia)/1.0d1
      a=0.0d0
      do 631 l=0,6
      if (ib+l.ne.i.and.ib+l.ne.j) a=a+sin(x(ib+l))-y(1)
  631 continue
      fa=(1.0d0+cos(c))*(x(j)-sin(x(i))-1.0d0)+ &
       5.0d0*(x(i)-2.0d0)*exp(sin(c*x(j)))+0.5d0*a
      f=f+fa**2
      a=cos(c)
      b=exp(sin(c*x(j)))
      g(i)=g(i)-(cos(x(i))*(1.0d0+a)-5.0d0*b)*fa
      g(j)=g(j)+((1.0d0+a)+5.0d0*(x(i)-2.0d0)*c*cos(c*x(j))*b)*fa
      do 632 l=0,6
      if (ib+l.ne.i.and.ib+l.ne.j) &
       g(ib+l)=g(ib+l)+0.5d0*cos(x(ib+l))*fa
  632 continue
  633 continue
      f=0.5d0*f
      return
  720 do 721 ka=1,na
      a1=0.414214d0
      if (ka.eq.1) then
      fa=x(1)-(1.0d0-x(1))*x(3)-a1*(1.0d0+4.0d0*x(2))
      g(1)=g(1)+(1.0d0+x(3))*fa
      g(2)=g(2)-4.0d0*a1*fa
      g(3)=g(3)-(1.0d0-x(1))*fa
      elseif (ka.eq.2) then
      fa=-(1.0d0-x(1))*x(4)-a1*(1.0d0+4.0d0*x(2))
      g(1)=g(1)+x(4)*fa
      g(2)=g(2)-4.0d0*a1*fa
      g(4)=g(4)-(1.0d0-x(1))*fa
      elseif (ka.eq.3) then
      fa=a1*x(1)-(1.0d0-x(1))*x(5)-x(3)*(1.0d0+4.0d0*x(2))
      g(1)=g(1)+(a1+x(5))*fa
      g(2)=g(2)-4.0d0*x(3)*fa
      g(3)=g(3)-(1.0d0+4.0d0*x(2))*fa
      g(5)=g(5)-(1.0d0-x(1))*fa
      elseif (ka.le.n-2) then
      fa=x(1)*x(ka-2)-(1.0d0-x(1))*x(ka+2)- &
       x(ka)*(1.0d0+4.0d0*x(ka-1))
      g(1)=g(1)+(x(ka-2)+x(ka+2))*fa
      g(ka-2)=g(ka-2)+x(1)*fa
      g(ka-1)=g(ka-1)-4.0d0*x(ka)*fa
      g(ka)=g(ka)-(1.0d0+4.0d0*x(ka-1))*fa
      g(ka+2)=g(ka+2)-(1.0d0-x(1))*fa
      elseif (ka.eq.n-1) then
      fa=x(1)*x(n-3)-x(n-1)*(1.0d0+4.0d0*x(n-2))
      g(1)=g(1)+x(n-3)*fa
      g(n-3)=g(n-3)+x(1)*fa
      g(n-2)=g(n-2)-4.0d0*x(n-1)*fa
      g(n-1)=g(n-1)-(1.0d0+4.0d0*x(n-2))*fa
      else
      fa=x(1)*x(n-2)-(1.0d0-x(1))-x(n)*(1.0d0+4.0d0*x(n-1))
      g(1)=g(1)+(x(n-2)+1.0d0)*fa
      g(n-2)=g(n-2)+x(1)*fa
      g(n-1)=g(n-1)-4.0d0*x(n)*fa
      g(n)=g(n)-(1.0d0+4.0d0*x(n-1))*fa
      endif
      f=f+fa**2
  721 continue
      f=0.5d0*f
      return
  740 do 741 ka=1,na
      if (ka.lt.2) then
      fa=3.0d0*x(ka)**3+2.0d0*x(ka+1)-5.0d0+ &
       sin(x(ka)-x(ka+1))*sin(x(ka)+x(ka+1))
      d1s=cos(x(ka)-x(ka+1))*sin(x(ka)+x(ka+1))
      d2s=sin(x(ka)-x(ka+1))*cos(x(ka)+x(ka+1))
      g(ka)=g(ka)+(9.0d0*x(ka)**2+d1s+d2s)*fa
      g(ka+1)=g(ka+1)+(2.0d0-d1s+d2s)*fa
      elseif (ka.lt.n) then
      fa=3.0d0*x(ka)**3+2.0d0*x(ka+1)-5.0d0+ &
       sin(x(ka)-x(ka+1))*sin(x(ka)+x(ka+1))+4.0d0*x(ka)- &
       x(ka-1)*exp(x(ka-1)-x(ka))-3.0d0
      d1s=cos(x(ka)-x(ka+1))*sin(x(ka)+x(ka+1))
      d2s=sin(x(ka)-x(ka+1))*cos(x(ka)+x(ka+1))
      ex=exp(x(ka-1)-x(ka))
      g(ka-1)=g(ka-1)-(ex+x(ka-1)*ex)*fa
      g(ka)=g(ka)+(9.0d0*x(ka)**2+d1s+d2s+4.0d0+x(ka-1)*ex)*fa
      g(ka+1)=g(ka+1)+(2.0d0-d1s+d2s)*fa
      else
      fa=4.0d0*x(ka)-x(ka-1)*exp(x(ka-1)-x(ka))-3.0d0
      ex=exp(x(ka-1)-x(ka))
      g(ka-1)=g(ka-1)-(ex+x(ka-1)*ex)*fa
      g(ka)=g(ka)+(4.0d0+x(ka-1)*ex)*fa
      endif
      f=f+fa**2
  741 continue
      f=0.5d0*f
      return
  750 do 751 ka=1,na
      if (mod(ka,2).eq.1) then
      fa=0.0d0
      if (ka.ne.1) fa=fa-6.0d0*(x(ka-2)-x(ka))**3+1.0d1- &
       4.0d0*x(ka-1)-2.0d0*sin(x(ka-2)-x(ka-1)-x(ka))* &
       sin(x(ka-2)+x(ka-1)-x(ka))
      if (ka.ne.n) fa=fa+3.0d0*(x(ka)-x(ka+2))**3-5.0d0+ &
       2.0d0*x(ka+1)+sin(x(ka)-x(ka+1)-x(ka+2))* &
       sin(x(ka)+x(ka+1)-x(ka+2))
      if (ka.ne.1) then
      d1s=cos(x(ka-2)-x(ka-1)-x(ka))*sin(x(ka-2)+x(ka-1)-x(ka))
      d2s=sin(x(ka-2)-x(ka-1)-x(ka))*cos(x(ka-2)+x(ka-1)-x(ka))
      g(ka-2)=g(ka-2)-(18.0d0*(x(ka-2)-x(ka))**2+2.0d0*(d1s+d2s))*fa
      g(ka-1)=g(ka-1)-(4.0d0-2.0d0*(d1s-d2s))*fa
      g(ka)=g(ka)+(18.0d0*(x(ka-2)-x(ka))**2+2.0d0*(d1s+d2s))*fa
      endif
      if (ka.ne.n) then
      d1s=cos(x(ka)-x(ka+1)-x(ka+2))*sin(x(ka)+x(ka+1)-x(ka+2))
      d2s=sin(x(ka)-x(ka+1)-x(ka+2))*cos(x(ka)+x(ka+1)-x(ka+2))
      g(ka)=g(ka)+(9.0d0*(x(ka)-x(ka+2))**2+d1s+d2s)*fa
      g(ka+1)=g(ka+1)+(2.0d0-d1s+d2s)*fa
      g(ka+2)=g(ka+2)-(9.0d0*(x(ka)-x(ka+2))**2+d1s+d2s)*fa
      endif
      else
      ex=exp(x(ka-1)-x(ka)-x(ka+1))
      fa=4.0d0*x(ka)-(x(ka-1)-x(ka+1))*ex-3.0d0
      w=x(ka-1)-x(ka+1)
      g(ka-1)=g(ka-1)-(ex+w*ex)*fa
      g(ka)=g(ka)+(4.0d0+w*ex)*fa
      g(ka+1)=g(ka+1)+(ex+w*ex)*fa
      endif
      f=f+fa**2
  751 continue
      f=0.5d0*f
      return
  760 do 761 ka=1,na
      h=2.0d0
      if (ka.eq.1) then
      fa=((3.0d0-h*x(1))*x(1)-2.0d0*x(2)+1.0d0)**2
      g(1)=g(1)+2.0d0*((3.0d0-h*x(1))*x(1)-2.0d0*x(2)+1.0d0)* &
       (3.0d0-2.0d0*h*x(1))*fa
      g(2)=g(2)-4.0d0*((3.0d0-h*x(1))*x(1)-2.0d0*x(2)+1.0d0)*fa
      elseif (ka.le.n-1) then
      fa=((3.0d0-h*x(ka))*x(ka)-x(ka-1)-2.0d0*x(ka+1)+1.0d0)**2
      g(ka-1)=g(ka-1)-2.0d0* &
       ((3.0d0-h*x(ka))*x(ka)-x(ka-1)-2.0d0*x(ka+1)+1.0d0)*fa
      g(ka)=g(ka)+2.0d0* &
       ((3.0d0-h*x(ka))*x(ka)-x(ka-1)-2.0d0*x(ka+1)+1.0d0)* &
       (3.0d0-2.0d0*h*x(ka))*fa
      g(ka+1)=g(ka+1)-4.0d0* &
       ((3.0d0-h*x(ka))*x(ka)-x(ka-1)-2.0d0*x(ka+1)+1.0d0)*fa
      else
      fa=((3.0d0-h*x(n))*x(n)-x(n-1)+1.0d0)**2
      g(n-1)=g(n-1)-2.0d0*((3.0d0-h*x(n))*x(n)-x(n-1)+1.0d0)*fa
      g(n)=g(n)+2.0d0*((3.0d0-h*x(n))*x(n)-x(n-1)+1.0d0)* &
       (3.0d0-2.0d0*h*x(n))*fa
      endif
      f=f+fa**2
  761 continue
      f=0.5d0*f
      return
  780 do 781 ka=1,na
      if (ka.lt.2) then
      fa=4.0d0*(x(ka)-x(ka+1)**2)+x(ka+1)-x(ka+2)**2
      g(ka)=g(ka)+4.0d0*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-2.0d0*x(ka+2)*fa
      elseif (ka.lt.3) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka+1)-x(ka+2)**2
      g(ka-1)=g(ka-1)-8.0d0*x(ka)*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-2.0d0*x(ka+2)*fa
      elseif (ka.lt.n-1) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka-1)**2-x(ka-2)+x(ka+1)- &
       x(ka+2)**2
      g(ka-2)=g(ka-2)-fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-2.0d0*x(ka+2)*fa
      elseif (ka.lt.n) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka-1)**2-x(ka-2)
      g(ka-2)=g(ka-2)-fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-8.0d0*x(ka+1)*fa
      else
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka)) &
       +x(ka-1)**2-x(ka-2)
      g(ka-2)=g(ka-2)-fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+2.0d0)*fa
      endif
      f=f+fa**2
  781 continue
      f=0.5d0*f
      return
  790 do 791 ka=1,na
      if (ka.lt.2) then
      fa=4.0d0*(x(ka)-x(ka+1)**2)+x(ka+1)-x(ka+2)**2+ &
       x(ka+2)-x(ka+3)**2
      g(ka)=g(ka)+4.0d0*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-(2.0d0*x(ka+2)-1.0d0)*fa
      g(ka+3)=g(ka+3)-2.0d0*x(ka+3)*fa
      elseif (ka.lt.3) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka-1)**2+x(ka+1)-x(ka+2)**2+ &
       x(ka+2)-x(ka+3)**2
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-(2.0d0*x(ka+2)-1.0d0)*fa
      g(ka+3)=g(ka+3)-2.0d0*x(ka+3)*fa
      elseif (ka.lt.4) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka-1)**2-x(ka-2)+ &
       x(ka+1)-x(ka+2)**2+x(ka-2)**2+x(ka+2)-x(ka+3)**2
      g(ka-2)=g(ka-2)+(2.0d0*x(ka-2)-1.0d0)*fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-(2.0d0*x(ka+2)-1.0d0)*fa
      g(ka+3)=g(ka+3)-2.0d0*x(ka+3)*fa
      elseif (ka.lt.n-2) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka-1)**2-x(ka-2)+ &
       x(ka+1)-x(ka+2)**2+x(ka-2)**2-x(ka-3)+x(ka+2)-x(ka+3)**2
      g(ka-3)=g(ka-3)-fa
      g(ka-2)=g(ka-2)+(2.0d0*x(ka-2)-1.0d0)*fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-(2.0d0*x(ka+2)-1.0d0)*fa
      g(ka+3)=g(ka+3)-2.0d0*x(ka+3)*fa
      elseif (ka.lt.n-1) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka-1)**2-x(ka-2)+ &
       x(ka+1)-x(ka+2)**2+x(ka-2)**2-x(ka-3)+x(ka+2)
      g(ka-3)=g(ka-3)-fa
      g(ka-2)=g(ka-2)+(2.0d0*x(ka-2)-1.0d0)*fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      g(ka+2)=g(ka+2)-(2.0d0*x(ka+2)-1.0d0)*fa
      elseif (ka.lt.n) then
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka))+ &
       4.0d0*(x(ka)-x(ka+1)**2)+x(ka-1)**2-x(ka-2)+ &
       x(ka+1)+x(ka-2)**2-x(ka-3)
      g(ka-3)=g(ka-3)-fa
      g(ka-2)=g(ka-2)+(2.0d0*x(ka-2)-1.0d0)*fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+6.0d0)*fa
      g(ka+1)=g(ka+1)-(8.0d0*x(ka+1)-1.0d0)*fa
      else
      fa=8.0d0*x(ka)*(x(ka)**2-x(ka-1))-2.0d0*(1.0d0-x(ka)) &
       +x(ka-1)**2-x(ka-2)+x(ka-2)**2-x(ka-3)
      g(ka-3)=g(ka-3)-fa
      g(ka-2)=g(ka-2)+(2.0d0*x(ka-2)-1.0d0)*fa
      g(ka-1)=g(ka-1)-(8.0d0*x(ka)-2.0d0*x(ka-1))*fa
      g(ka)=g(ka)+(24.0d0*x(ka)**2-8.0d0*x(ka-1)+2.0d0)*fa
      endif
      f=f+fa**2
  791 continue
      f=0.5d0*f
      return
  810 do 811 ka=1,na
      if(mod(ka,2).eq.1) then
      fa=x(ka)+((5.0d0-x(ka+1))*x(ka+1)-2.0d0)*x(ka+1)-1.3d1
      g(ka)=g(ka)+fa
      g(ka+1)=g(ka+1)+(10.0d0*x(ka+1)-3.0d0*x(ka+1)**2-2.0d0)*fa
      else
      fa=x(ka-1)+((x(ka)+1.0d0)*x(ka)-1.4d1)*x(ka)-2.9d1
      g(ka-1)=g(ka-1)+fa
      g(ka)=g(ka)+(3.0d0*x(ka)**2+2.0d0*x(ka)-1.4d1)*fa
      endif
      f=f+fa**2
  811 continue
      f=0.5d0*f
      return
  830 do 831 ka=1,na
      if(mod(ka,4).eq.1) then
      a=exp(x(ka))-x(ka+1)
      fa=a**2
      g(ka)=g(ka)+2.0d0*a*exp(x(ka))*fa
      g(ka+1)=g(ka+1)-2.0d0*a*fa
      elseif (mod(ka,4).eq.2) then
      fa=1.0d1*(x(ka)-x(ka+1))**3
      a=3.0d1*(x(ka)-x(ka+1))**2*fa
      g(ka)=g(ka)+a
      g(ka+1)=g(ka+1)-a
      elseif (mod(ka,4).eq.3) then
      a=x(ka)-x(ka+1)
      fa=(sin(a)/cos(a))**2
      b=2.0d0*sin(a)/(cos(a))**3*fa
      g(ka)=g(ka)+b
      g(ka+1)=g(ka+1)-b
      else
      fa=x(ka)-1.0d0
      g(ka)=g(ka)+fa
      endif
      f=f+fa**2
  831 continue
      f=0.5d0*f
      return
  840 do 841 ka=1,na
      if(ka.lt.2) then
      fa=x(ka)*(0.5d0*x(ka)-3.0d0)-1.0d0+2.0d0*x(ka+1)
      g(ka)=g(ka)+(x(ka)-3.0d0)*fa
      g(ka+1)=g(ka+1)+2.0d0*fa
      elseif (ka.lt.n) then
      fa=x(ka-1)+x(ka)*(0.5d0*x(ka)-3.0d0)-1.0d0+2.0d0*x(ka+1)
      g(ka-1)=g(ka-1)+fa
      g(ka)=g(ka)+(x(ka)-3.0d0)*fa
      g(ka+1)=g(ka+1)+2.0d0*fa
      else
      fa=x(ka-1)+x(ka)*(0.5d0*x(ka)-3.0d0)-1.0d0
      g(ka-1)=g(ka-1)+fa
      g(ka)=g(ka)+(x(ka)-3.0d0)*fa
      endif
      f=f+fa**2
  841 continue
      f=0.5d0*f
      return
  860 do 861 ka=1,na
      if(mod(ka,2).eq.1) then
      fa=1.0d4*x(ka)*x(ka+1)-1.0d0
      g(ka)=g(ka)+1.0d4*x(ka+1)*fa
      g(ka+1)=g(ka+1)+1.0d4*x(ka)*fa
      else
      fa=exp(-x(ka-1))+exp(-x(ka))-1.0001d0
      g(ka-1)=g(ka-1)-exp(-x(ka-1))*fa
      g(ka)=g(ka)-exp(-x(ka))*fa
      endif
      f=f+fa**2
  861 continue
      f=0.5d0*f
      return
  870 do 871 ka=1,na
      if(mod(ka,4).eq.1) then
      fa=-2.0d2*x(ka)*(x(ka+1)-x(ka)**2)-(1.0d0-x(ka))
      g(ka)=g(ka)-(2.0d2*(x(ka+1)-3.0d0*x(ka)**2)-1.0d0)*fa
      g(ka+1)=g(ka+1)-2.0d2*x(ka)*fa
      elseif(mod(ka,4).eq.2) then
      fa=2.0d2*(x(ka)-x(ka-1)**2)+2.02d1*(x(ka)-1.0d0)+ &
        1.98d1*(x(ka+2)-1.0d0)
      g(ka-1)=g(ka-1)-4.0d2*x(ka-1)*fa
      g(ka)=g(ka)+2.202d2*fa
      g(ka+2)=g(ka+2)+1.98d1*fa
      elseif(mod(ka,4).eq.3) then
      fa=-1.8d2*x(ka)*(x(ka+1)-x(ka)**2)-(1.0d0-x(ka))
      g(ka)=g(ka)-(1.8d2*(x(ka+1)-3.0d0*x(ka)**2)-1.0d0)*fa
      g(ka+1)=g(ka+1)-1.8d2*x(ka)*fa
      else
      fa=1.8d2*(x(ka)-x(ka-1)**2)+2.02d1*(x(ka)-1.0d0)+ &
       1.98d1*(x(ka-2)-1.0d0)
      g(ka-2)=g(ka-2)+1.98d1*fa
      g(ka-1)=g(ka-1)-3.6d2*x(ka-1)*fa
      g(ka)=g(ka)+2.002d2*fa
      endif
      f=f+fa**2
  871 continue
      f=0.5d0*f
      return
  880 do 881 ka=1,na
      if (ka.lt.2) then
      a=exp(cos(dble(ka)*(x(ka)+x(ka+1))))
      b=a*dble(ka)*sin(dble(ka)*(x(ka)+x(ka+1)))
      fa=x(ka)-a
      g(ka+1)=g(ka+1)+b*fa
      g(ka)=g(ka)+(b+1.0d0)*fa
      elseif (ka.lt.n) then
      a=exp(cos(dble(ka)*(x(ka-1)+x(ka)+x(ka+1))))
      b=a*sin(dble(ka)*(x(ka-1)+x(ka)+x(ka+1)))*dble(ka)
      fa=x(ka)-a
      g(ka-1)=g(ka-1)+b*fa
      g(ka+1)=g(ka+1)+b*fa
      g(ka)=g(ka)+(b+1.0d0)*fa
      else
      a=exp(cos(dble(ka)*(x(ka-1)+x(ka))))
      b=a*sin(dble(ka)*(x(ka-1)+x(ka)))*dble(ka)
      fa=x(ka)-a
      g(ka-1)=g(ka-1)+b*fa
      g(ka)=g(ka)+(b+1.0d0)*fa
      endif
      f=f+fa**2
  881 continue
      f=0.5d0*f
      return
  900 do 901 ka=1,na
      if(ka.eq.1) then
      fa=3.0d0*x(ka)*(x(ka+1)-2.0d0*x(ka))+0.25d0*x(ka+1)**2
      g(ka)=g(ka)+3.0d0*(x(ka+1)-4.0d0*x(ka))*fa
      g(ka+1)=g(ka+1)+(3.0d0*x(ka)+0.5d0*x(ka+1))*fa
      elseif(ka.eq.n) then
      fa=3.0d0*x(ka)*(2.0d1-2.0d0*x(ka)+x(ka-1))+ &
       0.25d0*(2.0d1-x(ka-1))**2
      g(ka-1)=g(ka-1)+(3.0d0*x(ka)-0.5d0*(2.0d1-x(ka-1)))*fa
      g(ka)=g(ka)+3.0d0*(2.0d1-4.0d0*x(ka)+x(ka-1))*fa
      else
      fa=3.0d0*x(ka)*(x(ka+1)-2.0d0*x(ka)+x(ka-1))+ &
       0.25d0*(x(ka+1)-x(ka-1))**2
      g(ka-1)=g(ka-1)+(3.0d0*x(ka)-0.5d0*(x(ka+1)-x(ka-1)))*fa
      g(ka)=g(ka)+3.0d0*(x(ka+1)-4.0d0*x(ka)+x(ka-1))*fa
      g(ka+1)=g(ka+1)+(3.0d0*x(ka)+0.5d0*(x(ka+1)-x(ka-1)))*fa
      endif
      f=f+fa**2
  901 continue
      f=0.5d0*f
      return
  910 do 911 ka=1,na
      h=1.0d0/dble(n+1)
      if (ka.lt.2) then
      fa=2.0d0*x(ka)+par*h**2*sinh(par*x(ka))-x(ka+1)
      g(ka)=g(ka)+(2.0d0+par**2*h**2*cosh(par*x(ka)))*fa
      g(ka+1)=g(ka+1)-fa
      else if (ka.lt.n) then
      fa=2.0d0*x(ka)+par*h**2*sinh(par*x(ka))-x(ka-1)-x(ka+1)
      g(ka-1)=g(ka-1)-fa
      g(ka)=g(ka)+(2.0d0+par**2*h**2*cosh(par*x(ka)))*fa
      g(ka+1)=g(ka+1)-fa
      else
      fa=2.0d0*x(ka)+par*h**2*sinh(par*x(ka))-x(ka-1)-1.0d0
      g(ka)=g(ka)+(2.0d0+par**2*h**2*cosh(par*x(ka)))*fa
      g(ka-1)=g(ka-1)-fa
      endif
      f=f+fa**2
  911 continue
      f=0.5d0*f
      return
  920 do 921 ka=1,na
      fa=6.0d0*x(ka)
      a1=0.0d0
      a2=0.0d0
      a3=0.0d0
      if (ka.gt.1) then
      fa=fa-4.0d0*x(ka-1)
      a1=a1-x(ka-1)
      a2=a2+x(ka-1)
      a3=a3+2.0d0*x(ka-1)
      endif
      if (ka.gt.2) then
      fa=fa+x(ka-2)
      a3=a3-x(ka-2)
      endif
      if (ka.lt.n-1) then
      fa=fa+x(ka+2)
      a3=a3+x(ka+2)
      endif
      if (ka.lt.n) then
      fa=fa-4.0d0*x(ka+1)
      a1=a1+x(ka+1)
      a2=a2+x(ka+1)
      a3=a3-2.0d0*x(ka+1)
      endif
      if (ka.ge.n-1) then
      fa=fa+1.0d0
      a3=a3+1.0d0
      endif
      if (ka.ge.n) then
      fa=fa-4.0d0
      a1=a1+1.0d0
      a2=a2+1.0d0
      a3=a3-2.0d0
      endif
      fa=fa-0.5d0*par*(a1*a2-x(ka)*a3)
      f=f+fa**2
      g(ka)=g(ka)+6.0d0*fa
      ga1(1)=0.0d0
      ga1(2)=0.0d0
      ga2(1)=0.0d0
      ga2(2)=0.0d0
      if (ka.gt.1) then
      g(ka-1)=g(ka-1)-(4.0d0-par*x(ka))*fa
      ga1(1)=-1.0d0
      ga2(1)= 1.0d0
      endif
      if (ka.gt.2) then
      g(ka-2)=g(ka-2)+(1.0d0-0.5d0*par*x(ka))*fa
      endif
      if (ka.lt.n-1) then
      g(ka+2)=g(ka+2)+(1.0d0+0.5d0*par*x(ka))*fa
      endif
      if (ka.lt.n) then
      g(ka+1)=g(ka+1)-(4.0d0+par*x(ka))*fa
      ga1(2)= 1.0d0
      ga2(2)= 1.0d0
      endif
      g(ka)=g(ka)+0.5d0*par*a3*fa
      if (ka.gt.1) &
       g(ka-1)=g(ka-1)-0.5d0*par*(ga1(1)*a2+a1*ga2(1))*fa
      if (ka.lt.n) &
       g(ka+1)=g(ka+1)-0.5d0*par*(ga1(2)*a2+a1*ga2(2))*fa
  921 continue
      f=0.5d0*f
      return
  930 do 931 ka=1,na
      h=1.0d0/dble(m+1)
      if(ka.le.m) then
      j=ka+m
      fa=6.0d0*x(ka)
      a1=0.0d0
      a2=0.0d0
      if (ka.eq.1) then
      a1=a1+1.0d0
      endif
      if (ka.gt.1) then
      fa=fa-4.0d0*x(ka-1)
      a1=a1-x(j-1)
      a2=a2+2.0d0*x(ka-1)
      endif
      if (ka.gt.2) then
      fa=fa+x(ka-2)
      a2=a2-x(ka-2)
      endif
      if (ka.lt.m-1) then
      fa=fa+x(ka+2)
      a2=a2+x(ka+2)
      endif
      if (ka.lt.m) then
      fa=fa-4.0d0*x(ka+1)
      a1=a1+x(j+1)
      a2=a2-2.0d0*x(ka+1)
      endif
      if (ka.eq.m) then
      a1=a1+1.0d0
      endif
      fa=fa+0.5d0*par*h*(x(ka)*a2+x(j)*a1*h**2)
      else
      j=ka-m
      fa=-2.0d0*x(ka)
      a1=0.0d0
      a2=0.0d0
      if (j.eq.1) then
      a2=a2+1.0d0
      endif
      if (j.gt.1) then
      fa=fa+x(ka-1)
      a1=a1-x(j-1)
      a2=a2-x(ka-1)
      endif
      if (j.lt.m) then
      fa=fa+x(ka+1)
      a1=a1+x(j+1)
      a2=a2+x(ka+1)
      endif
      if (j.eq.m) then
      a2=a2+1.0d0
      endif
      fa=fa+0.5d0*par*h*(x(ka)*a1+x(j)*a2)
      endif
      f=f+fa**2
      if(ka.le.m) then
      g(ka)=g(ka)+6.0d0*fa
      if (ka.gt.1) then
      g(ka-1)=g(ka-1)-(4.0d0-par*h*x(ka))*fa
      g(j-1)=g(j-1)-0.5d0*par*h**3*x(j)*fa
      endif
      if (ka.gt.2) then
      g(ka-2)=g(ka-2)+(1.0d0-0.5d0*par*h*x(ka))*fa
      endif
      if (ka.lt.m-1) then
      g(ka+2)=g(ka+2)+(1.0d0+0.5d0*par*h*x(ka))*fa
      endif
      if (ka.lt.m) then
      g(ka+1)=g(ka+1)-(4.0d0+par*h*x(ka))*fa
      g(j+1)=g(j+1)+0.5d0*par*h**3*x(j)*fa
      endif
      g(ka)=g(ka)+0.5d0*par*h*a2*fa
      g(j)=g(j)+0.5d0*par*h**3*a1*fa
      else
      g(ka)=g(ka)-2.0d0*fa
      if (j.gt.1) then
      g(ka-1)=g(ka-1)+(1.0d0-0.5d0*par*h*x(j))*fa
      g(j-1)=g(j-1)-0.5d0*par*h*x(ka)*fa
      endif
      if (j.lt.m) then
      g(ka+1)=g(ka+1)+(1.0d0+0.5d0*par*h*x(j))*fa
      g(j+1)=g(j+1)+0.5d0*par*h*x(ka)*fa
      endif
      g(ka)=g(ka)+0.5d0*par*h*a1*fa
      g(j)=g(j)+0.5d0*par*h*a2*fa
      endif
  931 continue
      f=0.5d0*f
      return
  940 do 941 ka=1,na
      fa=4.0d0*x(ka)-par*exp(x(ka))
      j=(ka-1)/m+1
      i=ka-(j-1)*m
      if(i.gt.1) fa=fa-x(ka-1)
      if(i.lt.m) fa=fa-x(ka+1)
      if(j.gt.1) fa=fa-x(ka-m)
      if(j.lt.m) fa=fa-x(ka+m)
      f=f+fa**2
      g(ka)=g(ka)+(4.0d0-par*exp(x(ka)))*fa
      if(j.gt.1) g(ka-m)=g(ka-m)-fa
      if(i.gt.1) g(ka-1)=g(ka-1)-fa
      if(i.lt.m) g(ka+1)=g(ka+1)-fa
      if(j.lt.m) g(ka+m)=g(ka+m)-fa
  941 continue
      f=0.5d0*f
      return
  950 do 951 ka=1,na
      j=(ka-1)/m+1
      i=ka-(j-1)*m
      fa=4.0d0*x(ka)+par*x(ka)**3/(1.0d0+par*dble(i)**2+ &
       par*dble(j)**2)
      if(i.eq.1) fa=fa-1.0d0
      if(i.gt.1) fa=fa-x(ka-1)
      if(i.lt.m) fa=fa-x(ka+1)
      if(i.eq.m) fa=fa-2.0d0+exp(dble(j)/dble(m+1))
      if(j.eq.1) fa=fa-1.0d0
      if(j.gt.1) fa=fa-x(ka-m)
      if(j.lt.m) fa=fa-x(ka+m)
      if(j.eq.m) fa=fa-2.0d0+exp(dble(i)/dble(m+1))
      f=f+fa**2
      g(ka)=g(ka)+(4.0d0+3.0d0*par*x(ka)**2/(1.0d0+par*dble(i)**2+ &
       par*dble(j)**2))*fa
      if(j.gt.1) g(ka-m)=g(ka-m)-fa
      if(i.gt.1) g(ka-1)=g(ka-1)-fa
      if(i.lt.m) g(ka+1)=g(ka+1)-fa
      if(j.lt.m) g(ka+m)=g(ka+m)-fa
  951 continue
      f=0.5d0*f
      return
  960 do 961 ka=1,na
      j=(ka-1)/m+1
      i=ka-(j-1)*m
      a1=dble(i)/dble(m+1)
      a2=dble(j)/dble(m+1)
      fa=4.0d0*x(ka)-par*sin(2.0d0*pi*x(ka))- &
       1.0d4*((a1-0.25d0)**2+(a2-0.75d0)**2)*par
      if(i.eq.1) fa=fa-x(ka+1)-par*sin(pi*x(ka+1)*dble(m+1))
      if(i.gt.1.and.i.lt.m) fa=fa-x(ka+1)-x(ka-1)- &
       par*sin(pi*(x(ka+1)-x(ka-1))*dble(m+1))
      if(i.eq.m) fa=fa-x(ka-1)+par*sin(pi*x(ka-1)*dble(m+1))
      if(j.eq.1) fa=fa-x(ka+m)-par*sin(pi*x(ka+m)*dble(m+1))
      if(j.gt.1.and.j.lt.m) fa=fa-x(ka+m)-x(ka-m)- &
       par*sin(pi*(x(ka+m)-x(ka-m))*dble(m+1))
      if(j.eq.m) fa=fa-x(ka-m)+par*sin(pi*x(ka-m)*dble(m+1))
      f=f+fa**2
      g(ka)=g(ka)+(4.0d0-2.0d0*pi*par*cos(2.0d0*pi*x(ka)))*fa
      if(i.eq.1) g(ka+1)=g(ka+1)- &
       (1.0d0+pi*dble(m+1)*par*cos(pi*x(ka+1)*dble(m+1)))*fa
      if(i.gt.1.and.i.lt.m) then
      g(ka-1)=g(ka-1)- &
       (1.0d0-pi*dble(m+1)*par*cos(pi*(x(ka+1)-x(ka-1))*dble(m+1)))*fa
      g(ka+1)=g(ka+1)- &
       (1.0d0+pi*dble(m+1)*par*cos(pi*(x(ka+1)-x(ka-1))*dble(m+1)))*fa
      endif
      if(i.eq.m) g(ka-1)=g(ka-1)- &
       (1.0d0-pi*dble(m+1)*par*cos(pi*x(ka-1)*dble(m+1)))*fa
      if(j.eq.1) g(ka+m)=g(ka+m)- &
       (1.0d0+pi*dble(m+1)*par*cos(pi*x(ka+m)*dble(m+1)))*fa
      if(j.gt.1.and.j.lt.m) then
      g(ka-m)=g(ka-m)- &
       (1.0d0-pi*dble(m+1)*par*cos(pi*(x(ka+m)-x(ka-m))*dble(m+1)))*fa
      g(ka+m)=g(ka+m)- &
       (1.0d0+pi*dble(m+1)*par*cos(pi*(x(ka+m)-x(ka-m))*dble(m+1)))*fa
      endif
      if(j.eq.m) g(ka-m)=g(ka-m)- &
       (1.0d0-pi*dble(m+1)*par*cos(pi*x(ka-m)*dble(m+1)))*fa
  961 continue
      f=0.5d0*f
      return
  970 do 971 ka=1,na
      j=(ka-1)/m+1
      i=ka-(j-1)*m
      fa=8.0d0*x(ka)**2
      if(i.eq.1) fa=fa-2.0d0*x(ka)*(x(ka+1)+1.0d0)- &
       0.5d0*(x(ka+1)-1.0d0)**2- &
       1.5d0*x(ka)**2*(x(ka+1)-1.0d0)*par
      if(i.gt.1.and.i.lt.m) fa=fa-2.0d0*x(ka)*(x(ka+1)+x(ka-1))- &
       0.5d0*(x(ka+1)-x(ka-1))**2- &
       1.5d0*x(ka)**2*(x(ka+1)-x(ka-1))*par
      if(i.eq.m) fa=fa-2.0d0*x(ka)*x(ka-1)- &
       0.5d0*x(ka-1)**2+ &
       1.5d0*x(ka)**2*x(ka-1)*par
      if(j.eq.1) fa=fa-2.0d0*x(ka)*(x(ka+m)+1.0d0)- &
       0.5d0*(x(ka+m)-1.0d0)**2
      if(j.gt.1.and.j.lt.m) fa=fa-2.0d0*x(ka)*(x(ka+m)+x(ka-m))- &
       0.5d0*(x(ka+m)-x(ka-m))**2
      if(j.eq.m) fa=fa-2.0d0*x(ka)*x(ka-m)- &
       0.5d0*x(ka-m)**2
      if (i.eq.1.and.j.eq.1) fa=fa-par/dble(m+1)
      f=f+fa**2
      g(ka)=g(ka)+1.6d1*x(ka)*fa
      if(i.eq.1) then
      g(ka)=g(ka)-(2.0d0*(x(ka+1)+1.0d0)+3.0d0*x(ka)* &
       (x(ka+1)-1.0d0)*par)*fa
      g(ka+1)=g(ka+1)-(2.0d0*x(ka)+(x(ka+1)-1.0d0)+ &
       1.5d0*x(ka)**2*par)*fa
      endif
      if(i.gt.1.and.i.lt.m) then
      g(ka)=g(ka)-(2.0d0*(x(ka+1)+x(ka-1))+3.0d0*x(ka)* &
       (x(ka+1)-x(ka-1))*par)*fa
      g(ka-1)=g(ka-1)-(2.0d0*x(ka)-(x(ka+1)-x(ka-1))- &
       1.5d0*x(ka)**2*par)*fa
      g(ka+1)=g(ka+1)-(2.0d0*x(ka)+(x(ka+1)-x(ka-1))+ &
       1.5d0*x(ka)**2*par)*fa
      endif
      if(i.eq.m) then
      g(ka)=g(ka)-(2.0d0*x(ka-1)-3.0d0*x(ka)*x(ka-1)*par)*fa
      g(ka-1)=g(ka-1)-(2.0d0*x(ka)+x(ka-1)-1.5d0*x(ka)**2*par)*fa
      endif
      if(j.eq.1) then
      g(ka)=g(ka)-2.0d0*(x(ka+m)+1.0d0)*fa
      g(ka+m)=g(ka+m)-(2.0d0*x(ka)+(x(ka+m)-1.0d0))*fa
      endif
      if(j.gt.1.and.j.lt.m) then
      g(ka)=g(ka)-2.0d0*(x(ka+m)+x(ka-m))*fa
      g(ka-m)=g(ka-m)-(2.0d0*x(ka)-(x(ka+m)-x(ka-m)))*fa
      g(ka+m)=g(ka+m)-(2.0d0*x(ka)+(x(ka+m)-x(ka-m)))*fa
      endif
      if(j.eq.m) then
      g(ka)=g(ka)-2.0d0*x(ka-m)*fa
      g(ka-m)=g(ka-m)-(2.0d0*x(ka)+x(ka-m))*fa
      endif
  971 continue
      f=0.5d0*f
      return
  980 do 981 ka=1,na
      a3=0.0d0
      j=(ka-1)/m+1
      i=ka-(j-1)*m
      a1=par*dble(i)
      a2=par*dble(j)
      fa=4.0d0*x(ka)-2.0d3*a1*a2*(1.0d0-a1)*(1.0d0-a2)*par**2
      if(i.gt.1) then
      fa=fa-x(ka-1)
      a3=a3-x(ka-1)
      endif
      if(i.lt.m) then
      fa=fa-x(ka+1)
      a3=a3+x(ka+1)
      endif
      if(j.gt.1) then
      fa=fa-x(ka-m)
      a3=a3-x(ka-m)
      endif
      if(j.lt.m) then
      fa=fa-x(ka+m)
      a3=a3+x(ka+m)
      endif
      fa=fa+2.0d1*par*a3*x(ka)
      f=f+fa**2
      g(ka)=g(ka)+4.0d0*fa
      if(i.gt.1) then
      g(ka-1)=g(ka-1)-(1.0d0+2.0d1*par*x(ka))*fa
      endif
      if(i.lt.m) then
      g(ka+1)=g(ka+1)-(1.0d0-2.0d1*par*x(ka))*fa
      endif
      if(j.gt.1) then
      g(ka-m)=g(ka-m)-(1.0d0+2.0d1*par*x(ka))*fa
      endif
      if(j.lt.m) then
      g(ka+m)=g(ka+m)-(1.0d0-2.0d1*par*x(ka))*fa
      endif
      g(ka)=g(ka)+2.0d1*par*a3*fa
  981 continue
      f=0.5d0*f
      return
  990 do 991 ka=1,na
      j=(ka-1)/m+1
      i=ka-(j-1)*m
      fa=2.0d1*x(ka)-par*max(0.0d0,x(ka))- &
       sign(par,(dble(i)/dble(m+2)-0.5d0))
      if (j.gt.2) then
        fa=fa+x(ka-m-m)
      endif
      if (j.gt.1) then
        if (i.gt.1) then
          fa=fa+2.0d0*x(ka-m-1)
        endif
        fa=fa-8.0d0*x(ka-m)
        if (i.lt.m) then
          fa=fa+2.0d0*x(ka-m+1)
        endif
      endif
      if (i.gt.1) then
        if (i.gt.2) then
          fa=fa+x(ka-2)
        endif
        fa=fa-8.0d0*x(ka-1)
      endif
      if (i.lt.m) then
        fa=fa-8.0d0*x(ka+1)
        if (i.lt.m-1) then
          fa=fa+x(ka+2)
        endif
      endif
      if (j.lt.m) then
        if (i.gt.1) then
          fa=fa+2.0d0*x(ka+m-1)
        endif
        fa=fa-8.0d0*x(ka+m)
        if (i.lt.m) then
          fa=fa+2.0d0*x(ka+m+1)
        endif
      endif
      if (j.lt.m-1) then
        fa=fa+x(ka+m+m)
      endif
      f=f+fa**2
      g(ka)=g(ka)+(2.0d1-par)*fa
      if (j.gt.2) then
        g(ka-m-m)=g(ka-m-m)+fa
      endif
      if (j.gt.1) then
        if (i.gt.1) then
          g(ka-m-1)=g(ka-m-1)+2.0d0*fa
        endif
        g(ka-m)=g(ka-m)-8.0d0*fa
        if (i.lt.m) then
          g(ka-m+1)=g(ka-m+1)+2.0d0*fa
        endif
      endif
      if (i.gt.1) then
        if (i.gt.2) then
          g(ka-2)=g(ka-2)+fa
        endif
        g(ka-1)=g(ka-1)-8.0d0*fa
      endif
      if (i.lt.m) then
        g(ka+1)=g(ka+1)-8.0d0*fa
        if (i.lt.m-1) then
          g(ka+2)=g(ka+2)+fa
        endif
      endif
      if (j.lt.m) then
        if (i.gt.1) then
          g(ka+m-1)=g(ka+m-1)+2.0d0*fa
        endif
        g(ka+m)=g(ka+m)-8.0d0*fa
        if (i.lt.m) then
          g(ka+m+1)=g(ka+m+1)+2.0d0*fa
        endif
      endif
      if (j.lt.m-1) then
        g(ka+m+m)=g(ka+m+m)+fa
      endif
  991 continue
      f=0.5d0*f
      return
  800 do 802 ka=1,na
      h=0.5d0/dble(m+2)
      j=(ka-1)/m+1
      i=ka-(j-1)*m
      fa=2.0d1*x(ka)
      a1=0.0d0
      a2=0.0d0
      a3=0.0d0
      a4=0.0d0
      if (j.gt.2) then
        fa=fa+x(ka-m-m)
        a4=a4+x(ka-m-m)
      endif
      if (j.gt.1) then
        if (i.gt.1) then
          fa=fa+2.0d0*x(ka-m-1)
          a3=a3+x(ka-m-1)
          a4=a4+x(ka-m-1)
        endif
        fa=fa-8.0d0*x(ka-m)
        a1=a1-x(ka-m)
        a4=a4-4.0d0*x(ka-m)
        if (i.lt.m) then
          fa=fa+2.0d0*x(ka-m+1)
          a3=a3-x(ka-m+1)
          a4=a4+x(ka-m+1)
        endif
      endif
      if (i.gt.1) then
        if (i.gt.2) then
          fa=fa+x(ka-2)
          a3=a3+x(ka-2)
        endif
        fa=fa-8.0d0*x(ka-1)
        a2=a2-x(ka-1)
        a3=a3-4.0d0*x(ka-1)
      endif
      if (i.lt.m) then
        fa=fa-8.0d0*x(ka+1)
        a2=a2+x(ka+1)
        a3=a3+4.0d0*x(ka+1)
        if (i.lt.m-1) then
          fa=fa+x(ka+2)
          a3=a3-x(ka+2)
        endif
      endif
      if (j.lt.m) then
        if (i.gt.1) then
          fa=fa+2.0d0*x(ka+m-1)
          a3=a3+x(ka+m-1)
          a4=a4-x(ka+m-1)
        endif
        fa=fa-8.0d0*x(ka+m)
        a1=a1+x(ka+m)
        a4=a4+4.0d0*x(ka+m)
        if (i.lt.m) then
          fa=fa+2.0d0*x(ka+m+1)
          a3=a3-x(ka+m+1)
          a4=a4-x(ka+m+1)
        endif
      endif
      if (j.lt.m-1) then
        fa=fa+x(ka+m+m)
        a4=a4-x(ka+m+m)
      endif
      if (j.eq.m) then
        if (i.gt.1) then
          fa=fa-h-h
          a3=a3-h
          a4=a4+h
        endif
        fa=fa+8.0d0*h
        a1=a1-h
        a4=a4-4.0d0*h
        if (i.lt.m) then
          fa=fa-2.0d0*h
          a3=a3+h
          a4=a4+h
        endif
        fa=fa+h
        a4=a4-h
      endif
      if (j.eq.m-1) then
        fa=fa-h
        a4=a4+h
      endif
      fa=fa+0.25d0*par*(a1*a3-a2*a4)
      f=f+fa**2
      g(ka)=g(ka)+2.0d1*fa
      a1=0.0d0
      a2=0.0d0
      a3=0.0d0
      a4=0.0d0
      ga1(1)=0.0d0
      ga1(2)=0.0d0
      ga2(1)=0.0d0
      ga2(2)=0.0d0
      do 801 k=1,6
      ga3(k)=0.0d0
      ga4(k)=0.0d0
  801 continue
      if (j.gt.2) then
        g(ka-m-m)=g(ka-m-m)+fa
        ga4(1)=ga4(1)+1.0d0
        a4=a4+x(ka-m-m)
      endif
      if (j.gt.1) then
        if (i.gt.1) then
          g(ka-m-1)=g(ka-m-1)+2.0d0*fa
          ga3(1)=ga3(1)+1.0d0
          ga4(2)=ga4(2)+1.0d0
          a3=a3+x(ka-m-1)
          a4=a4+x(ka-m-1)
        endif
        g(ka-m)=g(ka-m)-8.0d0*fa
        ga1(1)=ga1(1)-1.0d0
        a1=a1-x(ka-m)
        if (i.lt.m) then
          g(ka-m+1)=g(ka-m+1)+2.0d0*fa
          ga3(2)=ga3(2)-1.0d0
          ga4(3)=ga4(3)+1.0d0
          a3=a3-x(ka-m+1)
          a4=a4+x(ka-m+1)
        endif
      endif
      if (i.gt.1) then
        if (i.gt.2) then
          g(ka-2)=g(ka-2)+fa
          ga3(3)=ga3(3)+1.0d0
          a3=a3+x(ka-2)
        endif
        g(ka-1)=g(ka-1)-8.0d0*fa
        ga2(1)=ga2(1)-1.0d0
        a2=a2-x(ka-1)
      endif
      if (i.lt.m) then
        g(ka+1)=g(ka+1)-8.0d0*fa
        ga2(2)=ga2(2)+1.0d0
        a2=a2+x(ka+1)
        if (i.lt.m-1) then
          g(ka+2)=g(ka+2)+fa
          ga3(4)=ga3(4)-1.0d0
          a3=a3-x(ka+2)
        endif
      endif
      if (j.lt.m) then
        if (i.gt.1) then
          g(ka+m-1)=g(ka+m-1)+2.0d0*fa
          ga3(5)=ga3(5)+1.0d0
          ga4(4)=ga4(4)-1.0d0
          a3=a3+x(ka+m-1)
          a4=a4-x(ka+m-1)
        endif
        g(ka+m)=g(ka+m)-8.0d0*fa
        ga1(2)=ga1(2)+1.0d0
        a1=a1+x(ka+m)
        if (i.lt.m) then
          g(ka+m+1)=g(ka+m+1)+2.0d0*fa
          ga3(6)=ga3(6)-1.0d0
          ga4(5)=ga4(5)-1.0d0
          a3=a3-x(ka+m+1)
          a4=a4-x(ka+m+1)
        endif
      endif
      if (j.lt.m-1) then
        g(ka+m+m)=g(ka+m+m)+fa
        ga4(6)=ga4(6)-1.0d0
        a4=a4-x(ka+m+m)
      endif
      if (j.eq.m) then
        if (i.gt.1) then
          a3=a3-h
          a4=a4+h
        endif
        a1=a1-h
        if (i.lt.m) then
          a3=a3+h
          a4=a4+h
        endif
        a4=a4-h
      endif
      if (j.eq.m-1) then
        a4=a4+h
      endif
      if (ka.gt.m+m) &
       g(ka-m-m)=g(ka-m-m)+0.25d0*par*(-a2*ga4(1))*fa
      if (ka.gt.m+1) &
       g(ka-m-1)=g(ka-m-1)+0.25d0*par*(+a1*ga3(1)-a2*ga4(2))*fa
      if (ka.gt.m) &
       g(ka-m)=g(ka-m)+0.25d0*par*(ga1(1)*a3)*fa
      if (ka.gt.m-1) &
       g(ka-m+1)=g(ka-m+1)+0.25d0*par*(+a1*ga3(2)-a2*ga4(3))*fa
      if (ka.gt.2) &
       g(ka-2)=g(ka-2)+0.25d0*par*(+a1*ga3(3))*fa
      if (ka.gt.1) &
       g(ka-1)=g(ka-1)+0.25d0*par*(-ga2(1)*a4)*fa
      if (ka.le.n-1) &
       g(ka+1)=g(ka+1)+0.25d0*par*(-ga2(2)*a4)*fa
      if (ka.le.n-2) &
       g(ka+2)=g(ka+2)+0.25d0*par*(+a1*ga3(4))*fa
      if (ka.le.n-m+1) &
       g(ka+m-1)=g(ka+m-1)+0.25d0*par*(+a1*ga3(5)-a2*ga4(4))*fa
      if (ka.le.n-m) &
       g(ka+m)=g(ka+m)+0.25d0*par*(ga1(2)*a3)*fa
      if (ka.le.n-m-1) &
       g(ka+m+1)=g(ka+m+1)+0.25d0*par*(+a1*ga3(6)-a2*ga4(5))*fa
      if (ka.le.n-m-m) &
       g(ka+m+m)=g(ka+m+m)+0.25d0*par*(-a2*ga4(6))*fa
  802 continue
      f=0.5d0*f
      return
  240 do 243 ka=1,na
      w=0.0d0
      do 241 i=1,n-1
      w=w+(dble(ka)/dble(ka+i))*x(i)
  241 continue
      fa=x(ka)-(1.0d0+(0.4d0/dble(n))*x(ka)*(0.5d0+w+ &
       0.5d0*(dble(ka)/dble(ka+n))*x(n)))
      f=f+fa**2
      w=w+0.5d0+0.5d0*dble(ka)/dble(ka+n)*x(n)
      do 242 i=1,n-1
      g(i)=g(i)-0.4d0/dble(n)*x(ka)*dble(ka)/dble(ka+i)*fa
  242 continue
      g(n)=g(n)-0.2d0/dble(n)*x(ka)*dble(ka)/dble(ka+n)*fa
      g(ka)=g(ka)+(1.0d0-0.4d0*w/dble(n))*fa
  243 continue
      f=0.5d0*f
      return
  410 do 411 ka=1,na
      if(ka.eq.1) then
        fa=1.0d0-x(1)
      g(1)=g(1)-fa
      else
        fa=10.0d0*dble(ka-1)*(x(ka)-x(ka-1))**2
        g(ka)=g(ka)+20.0d0*dble(ka-1)*(x(ka)-x(ka-1))*fa
        g(ka-1)=g(ka-1)-20.0d0*dble(ka-1)*(x(ka)-x(ka-1))*fa
      endif
      f=f+fa**2
  411 continue
      f=0.5d0*f
      return
  420 do 421 ka=1,na
      if(ka.eq.n) then
        fa=x(ka)-0.1d0*x(1)**2
        g(1)=g(1)-0.20d0*x(1)*fa
        g(n)=g(n)+fa
      else
        fa=x(ka)-0.1d0*x(ka+1)**2
        g(ka)=g(ka)+fa
        g(ka+1)=g(ka+1)-0.2d0*x(ka+1)*fa
      endif
      f=f+fa**2
  421 continue
      f=0.5d0*f
      return
  650 do 653 ka=1,na
      s=0.0d0
      do 651 j=1,n
        s=s+x(j)**3
  651 continue
      fa=x(ka)-1.0d0/dble(2*n)*(s+dble(ka))
      f=f+fa**2
      do 652 j=1,n
        if(j.eq.ka) then
           g(j)=g(j)+(1.0d0-3.0d0*x(j)**2/(2.0d0*dble(n)))*fa
        else
           g(j)=g(j)-3.0d0*x(j)**2/(2.0d0*dble(n))*fa
        endif
  652 continue
  653 continue
      f=0.5d0*f
      return
  660 do 661 ka=1,na
      s=(1.0d0/dble(n+1))**2*exp(x(ka))
      if(n.eq.1) then
        fa=-2.0d0*x(ka)-s
        g(ka)=g(ka)-(2.0d0+s)*fa
      else if(ka.eq.1) then
        fa=-2.0d0*x(ka)+x(ka+1)-s
        g(ka)=g(ka)-(2.0d0+s)*fa
        g(ka+1)=g(ka+1)+fa
      else if(ka.eq.n) then
        fa=x(ka-1)-2.0d0*x(ka)-s
        g(ka)=g(ka)-(2.0d0+s)*fa
        g(ka-1)=g(ka-1)+fa
      else
        fa=x(ka-1)-2.0d0*x(ka)+x(ka+1)-s
        g(ka)=g(ka)-(2.0d0+s)*fa
        g(ka-1)=g(ka-1)+fa
        g(ka+1)=g(ka+1)+fa
      endif
      f=f+fa**2
  661 continue
      f=0.5d0*f
      return
  670 do 671 ka=1,na
      s=0.1d0
      if(n.eq.1) then
        fa=(3.0d0-s*x(ka))*x(ka)+1.0d0
        g(ka)=g(ka)+(3.0d0-2.0d0*s*x(ka))*fa
      else if(ka.eq.1) then
        fa=(3.0d0-s*x(ka))*x(ka)+1.0d0-2.0d0*x(ka+1)
        g(ka)=g(ka)+(3.0d0-2.0d0*s*x(ka))*fa
        g(ka+1)=g(ka+1)-2.0d0*fa
      else if(ka.eq.n) then
        fa=(3.0d0-s*x(ka))*x(ka)+1.0d0-x(ka-1)
        g(ka)=g(ka)+(3.0d0-2.0d0*s*x(ka))*fa
        g(ka-1)=g(ka-1)-fa
      else
        fa=(3.0d0-s*x(ka))*x(ka)+1.0d0-x(ka-1)-2.0d0*x(ka+1)
        g(ka)=g(ka)+(3.0d0-2.0d0*s*x(ka))*fa
        g(ka-1)=g(ka-1)-fa
        g(ka+1)=g(ka+1)-2.0d0*fa
      endif
      f=f+fa**2
  671 continue
      f=0.5d0*f
      return
  680 do 683 ka=1,na
      s1=1.0d0
      s2=1.0d0
      s3=1.0d0
      j1=3
      j2=3
      if(ka-j1.gt.1) then
        i1=ka-j1
      else
        i1=1
      endif
      if(ka+j2.lt.n) then
        i2=ka+j2
      else
        i2=n
      endif
      s=0.0d0
      do 681 j=i1,i2
        if(j.ne.ka) s=s+x(j)+x(j)**2
  681 continue
      fa=(s1+s2*x(ka)**2)*x(ka)+1.d0-s3*s
      g(ka)=g(ka)+(s1+3.0d0*s2*x(ka)**2)*fa
      do 682 j=i1,i2
        g(j)=g(j)-s3*(1.0d0+2.0d0*x(j))*fa
  682 continue
      f=f+fa**2
  683 continue
      f=0.5d0*f
      return
  690 do 691 ka=1,na
      if(ka.eq.1) then
        fa=x(1)**2-1.0d0
        g(1)=g(1)+2.0d0*x(1)*fa
      else
        fa=x(ka-1)**2+log(x(ka))-1.0d0
        g(ka-1)=g(ka-1)+2.0d0*x(ka-1)*fa
        g(ka)=g(ka)+1.0d0/x(ka)*fa
      endif
      f=f+fa**2
  691 continue
      f=0.5d0*f
      return
  340 do 341 ka=1,na
      if(ka.eq.1) then
        fa=x(1)
        g(1)=g(1)+fa
      else
        fa=cos(x(ka-1))+x(ka)-1.0d0
        g(ka)=g(ka)+fa
        g(ka-1)=g(ka-1)-sin(x(ka-1))*fa
      endif
      f=f+fa**2
  341 continue
      f=0.5d0*f
      return
  360 do 361 ka=1,na
      s=(1.0d0/dble(n+1))**2
      if(n.eq.1) then
        fa=2.0d0*x(ka)-1.0d0+s*(x(ka)+sin(x(ka)))
        g(ka)=g(ka)+(2.0d0+s*(1.0d0+cos(x(ka))))*fa
      else if(ka.eq.1) then
        fa=2.0d0*x(ka)-x(ka+1)+s*(x(ka)+sin(x(ka)))
        g(ka)=g(ka)+(2.0d0+s*(1.0d0+cos(x(ka))))*fa
        g(ka+1)=g(ka+1)-fa
      else if(ka.eq.n) then
        fa=-x(ka-1)+2.0d0*x(ka)-1.0d0+s*(x(ka)+sin(x(ka)))
        g(ka)=g(ka)+(2.0d0+s*(1.0d0+cos(x(ka))))*fa
        g(ka-1)=g(ka-1)-fa
      else
        fa=-x(ka-1)+2.0d0*x(ka)-x(ka+1)+s*(x(ka)+sin(x(ka)))
        g(ka)=g(ka)+(2.0d0+s*(1.0d0+cos(x(ka))))*fa
        g(ka-1)=g(ka-1)-fa
        g(ka+1)=g(ka+1)-fa
      endif
      f=f+fa**2
  361 continue
      f=0.5d0*f
      return
  380 do 383 ka=1,na
      if(ka-5.gt.1) then
        i1=ka-5
      else
        i1=1
      endif
      if(ka+1.lt.n) then
        i2=ka+1
      else
        i2=n
      endif
      s=0.0d0
      do 381 j=i1,i2
        if(j.ne.ka) s=s+x(j)*(1.0d0+x(j))
  381 continue
      fa=x(ka)*(2.0d0+5.0d0*x(ka)**2)+1.0d0-s
      f=f+fa**2
      g(ka)=g(ka)+(2.0d0+15.0d0*x(ka)**2)*fa
      do 382 j=i1,i2
      if(j.ne.ka) g(j)=g(j)-(1.0d0+2.0d0*x(j))*fa
  382 continue
  383 continue
      f=0.5d0*f
      return
  430 do 433 ka=1,na
      alf=5
      bet=14
      gam=3
      fa=dble(bet*n)*x(ka)+(dble(ka)-dble(n)/2.0d0)**gam
      do 431 j=1,n
        if(j.ne.ka) then
        t=sqrt(x(j)**2+dble(ka)/dble(j))
        s1=log(t)
        fa=fa+t*(sin(s1)**alf+cos(s1)**alf)
      endif
  431 continue
      f=f+fa**2
      do 432 j=1,n
        if(j.ne.ka) then
        t=sqrt(x(j)**2+dble(ka)/dble(j))
        s1=log(t)
        g(j)=g(j)+(x(j)*(sin(s1)**alf+cos(s1)**alf+ &
           alf*sin(s1)**(alf-1)*cos(s1)- &
           alf*sin(s1)*cos(s1)**(alf-1))/t)*fa
        else
          g(j)=g(j)+dble(bet*n)*fa
        endif
  432 continue
  433 continue
      f=0.5d0*f
      return
  440 do 445 ka=1,na
      c=0.5d0
      h=1.0d0/dble(n)
      fa=(1.0d0-c*h/4.0d0)
      do 441 j=1,n
        s=c*h*dble(ka)/dble(2*(ka+j))
        if(j.eq.n) s=s/2.0d0
        fa=fa-s*x(j)
  441 continue
      fa=-1.0d0+x(ka)*fa
      f=f+fa**2
      do 442 j=1,n
        sx(j)=c*h*dble(ka)/dble(2*(ka+j))
  442 continue
      sx(n)=0.5d0*sx(n)
      do 444 j=1,n
        if(ka.ne.j) then
           g(j)=g(j)-sx(j)*x(ka)*fa
        else
           t=1.0d0-c*h/4.0d0
           do 443 l=1,n
             if(l.eq.ka) then
               t=t-2.0d0*sx(ka)*x(ka)
             else
               t=t-sx(l)*x(l)
             endif
  443      continue
           g(j)=g(j)+t*fa
        endif
  444 continue
  445 continue
      f=0.5d0*f
      return
  270 do 271 ka=1,na
      s=0.5d0
      h=1.0d0/dble(n+1)
      t=2.0d0*h**2
      if(ka.eq.1) then
      fa=2.0d0*x(ka)-x(ka+1)-t*x(ka)**2+h*x(ka+1)
      g(ka)=g(ka)+2.0d0*(1.0d0-h**2*x(ka)/s)*fa
      g(ka+1)=g(ka+1)-(1.0d0-h)*fa
      else if(1.lt.ka.and.ka.lt.n) then
      fa=-x(ka-1)+2.0d0*x(ka)-x(ka+1)-t*x(ka)**2+h*(x(ka+1)-x(ka-1))
      g(ka)=g(ka)+2.0d0*(1.0d0-t*x(ka))*fa
      g(ka-1)=g(ka-1)-(1.0d0+h)*fa
      g(ka+1)=g(ka+1)-(1.0d0-h)*fa
      else if(ka.eq.n) then
      fa=-x(ka-1)+2.0d0*x(ka)-0.5d0-t*x(ka)**2+h*(0.5d0-x(ka-1))
      g(ka)=g(ka)+2.0d0*(1.0d0-h**2*x(ka)/s)*fa
      g(ka-1)=g(ka-1)-(1.0d0+h)*fa
      endif
      f=f+fa**2
  271 continue
      f=0.5d0*f
      return
  280 do 288 ka=1,na
      s=0.5d0
      h=1.0d0/dble(n+1)
      t=h**2/s
      t1=2.0d0*h
      al=0.0d0
      be=0.5d0
      s1=0.0d0
      do 281 j=1,ka
        if(j.eq.1) then
           s1=s1+dble(j)*(x(j)**2+(x(j+1)-al)/t1)
        endif
        if(1.lt.j.and.j.lt.n) then
           s1=s1+dble(j)*(x(j)**2+(x(j+1)-x(j-1))/t1)
        endif
        if(j.eq.n) then
           s1=s1+dble(j)*(x(j)**2+(be-x(j-1))/t1)
        endif
  281 continue
      s1=(1.0d0-dble(ka)*h)*s1
      if(ka.eq.n) go to 283
      s2=0.0d0
      do 282 j=ka+1,n
        if(j.lt.n) then
           s2=s2+(1.0d0-dble(j)*h)*(x(j)**2+(x(j+1)-x(j-1))/t1)
        else
           s2=s2+(1.0d0-dble(j)*h)*(x(j)**2+(be-x(j-1))/t1)
        endif
  282 continue
      s1=s1+dble(ka)*s2
  283 fa=x(ka)-0.5d0*dble(ka)*h-t*s1
      f=f+fa**2
      s1=h**2/s
      s2=1.0d0-dble(ka)*h
      do 284 j=1,ka
        sx(j)=dble(j)*s2
  284 continue
      if(ka.eq.n) go to 286
      do 285 j=ka+1,n
        sx(j)=dble(ka)*(1.0d0-dble(j)*h)
  285 continue
  286 g(1)=g(1)-s1*(sx(1)*2.0d0*x(1)-sx(2)/t1)*fa
      g(n)=g(n)-s1*(sx(n-1)/t1+sx(n)*2.0d0*x(n))*fa
      do 287 j=2,n-1
        g(j)=g(j)-s1*((sx(j-1)-sx(j+1))/t1+sx(j)*2.0d0*x(j))*fa
  287 continue
      g(ka)=g(ka)+fa
      return
  288 continue
      f=0.5d0*f
      return
  290 do 291 ka=1,na
      a=-9.0d-3
      b=1.0d-3
      al=0.0d0
      be=25.0d0
      ga=20.0d0
      ca=0.3d0
      cb=0.3d0
      h=(b-a)/dble(n+1)
      t=a+dble(ka)*h
      h=h**2
      s=dble(ka)/dble(n+1)
      u=al*(1.0d0-s)+be*s+x(ka)
      ff=cb*exp(ga*(u-be))-ca*exp(ga*(al-u))
      fg=ff*ga
      if(t.le.0) then
        ff=ff+ca
      else
        ff=ff-cb
      endif
      if(n.eq.1) then
        fa=-al+2.0d0*x(ka)-be+h*ff
        g(ka)=g(ka)+(2.0d0+h*fg)*fa
      elseif(ka.eq.1) then
        fa=-al+2.0d0*x(ka)-x(ka+1)+h*ff
        g(ka)=g(ka)+(2.0d0+h*fg)*fa
        g(ka+1)=g(ka+1)-fa
      elseif(ka.lt.n) then
        fa=-x(ka-1)+2.0d0*x(ka)-x(ka+1)+h*ff
        g(ka)=g(ka)+(2.0d0+h*fg)*fa
        g(ka-1)=g(ka-1)-fa
        g(ka+1)=g(ka+1)-fa
      else
        fa=-x(ka-1)+2.0d0*x(ka)-be+h*ff
        g(ka)=g(ka)+(2.0d0+h*fg)*fa
        g(ka-1)=g(ka-1)-fa
      endif
      f=f+fa**2
  291 continue
      f=0.5d0*f
      return
  300 do 301 ka=1,na
      al1=0.0d0
      al2=0.0d0
      be1=0.0d0
      be2=0.0d0
      n1=n/2
      h=1.0d0/dble(n1+1)
      t=dble(ka)*h
      if(ka.eq.1) then
        s1=2.0d0*x(ka)-x(ka+1)
        b=al1
      else if(ka.eq.n1+1) then
        s1=2.0d0*x(ka)-x(ka+1)
        b=al2
      else if(ka.eq.n1) then
        s1=-x(ka-1)+2.0d0*x(ka)
        b=be1
      else if(ka.eq.n) then
        s1=-x(ka-1)+2.0d0*x(ka)
        b=be2
      else
        s1=-x(ka-1)+2.0d0*x(ka)-x(ka+1)
        b=0.0d0
      endif
      if(ka.le.n1) then
        s2=x(ka)**2+x(ka)+0.1d0*x(n1+ka)**2-1.2d0
      else
        s2=0.2d0*x(ka-n1)**2+x(ka)**2+2.0d0*x(ka)-0.6d0
      endif
      fa=s1+h**2*s2-b
      f=f+fa**2
      h=1.0d0/dble(n1+1)**2
      if(ka.eq.1) then
        g(ka)=g(ka)+(2.0d0+h*(2.0d0*x(ka)+1.0d0))*fa
        g(ka+1)=g(ka+1)-fa
        g(n1+ka)=g(n1+ka)+h*0.2d0*x(n1+ka)*fa
      else if(ka.eq.n1+1) then
        g(1)=g(1)+h*0.4d0*x(1)*fa
        g(ka)=g(ka)+(2.0d0+h*(2.0d0*x(ka)+2.0d0))*fa
        g(ka+1)=g(ka+1)-fa
      else if(ka.eq.n1) then
        g(ka-1)=g(ka-1)-fa
        g(ka)=g(ka)+(2.0d0+h*(2.0d0*x(ka)+1.0d0))*fa
        g(n1+ka)=g(n1+ka)+h*0.2d0*x(n1+ka)*fa
      else if(ka.eq.n) then
        g(n1)=g(n1)+h*0.4d0*x(n1)*fa
        g(ka-1)=g(ka-1)-fa
        g(ka)=g(ka)+(2.0d0+h*(2.0d0*x(ka)+2.0d0))*fa
      else if(ka.lt.n1) then
        g(ka-1)=g(ka-1)-fa
        g(ka)=g(ka)+(2.0d0+h*(2.0d0*x(ka)+1.0d0))*fa
        g(ka+1)=g(ka+1)-fa
        g(n1+ka)=g(n1+ka)+h*0.2d0*x(n1+ka)*fa
      else
         g(ka-n1)=g(ka-n1)+h*0.4d0*x(ka-n1)*fa
         g(ka-1)=g(ka-1)-fa
         g(ka)=g(ka)+(2.0d0+h*(2.0d0*x(ka)+2.0d0))*fa
         g(ka+1)=g(ka+1)-fa
      endif
  301 continue
      f=0.5d0*f
      return
  710 do 711 ka=1,na
      nd=int(sqrt(dble(n)))
      l=mod(ka,nd)
      if(l.eq.0) then
         k=ka/nd
         l=nd
      else
         k=int(ka/nd)+1
      endif
      la=1.0d0
      h=1.0d0/dble(nd+1)
      h2=la*h*h
      if(l.eq.1.and.k.eq.1) then
         fa=4.0d0*x(1)-x(2)-x(nd+1)+h2*exp(x(1))
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         fa=4.0d0*x(l)-x(l-1)-x(l+1)-x(l+nd)+h2*exp(x(l))
      endif
      if(l.eq.nd.and.k.eq.1) then
         fa=4.0d0*x(nd)-x(nd-1)-x(nd+nd)+h2*exp(x(nd))
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         fa=4.0d0*x(ka)-x(ka+1)-x(ka-nd)-x(ka+nd)+h2*exp(x(ka))
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         fa=4.0d0*x(ka)-x(ka-nd)-x(ka-1)-x(ka+nd)+h2*exp(x(ka))
      endif
      if(l.eq.1.and.k.eq.nd) then
         fa=4.0d0*x(ka)-x(ka+1)-x(ka-nd)+h2*exp(x(ka))
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         fa=4.0d0*x(ka)-x(ka-1)-x(ka+1)-x(ka-nd)+h2*exp(x(ka))
      endif
      if(l.eq.nd.and.k.eq.nd) then
        fa=4.0d0*x(ka)-x(ka-1)-x(ka-nd)+h2*exp(x(ka))
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
      fa=4.0d0*x(ka)-x(ka-1)-x(ka+1)-x(ka-nd)-x(ka+nd)+h2*exp(x(ka))
      endif
      f=f+fa**2
      if(l.eq.1.and.k.eq.1) then
         g(1)=g(1)+(4.0d0+h2*exp(x(1)))*fa
         g(2)=g(2)-fa
         g(nd+1)=g(nd+1)-fa
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         g(l)=g(l)+(4.0d0+h2*exp(x(l)))*fa
         g(l-1)=g(l-1)-fa
         g(l+1)=g(l+1)-fa
         g(l+nd)=g(l+nd)-fa
      endif
      if(l.eq.nd.and.k.eq.1) then
         g(nd)=g(nd)+(4.0d0+h2*exp(x(nd)))*fa
         g(nd-1)=g(nd-1)-fa
         g(nd+nd)=g(nd+nd)-fa
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         g(ka)=g(ka)+(4.0d0+h2*exp(x(ka)))*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka+1)=g(ka+1)-fa
         g(ka+nd)=g(ka+nd)-fa
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         g(ka)=g(ka)+(4.0d0+h2*exp(x(ka)))*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
         g(ka+nd)=g(ka+nd)-fa
      endif
      if(l.eq.1.and.k.eq.nd) then
         g(ka)=g(ka)+(4.0d0+h2*exp(x(ka)))*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka+1)=g(ka+1)-fa
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         g(ka)=g(ka)+(4.0d0+h2*exp(x(ka)))*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
         g(ka+1)=g(ka+1)-fa
      endif
      if(l.eq.nd.and.k.eq.nd) then
         g(ka)=g(ka)+(4.0d0+h2*exp(x(ka)))*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
         g(ka)=g(ka)+(4.0d0+h2*exp(x(ka)))*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
         g(ka+1)=g(ka+1)-fa
         g(ka+nd)=g(ka+nd)-fa
      endif
  711 continue
      f=0.5d0*f
      return
  820 do 821 ka=1,na
      nd=int(sqrt(dble(n)))
      l=mod(ka,nd)
      if(l.eq.0) then
         k=ka/nd
         l=nd
      else
         k=int(ka/nd)+1
      endif
      h=1.0d0/dble(nd+1)
      h2=h*h
      if(l.eq.1.and.k.eq.1) then
         fa=4.0d0*x(1)-x(2)-x(nd+1)+h2*x(1)**2-24.0d0/(h+1.0d0)**2
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         fa=4.0d0*x(l)-x(l-1)-x(l+1)-x(l+nd)+h2*x(l)**2 &
      -12.0d0/(dble(l)*h+1.0d0)**2
      endif
      if(l.eq.nd.and.k.eq.1) then
         fa=4.0d0*x(nd)-x(nd-1)-x(nd+nd)+h2*x(nd)**2 &
      -12.0d0/(dble(nd)*h+1.0d0)**2-12.0d0/(h+2.0d0)**2
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         fa=4.0d0*x(ka)-x(ka+1)-x(ka-nd)-x(ka+nd)+h2*x(ka)**2 &
      -12.0d0/(dble(k)*h+1.0d0)**2
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         fa=4.0d0*x(ka)-x(ka-nd)-x(ka-1)-x(ka+nd)+h2*x(ka)**2 &
      -12.0d0/(dble(k)*h+2.0d0)**2
      endif
      if(l.eq.1.and.k.eq.nd) then
         fa=4.0d0*x(ka)-x(ka+1)-x(ka-nd)+h2*x(ka)**2 &
      -12.0d0/(dble(nd)*h+1.0d0)**2-12.0d0/(h+2.0d0)**2
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         fa=4.0d0*x(ka)-x(ka-1)-x(ka+1)-x(ka-nd)+h2*x(ka)**2 &
      -12.0d0/(dble(l)*h+2.0d0)**2
      endif
      if(l.eq.nd.and.k.eq.nd) then
        fa=4.0d0*x(ka)-x(ka-1)-x(ka-nd)+h2*x(ka)**2 &
      -24.0d0/(dble(nd)*h+2.0d0)**2
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
         fa=4.0d0*x(ka)-x(ka-1)-x(ka+1)-x(ka-nd)-x(ka+nd)+h2*x(ka)**2
      endif
      f=f+fa**2
      if(l.eq.1.and.k.eq.1) then
         g(1)=g(1)+(4.0d0+h2*x(1)*2.0d0)*fa
         g(2)=g(2)-fa
         g(nd+1)=g(nd+1)-fa
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.1) then
         g(l)=g(l)+(4.0d0+h2*x(l)*2.0d0)*fa
         g(l-1)=g(l-1)-fa
         g(l+1)=g(l+1)-fa
         g(l+nd)=g(l+nd)-fa
      endif
      if(l.eq.nd.and.k.eq.1) then
         g(nd)=g(nd)+(4.0d0+h2*x(nd)*2.0d0)*fa
         g(nd-1)=g(nd-1)-fa
         g(nd+nd)=g(nd+nd)-fa
      endif
      if(l.eq.1.and.1.lt.k.and.k.lt.nd) then
         g(ka)=g(ka)+(4.0d0+h2*x(ka)*2.0d0)*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka+1)=g(ka+1)-fa
         g(ka+nd)=g(ka+nd)-fa
      endif
      if(l.eq.nd.and.1.lt.k.and.k.lt.nd) then
         g(ka)=g(ka)+(4.0d0+h2*x(ka)*2.0d0)*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
         g(ka+nd)=g(ka+nd)-fa
      endif
      if(l.eq.1.and.k.eq.nd) then
         g(ka)=g(ka)+(4.0d0+h2*x(ka)*2.0d0)*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka+1)=g(ka+1)-fa
      endif
      if(1.lt.l.and.l.lt.nd.and.k.eq.nd) then
         g(ka)=g(ka)+(4.0d0+h2*x(ka)*2.0d0)*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
         g(ka+1)=g(ka+1)-fa
      endif
      if(l.eq.nd.and.k.eq.nd) then
         g(ka)=g(ka)+(4.0d0+h2*x(ka)*2.0d0)*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
      endif
      if(1.lt.l.and.l.lt.nd.and.1.lt.k.and.k.lt.nd) then
         g(ka)=g(ka)+(4.0d0+h2*x(ka)*2.0d0)*fa
         g(ka-nd)=g(ka-nd)-fa
         g(ka-1)=g(ka-1)-fa
         g(ka+1)=g(ka+1)-fa
         g(ka+nd)=g(ka+nd)-fa
      endif
  821 continue
      f=0.5d0*f
      return
      end subroutine tfbu28

! subroutine tild22             all systems                 99/12/01
! portability : all systems
! 94/12/01 lu : original version
!
! purpose :
!  initiation of variables for nonlinear minimax approximation.
!  linearly constrained dense version.
!
! parameters :
!  io  n  number of variables.
!  io  na  number of partial functions.
!  io  nb  number of box constraints.
!  io  nc  number of general linear constraints.
!  ro  x(n)  vector of variables.
!  io  ix(nf)  vector containing types of bounds.
!  ro  xl(nf)  vector containing lower bounds for variables.
!  ro  xu(nf)  vector containing upper bounds for variables.
!  io  ic(nc)  vector containing types of constraints.
!  ro  cl(nc)  vector containing lower bounds for constraint functions.
!  ro  cu(nc)  vector containing upper bounds for constraint functions.
!  ro  cg(nf*nc) matrix whose columns are normals of the linear
!         constraints.
!  ro  fmin  lower bound for value of the objective function.
!  ro  xmax  maximum stepsize.
!  io  next  number of the test problem.
!  io  iext  type of objective function. iext<0-maximum of values.
!         iext=0-maximum of absolute values.
!  io  ierr  error indicator.
!
      subroutine tild22(n,na,nb,nc,x,ix,xl,xu,ic,cl,cu,cg,fmin,xmax,&
                        next,iext,ierr)
      double precision pi
      parameter (pi=3.14159265358979323846d0)
      double precision fmin,xmax
      integer ierr,iext,n,na,nb,nc,next
      double precision cg(n*nc),cl(nc),cu(nc),x(n),xl(n),xu(n)
      integer ic(nc),ix(n)
      double precision y(163)
      integer i,j,k,l
      common /empr22/y
      fmin = -1.0d60
      xmax = 1.0d3
      iext = -1
      ierr = 0
      nb = 0
      go to (10,20,30,40,50,100,160) next
   10 if (n.ge.2 .and. na.ge.3) then
          n = 2
          na = 3
          nc = 1
          x(1) = 1.0d0
          x(2) = 2.0d0
          ic(1) = 1
          cl(1) = 0.5d0
          cg(1) = 1.0d0
          cg(2) = 1.0d0
      else
          ierr = 1
      end if
      return
   20 if (n.ge.2 .and. na.ge.3) then
          n = 2
          na = 3
          nc = 1
          x(1) = -2.0d0
          x(2) = -1.0d0
          ic(1) = 2
          cu(1) = -2.5d0
          cg(1) = 3.0d0
          cg(2) = 1.0d0
      else
          ierr = 1
      end if
      return
   30 if (n.ge.2 .and. na.ge.3) then
          n = 2
          na = 3
          nb = n
          nc = 1
          x(1) = -1.0d0
          x(2) = 1.0d-2
          ix(1) = 0
          ix(2) = 1
          xl(2) = 1.0d-2
          ic(1) = 1
          cl(1) = -5.0d-1
          cg(1) = 5.0d-2
          cg(2) = -1.0d0
      else
          ierr = 1
      end if
      return
   40 if (n.ge.2 .and. na.ge.3) then
          n = 2
          na = 3
          nb = n
          nc = 1
          x(1) = -1.0d0
          x(2) = 3.0d0
          ix(1) = 0
          ix(2) = 1
          xl(2) = 1.0d-2
          ic(1) = 1
          cl(1) = 1.0d0
          cg(1) = -9.0d-1
          cg(2) = 1.0d0
      else
          ierr = 1
      end if
      return
   50 if (n.ge.6 .and. na.ge.3) then
          n = 6
          na = 3
          x(1) = -1.0d0
          x(2) = 0.0d0
          x(3) = 0.0d0
          x(4) = -1.0d0
          x(5) = 1.0d0
          x(6) = 1.0d0
          nc = 5*na
          do 60 i = 1,nc
              cu(i) = 1.0d0
              ic(i) = 2
   60     continue
          do 70 i = 1,n*nc
              cg(i) = 0.0d0
   70     continue
          k = 1
          do 90 i = 1,na
              l = 2* (i-1)
              do 80 j = 1,5
                  cg(k+l) = sin(2.0d0*pi*dble(j-1)/5.0d0)
                  cg(k+l+1) = cos(2.0d0*pi*dble(j-1)/5.0d0)
                  k = k + n
   80         continue
   90     continue
      else
          ierr = 1
      end if
      return
  100 if (n.ge.7 .and. na.ge.163) then
          n = 7
          na = 163
          nb = n
          do 110 i = 1,n
              x(i) = dble(i)*0.5d0
              ix(i) = 0
  110     continue
          xl(1) = 0.4d0
          ix(1) = 1
          ix(7) = 5
          do 120 i = 1,na
              y(i) = 2.0d0*pi*sin(pi* (8.5d0+dble(i)*0.5d0)/180.0d0)
  120     continue
          nc = 7
          do 130 i = 1,6
              cl(i) = 0.4d0
              ic(i) = 1
  130     continue
          cl(7) = 1.0d0
          cu(7) = 1.0d0
          ic(7) = 5
          do 140 i = 1,n*nc
              cg(i) = 0.0d0
  140     continue
          k = 0
          do 150 i = 1,6
              cg(k+i) = -1.0d0
              cg(k+i+1) = 1.0d0
              k = k + n
  150     continue
          cg(46) = -1.0d0
          cg(48) = 1.0d0
          iext = 0
          fmin = 0.0d0
      else
          ierr = 1
      end if
      return
  160 if (n.ge.8 .and. na.ge.8) then
          n = 8
          na = 8
          nb = n
          do 170 i = 1,n
              x(i) = 0.125d0
              xl(i) = 1.0d-8
              ix(i) = 1
  170     continue
          do 180 i = 1,40
              y(i) = 1.0d0
              y(i+40) = 0.1d0
  180     continue
          y(9) = 2.0d0
          y(10) = 0.8d0
          y(12) = 0.5d0
          y(18) = 1.2d0
          y(19) = 0.8d0
          y(20) = 1.2d0
          y(21) = 1.6d0
          y(22) = 2.0d0
          y(23) = 0.6d0
          y(24) = 0.1d0
          y(25) = 2.0d0
          y(26) = 0.1d0
          y(27) = 0.6d0
          y(28) = 2.0d0
          y(32) = 2.0d0
          y(33) = 1.2d0
          y(34) = 1.2d0
          y(35) = 0.8d0
          y(37) = 1.2d0
          y(38) = 0.1d0
          y(39) = 3.0d0
          y(40) = 4.0d0
          y(41) = 3.0d0
          y(42) = 1.0d0
          y(45) = 5.0d0
          y(48) = 6.0d0
          y(50) = 1.0d1
          y(53) = 5.0d0
          y(58) = 9.0d0
          y(59) = 1.0d1
          y(61) = 4.0d0
          y(63) = 7.0d0
          y(68) = 1.0d1
          y(70) = 3.0d0
          y(80) = 1.1d1
          y(81) = 0.5d0
          y(82) = 1.2d0
          y(83) = 0.8d0
          y(84) = 2.0d0
          y(85) = 1.5d0
          nc = 1
          cl(1) = 1.0d0
          cu(1) = 1.0d0
          ic(1) = 5
          do 190 i = 1,n
              cg(i) = 1.0d0
  190     continue
      else
          ierr = 1
      end if
      return
      end subroutine tild22

! subroutine tafu22             all systems                 99/12/01
! portability : all systems
! 94/12/01 lu : original version
!
! purpose :
!  values of partial functions in the minimax criterion.
!
! parameters :
!  ii  n  number of variables.
!  ii  ka  index of the partial function.
!  ri  x(n)  vector of variables.
!  ro  fa  value of the partial function at the
!          selected point.
!  ii  next  number of the test problem.
!
      subroutine tafu22(n,ka,x,fa,next)
      double precision fa
      integer ka,n,next
      double precision x(n)
      double precision y(163)
      double precision a,p
      integer i,j,k
      common /empr22/y
      go to (10,10,50,50,90,130,150) next
   10 go to (20,30,40) ka
   20 fa = x(1)**2 + x(2)**2 + x(1)*x(2) - 1.0d0
      return
   30 fa = sin(x(1))
      return
   40 fa = -cos(x(2))
      return
   50 go to (60,70,80) ka
   60 fa = -exp(x(1)-x(2))
      return
   70 fa = sinh(x(1)-1.0d0) - 1.0d0
      return
   80 fa = -log(x(2)) - 1.0d0
      return
   90 go to (100,110,120) ka
  100 fa = -sqrt((x(1)-x(3))**2+ (x(2)-x(4))**2)
      return
  110 fa = -sqrt((x(3)-x(5))**2+ (x(4)-x(6))**2)
      return
  120 fa = -sqrt((x(5)-x(1))**2+ (x(6)-x(2))**2)
      return
  130 a = 0.0d0
      do 140 i = 1,n
          a = a + cos(y(ka)*x(i))
  140 continue
      fa = (1.0d0+2.0d0*a)/1.5d1
      return
  150 fa = 0.0d0
      k = 0
      do 170 i = 1,5
          a = 0.0d0
          p = 0.0d0
          do 160 j = 1,n
              a = a + y(k+j)*x(j)** (1.0d0-y(i+80))
              p = p + y(k+j+40)*x(j)
  160     continue
          fa = fa + y(k+ka)*p/ (x(ka)**y(i+80)*a) - y(k+ka+40)
          k = k + n
  170 continue
      return
      end subroutine tafu22

! subroutine tagu22             all systems                 99/12/01
! portability : all systems
! 94/12/01 lu : original version
!
! purpose :
!  gradients of partial functions in the minimax criterion.
!
! parameters :
!  ii  n  number of variables.
!  ii  ka  index of the partial function.
!  ri  x(n)  vector of variables.
!  ro  ga(n)  gradient of the partial function at the
!          selected point.
!  ii  next  number of the test problem.
!
      subroutine tagu22(n,ka,x,ga,next)
      integer ka,n,next
      double precision ga(n),x(n)
      double precision y(163)
      double precision a,b,c,p
      integer i,j,k
      common /empr22/y
      go to (10,10,50,50,90,130,150) next
   10 go to (20,30,40) ka
   20 ga(1) = 2.0d0*x(1) + x(2)
      ga(2) = 2.0d0*x(2) + x(1)
      return
   30 ga(1) = cos(x(1))
      ga(2) = 0.0d0
      return
   40 ga(1) = 0.0d0
      ga(2) = sin(x(2))
      return
   50 go to (60,70,80) ka
   60 ga(1) = -exp(x(1)-x(2))
      ga(2) = exp(x(1)-x(2))
      return
   70 ga(1) = cosh(x(1)-1.0d0)
      ga(2) = 0.0d0
      return
   80 ga(1) = 0.0d0
      ga(2) = -1.0d0/x(2)
      return
   90 go to (100,110,120) ka
  100 a = sqrt((x(1)-x(3))**2+ (x(2)-x(4))**2)
      ga(1) = - (x(1)-x(3))/a
      ga(2) = - (x(2)-x(4))/a
      ga(3) = -ga(1)
      ga(4) = -ga(2)
      ga(5) = 0.0d0
      ga(6) = 0.0d0
      return
  110 a = sqrt((x(3)-x(5))**2+ (x(4)-x(6))**2)
      ga(1) = 0.0d0
      ga(2) = 0.0d0
      ga(3) = - (x(3)-x(5))/a
      ga(4) = - (x(4)-x(6))/a
      ga(5) = -ga(3)
      ga(6) = -ga(4)
      return
  120 a = sqrt((x(5)-x(1))**2+ (x(6)-x(2))**2)
      ga(1) = (x(5)-x(1))/a
      ga(2) = (x(6)-x(2))/a
      ga(3) = 0.0d0
      ga(4) = 0.0d0
      ga(5) = -ga(1)
      ga(6) = -ga(2)
      return
  130 do 140 i = 1,n
          ga(i) = -2.0d0*y(ka)*sin(y(ka)*x(i))/1.5d1
  140 continue
      return
  150 do 160 i = 1,n
          ga(i) = 0.0d0
  160 continue
      k = 0
      do 190 i = 1,5
          a = 0.0d0
          p = 0.0d0
          do 170 j = 1,n
              a = a + y(k+j)*x(j)** (1.0d0-y(i+80))
              p = p + y(k+j+40)*x(j)
  170     continue
          b = y(k+ka)/ (x(ka)**y(i+80)*a)
          do 180 j = 1,n
              c = y(k+j)* (1.0d0-y(i+80))/ (x(j)**y(i+80)*a)
              ga(j) = ga(j) + b* (y(k+j+40)-c*p)
  180     continue
          ga(ka) = ga(ka) - b*y(i+80)*p/x(ka)
          k = k + n
  190 continue
      return
      end subroutine tagu22

! subroutine tytim1                ms dos                     91/12/01
! portability : ms dos / ms fortran v.5.0
! 91/12/01 si : original version
!
! purpose :
!  get time in 100th of sec.
!
      subroutine tytim1(itime)
      integer :: itime
      real :: time
      call cpu_time(time)
      itime=1.0d2*time
      end subroutine tytim1
! subroutine tytim2                all systems                91/12/01
! portability : all systems
! 91/12/01 si : original version
!
! purpose :
!  print time elapsed.
!
      subroutine tytim2(itime)
      integer itime
      integer ihr,it,imin,isec
      call tytim1(it)
      it=it-itime
      ihr=it/(60*60*100)
      it=it-ihr*60*60*100
      imin=it/(60*100)
      it=it-imin*60*100
      isec=it/100
      it=it-isec*100
      write(6,10) ihr,imin,isec,it
   10 format(' time=',i2,':',i2.2,':',i2.2,'.',i2.2)
      end subroutine tytim2

end module tqsubs_module

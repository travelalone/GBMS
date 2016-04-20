program main
  implicit none
  integer,parameter :: n=1
  integer :: NRt=0,NRc=0
  integer :: stat,i,j,k,tr1,tr2,tr3
  real*8 :: area,s,sumarea,dot_r_n,temp
  real*8 :: cx,cy,cz,atomx(n),atomy(n),atomz(n),rdist,radii(n),ncx,ncy,ncz
  real*8 :: xtri(3),ytri(3),ztri(3),intsum(n)=0.0d0,bornradii(n),tpvec1(3),tpvec2(3),tpvec3(3)
  real*8,allocatable :: x(:),y(:),z(:)
  real*8 :: ra,rb,rc,theta
  character*50 :: buffer
  open(10,file="triangle")
  open(11,file="coordinate")
!============== calculate the line numbers ================
  do while(.true.)
     read(10,fmt="(A50)",iostat=stat) buffer
     if (stat==0) then
        NRt = NRt + 1
     else
        exit
     endif
  end do
  close(10)

  do while(.true.)
     read(11,fmt='(A)',iostat=stat) buffer
     if (stat==0) then
        NRc = NRc + 1
     else
        exit
     endif
  end do
  close(11)
  open(20,file="tinker.xyzr",action='read')
  
  do i = 1,n
     read(20,*,iostat=stat) atomx(i),atomy(i),atomz(i),radii(i)
  enddo
  
  close(20)
  
  allocate(x(NRc))
  allocate(y(NRc))
  allocate(z(NRc))
      
  !========== change the text file to bin file ===================
  open(11,file='coordinate')
  open(1,file='bin',access='direct',recl=24)
  do i=1,NRc
     read(11,*)x(i),y(i),z(i)
     write(1,rec=i)x(i),y(i),z(i)
  end do
  close(11)

!==================calculate the triangle area ===================
  open(10,file='triangle')
  sumarea=0.0d0
  do i=1,NRt
     read(10,*)tr1,tr2,tr3
     read(1,rec=tr1)xtri(1),ytri(1),ztri(1)
     read(1,rec=tr2)xtri(2),ytri(2),ztri(2)
     read(1,rec=tr3)xtri(3),ytri(3),ztri(3)

     !  Gravity Center of a Triangle
     cx = (xtri(1) + xtri(2) + xtri(3))/3.0d0
     cy = (ytri(1) + ytri(2) + ytri(3))/3.0d0
     cz = (ztri(1) + ztri(2) + ztri(3))/3.0d0

     tpvec1(1) = xtri(2) - xtri(1); tpvec1(2) = ytri(2) - ytri(1); tpvec1(3) = ztri(2) - ztri(1) 
     tpvec2(1) = xtri(3) - xtri(1); tpvec2(2) = ytri(3) - ytri(1); tpvec2(3) = ztri(3) - ztri(1)

     ! using the cross_product to find the normal vector 
     call cross_product(tpvec1,tpvec2,tpvec3)
     temp = sqrt(dot_product(tpvec3,tpvec3));

     tpvec3(1:3)=tpvec3(1:3)/temp
     ncx=tpvec3(1); ncy=tpvec3(2); ncz=tpvec3(3)
     tpvec3(1) = xtri(3) - xtri(2); tpvec3(2) = ytri(3) - ytri(2); tpvec3(3) = ztri(3) - ztri(2)
     ra = sqrt(dot_product(tpvec1,tpvec1)); rb = sqrt(dot_product(tpvec2,tpvec2)); rc=sqrt(dot_product(tpvec3,tpvec3));

     theta = (ra**2+rb**2-rc**2)/(2d0*ra*rb); 
     theta=sqrt(1.0d0-theta**2)
     area = 0.5d0*ra*rb*theta
     sumarea = sumarea + area

     if ( area == 0d0 ) then
        print *,"error ",i, ncx,ncy,ncz
        print "(a,3f8.3)","point 1 ",xtri(1),ytri(1),ztri(1)
        print "(a,3f8.3)","point 2 ",xtri(2),ytri(2),ztri(2)
        print "(a,3f8.3)","point 3 ",xtri(3),ytri(3),ztri(3)
        cycle
     end if

     do j = 1 ,n
        rdist = (cx - atomx(j))**2 + (cy - atomy(j))**2 + (cz -atomz(j))**2        
        dot_r_n = (cx - atomx(j)) * ncx + ( cy - atomy(j) )* ncy + (cz - atomz(j) )*ncz
        intsum(j) = intsum(j) + dot_r_n*area/rdist**2
     enddo
  end do
  
  do i = 1,n
     bornradii(i) = (4.0d0*3.14159d0)/intsum(i)
     write(*,'(I4,F7.3)')i,bornradii(i)
  end do
  
  print *, sumarea

end program main

subroutine cross_product(a,b,c)
  integer i,j,k
  real*8 a(3),b(3),c(3)
  c(1)=a(2)*b(3) - a(3)*b(2)
  c(2)=a(3)*b(1) - a(1)*b(3)
  c(3)=a(1)*b(2) - a(2)*b(1)
end subroutine cross_product



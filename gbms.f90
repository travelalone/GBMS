program main
  implicit none
  integer,parameter :: n=642
  integer :: NRt=0,NRc=0
  integer :: stat,i,j,tr1,tr2,tr3
  real*8 :: side1,side2,side3,area,s,sumarea,dot_r_n
  real*8 :: cx,cy,cz,atomx,atomy,atomz,rdist,radii,ncx,ncy,ncz
  real*8 :: xtri(3),ytri(3),ztri(3),intsum(n)=0.0d0,bornradii(n),nxt(3),nyt(3),nzt(3)
  real*8,allocatable :: x(:),y(:),z(:),nx(:),ny(:),nz(:)
  character*50 :: buffer
  open(10,file="triangle")
  open(11,file="coor_normal")
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

  allocate(x(NRc))
  allocate(y(NRc))
  allocate(z(NRc))
  allocate(nx(NRc))
  allocate(ny(NRc))
  allocate(nz(NRc))
  
  
  !========== change the text file to bin file ===================
  open(11,file='coor_normal')
  open(1,file='bin',access='direct',recl=48)
  do i=1,NRc
     read(11,*)x(i),y(i),z(i),nx(i),ny(i),nz(i)
     write(1,rec=i)x(i),y(i),z(i),nx(i),ny(i),nz(i)
  end do
  close(11)

!==================calculate the triangle area ===================
  open(10,file='triangle')
  sumarea=0.0d0
  do i=1,NRt
     read(10,*)tr1,tr2,tr3
     read(1,rec=tr1)xtri(1),ytri(1),ztri(1),nxt(1),nyt(1),nzt(1)
     read(1,rec=tr2)xtri(2),ytri(2),ztri(2),nxt(2),nyt(2),nzt(2)
     read(1,rec=tr3)xtri(3),ytri(3),ztri(3),nxt(3),nyt(3),nzt(3)
     side1 = sqrt( (xtri(1)-xtri(2)) **2 + (ytri(1)-ytri(2)) **2 + (ztri(1)-ztri(2)) **2 )
     side2 = sqrt( (xtri(1)-xtri(3))**2 + (ytri(1)-ytri(3))**2 + (ztri(1)-ztri(3))**2 )
     side3 = sqrt( (xtri(2)-xtri(3))**2 + (ytri(2)-ytri(3))**2 + (ztri(2)-ztri(3))**2 )
     s =( side1 + side2 + side3) /2.0d0
     area = sqrt(s*(s-side1)*(s-side2)*(s-side3))
     sumarea = sumarea + area

     !  Gravity Center of a Triangle
     cx = (xtri(1) + xtri(2) + xtri(3))/3.0d0
     cy = (ytri(1) + ytri(2) + ytri(3))/3.0d0
     cz = (ztri(1) + ztri(2) + ztri(3))/3.0d0
     ! the normal of gravtity of a triangle
     ncx = (nxt(1)+nxt(2)+nxt(3))/3.0d0
     ncy = (nyt(1)+nyt(2)+nyt(3))/3.0d0
     ncz = (nzt(1)+nzt(2)+nzt(3))/3.0d0
     
     open(20,file="tinker.xyzr",action='read')
     do j = 1 ,n
        read(20,*,iostat=stat) atomx,atomy,atomz,radii
        rdist = (cx - atomx)**2 + (cy - atomy)**2 + (cz -atomz)**2
        dot_r_n = (cx - atomx) * ncx + ( cy - atomy )* ncy + (cz - atomz )*ncz
        intsum(j) = intsum(j) + 1.0d0/(rdist**2)*area*dot_r_n
     enddo
     close(20)
  end do
  
  do i = 1,n
     bornradii(i) = 1.0d0/(intsum(i)/(4*3.14159) )
     write(*,'(I4,F7.3)')i,bornradii(i)
  end do
  print *, sumarea

end program main

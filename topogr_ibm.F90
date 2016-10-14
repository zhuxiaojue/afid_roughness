!
!     This routine finds the indices in the computational grid
!     close to the physical location of the body.
!
      subroutine topogr
      use param
      use decomp_2d, only: xstart,xend
      use ibm_param
      implicit none
      integer :: i,j,k,l,kstartp
      integer :: km,kp,jm,jp,im,ip

      real    :: xe, xem, xep
      real    :: ye, yem, yep
      real    :: ze, zem, zep
      real    :: delta1x, delta2x

      infig=1
      q1bo=0.d0
      q2bo=0.d0
      q3bo=0.d0
      densb=0.d0 
      allocate(plth1(1:nzm,1:nym))
      allocate(plth2(1:nzm,1:nym))

      do i=1,nzm
       do j=1,nym
!        plth1(j,i)=abs(dcos(10*pi*zm(i)/zlen) &
!                     +dcos(10*pi*ym(j)/ylen))*0.01

         plth1(j,i)=abs(dcos(10*pi*ym(j)/ylen))*0.05
!        plth2(j,i)=abs(dcos(10*pi*zm(i)/zlen)  &
!                     +dcos(10*pi*ym(j)/ylen))*0.01
         plth2(j,i)=abs(dcos(10*pi*ym(j)/ylen))*0.05
       end do
      end do


      npunx=0
      npuny=0
      npunz=0
!
!     IDENTIFICATION OF THE GRID POINTS IN THE BODY
!     (BOUNDARY + INNER PART)
!
!     X_3 vertical direction
!     X_2 radial direction
!     X_1 azimuthal direction
!
      
        
        allocate(forclo(1:nx,xstart(2):xend(2),xstart(3):xend(3)))

      do l = 1,3 !{ start do over the 3 velocity components
      n=0

!     l = 1   Q_1 vel. component
!     l = 2   Q_2 vel. component
!     l = 3   Q_3 vel. component
!

      do i=xstart(3),xend(3)
       do j=xstart(2),xend(2)
        do k=1,nxm
         km=kmv(k)
         kp=kpv(k)
         xe=xm(k)
         xem=xm(km)
         xep=xm(kp)
         if(l.eq.3) then
           xe=xc(k)
           xem=xc(km)
           xep=xc(kp)
         end if
!
!    SOLID PART
!
         if((xe.lt.plth1(j,i))) then
          n=n+1
          indgeo(l,n,1)=i
          indgeo(l,n,2)=j
          indgeo(l,n,3)=k
          indgeoe(l,n,1)=i
          indgeoe(l,n,2)=j
          indgeoe(l,n,3)=k
          distb(l,n)= 0.
         elseif(xe.gt.(alx3-plth2(j,i))) then

          n=n+1
          indgeo(l,n,1)=i
          indgeo(l,n,2)=j
          indgeo(l,n,3)=k
          indgeoe(l,n,1)=i
          indgeoe(l,n,2)=j
          indgeoe(l,n,3)=k
          distb(l,n)= 0.
        
!    LOWER FLUID/PLATE BOUNDARY
!
        elseif((xe.ge.plth1(j,i)).and.(xem.lt.plth1(j,i))) then
          n=n+1
          indgeo(l,n,1)=i
          indgeo(l,n,2)=j
          indgeo(l,n,3)=k
          indgeoe(l,n,1)=i
          indgeoe(l,n,2)=j 
          indgeoe(l,n,3)=kp
          delta1x=(xep-xe)
          delta2x=(xe-plth1(j,i))
          distb(l,n)= delta2x/(delta1x+delta2x)

!
!    UPPER FLUID/PLATE BOUNDARY
!
      elseif((xe.le.(alx3-plth2(j,i))).and.(xep.gt.(alx3-plth2(j,i)))) &
         then
           n=n+1
           indgeo(l,n,1)=i
           indgeo(l,n,2)=j
           indgeo(l,n,3)=k
           indgeoe(l,n,1)=i
           indgeoe(l,n,2)=j 
           indgeoe(l,n,3)=km
          delta1x=(xe-xem)
           delta2x=((alx3-plth2(j,i))-xe)
           distb(l,n)= delta2x/(delta1x+delta2x)
          end if

          end do
        end do
      end do

      if(l.eq.1) then
        if(n.gt.mpun)  &
       write(*,*) 'Dim max di indgeot e'' stata superata n=',n
        npunz= n
        write(6,332)npunz
 332  format(5x,'For Q_1 N ='i7)
      end if
      if(l.eq.2) then
        if(n.gt.mpun)   &
       write(*,*) 'Dim max di indgeor e'' stata superata n=',n
        npuny= n
        write(6,331)npuny
 331  format(5x,'For Q_2 N ='i7)
      end if
      if(l.eq.3) then
        if(n.gt.mpun)  &
       write(*,*) 'Dim max di indgeoz e'' stata superata n=',n
        npunx= n
        write(6,330)npunx
 330  format(5x,'For Q_3 N ='i7)
      end if
      end do   !} end do over the 3 velocity components 

  

!
!     INDICES FOR TEMPERATURE
!
      n=0
      forclo =0.0d0
!

      do i=xstart(3),xend(3)
       do j=xstart(2),xend(2)
        do k=1,nxm
        km=kmv(k)
        kp=kpv(k)
        xe=xc(k)
        xem=xc(km)
        xep=xc(kp)
!
!    SOLID PART
!
            if(xe.lt.plth1(j,i)) then
             n=n+1
             indgeot(n,1)=i
             indgeot(n,2)=j
             indgeot(n,3)=k
             indgeoet(n,1)=i
             indgeoet(n,2)=j
             indgeoet(n,3)=k
             distbt(n)= 0.
             temb(n) = tempbp(j,i)
             forclo(k,j,i) = 1.
              
            elseif(xe.gt.(alx3-plth2(j,i))) then
             n=n+1
             indgeot(n,1)=i
             indgeot(n,2)=j
             indgeot(n,3)=k
             indgeoet(n,1)=i
             indgeoet(n,2)=j
             indgeoet(n,3)=k
             distbt(n)= 0.
             temb(n) = temptp(j,i)
             forclo(k,j,i) = 1.
!            end if

!
!    LOWER FLUID/PLATE BOUNDARY
!
            elseif((xe.ge.plth1(j,i)).and.(xem.lt.plth1(j,i))) then
              n=n+1
              indgeot(n,1)=i
              indgeot(n,2)=j
              indgeot(n,3)=k
              indgeoet(n,1)=i
              indgeoet(n,2)=j 
              indgeoet(n,3)=kp
              delta1x=(xep-xe)
              delta2x=(xe-plth1(j,i))
              distbt(n)= delta2x/(delta1x+delta2x)
              temb(n) = tempbp(j,i)
!
!    UPPER FLUID/PLATE BOUNDARY
!
      elseif((xe.le.(alx3-plth2(j,i))).and.(xep.gt.(alx3-plth2(j,i)))) &
          then
            n=n+1
            indgeot(n,1)=i
            indgeot(n,2)=j
            indgeot(n,3)=k
            indgeoet(n,1)=i
            indgeoet(n,2)=j 
            indgeoet(n,3)=km
            delta1x=(xe-xem)
            delta2x=((alx3-plth2(j,i))-xe)
            distbt(n)= delta2x/(delta1x+delta2x)
            temb(n) = temptp(j,i)
           end if

          end do
        end do
      end do
        if(n.gt.mpun) &
       write(*,*) 'Dim max di indgeote e'' stata superata n=',n
        npunte= n
        write(6,329)npunte
 329  format(5x,'For Temperature N ='i7)
       if(allocated(plth1)) deallocate(plth1)
       if(allocated(plth2)) deallocate(plth2)

      return
      end


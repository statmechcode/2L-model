      program car_lan                                           
c     *************************************************************       
c     * Nakanishi algorithm for Carlson-Langer                    *       
c     * Last revision Dic 2009                                    *       
c     *************************************************************       
	implicit none                                                           
                                                                          
c       -----------------------------------------------------             
c       First Variables Declaration                                       
                                                                          
        integer Lx,Ly,final_time,realtot,loopmax,nimax                
        real*8 vp,sigmath,sigmac,D
        real k0,kk,theta,zetab,lambda
        real*8 eta

        CHARACTER(len=100)::nomef,stringa !Nome files uscita
        CHARACTER(len=4)::number_Ly      
        CHARACTER(len=4)::number_Lx
        CHARACTER(len=5)::number_sigmath
        CHARACTER(len=5)::number_theta   
        CHARACTER(len=6)::number_k0
        CHARACTER(len=6)::number_zetab
                                                                          
c........input parameters.......................................          
                                                                          
c       Remember to change random number seed below !                     
c       Check that the interface does not reach the boundary !            
                                                                          
        parameter (Lx=100)        ! system size along x 
        parameter (Ly=100)        ! system size along y
        parameter (final_time=8E8)      ! total number of iterations      
        parameter (realtot=1)
        parameter(vp=.01)               !driven velocity 
        parameter(sigmath=5.0)           !st.dev. in the threshold 
        parameter(k0=.02)     ! dissipation (no dissipation for alpha=1)  
        parameter(theta=0.5)    ! dissipation (no dissipation for alpha=1)
        parameter(kk=0.3)  ! dissipation (no dissipation for alpha=1)
        parameter(eta=0.0000001) ! dissipation (no dissipation for alpha=1)
        parameter(zetab=2)          ! st. dev. of afterslip values 0.57
        parameter(loopmax=1E8)
        parameter(lambda=0.01)
        parameter(nimax=int(0.5*loopmax))
                                                                          
c***********************************************************
        integer realiz,i,j,m,ni,kscrive,nio,nif,jlist(loopmax)
        integer nevent,nn_i(0:lx+1,0:Ly+1,4),nn_j(0:lx+1,0:Ly+1,4)
        integer every,itime,mag1,n2(1000000),n1(1000000),nevent2
        integer seed1,nm(0:1000000),imag,flags,ran,ilist(loopmax)
        real*8 f(0:Lx+1,0:Ly+1),ff, fth(0:Lx+1,0:Ly+1),diff
        real*8 g(0:Lx+1,0:Ly+1),diffp,fff,sigma,mags,zt(loopmax)
        real ran2,bin0,rr1,rr2,sum,isoglia,amp,tscrive
        real*8 ffmin,deltat,deltatm,z,time,aveg2,temp
   
        real rate,ratep,qmax,Delta,size,gamma
        integer mprima,iiconta,mmax,imax,jmax,jloop,track(Lx,Ly)
        integer iloop,iepi,jepi,final_loop,u(0:Lx,0:Ly+1),i1,j1,ii
        integer bor(Lx,Ly),iflag(loopmax),nseq,k,indx(0:Lx*Ly)
        integer indy(0:Lx*Ly)
        integer*8 mag

        !print*,kk,theta,kk*theta,kk*(1-theta)
        !pause
        write(number_theta,63)theta
	write(number_sigmath,63)sigmath
        write(number_Ly,73)Ly
        write(number_Lx,73)Lx
        write(number_k0,66)k0
        write(number_zetab,66)zetab

        nomef='fast_Lx_'//number_Lx//'Ly_'//number_Ly//
     &       'theta_'//number_theta//'sigma_'//number_sigmath//
     &       'dh_'//number_zetab//'.dat'
        
 61     FORMAT(f5.3)      
 62     FORMAT(f4.2)
 63     FORMAT(f5.3)
 73     FORMAT(i4.4)
 74     FORMAT(i3.3)
 64     FORMAT(f8.6)
 66     FORMAT(f6.3)
        
      open(2101,file=nomef,status='unknown')
c      open(800,file='soglia.dat',status='unknown')

      mags=0d0
      aveg2=0d0
      amp=-0.15
      seed1=-2128719
      gamma=1.0d0

      tscrive=20000
      kscrive=0
      
          
c.......initial settings......................................            
c*************************************************************************
      do j=1,ly-1,1
         do i=1,Lx
            nn_i(i,j,1)=i
            nn_j(i,j,1)=j-1
            nn_i(i,j,2)=i-1
            nn_j(i,j,2)=j
            nn_i(i,j,3)=i
            nn_j(i,j,3)=j+1
            nn_i(i,j,4)=i+1
            nn_j(i,j,4)=j
            u(i,j)=1
         enddo
      enddo
      do i=1,Lx-1
         u(i,0)=0
         u(i,Ly)=0
      enddo
      do j=1,Ly
         u(0,j)=0
         u(Lx,j)=0
      enddo

c******************************************************************                                                                          
c     ----- Beginning of loop on realizations                             
c**********************************************************************
      
        do realiz=1,realtot
           nevent=0
           nseq=0
           do iloop=1,loopmax
              ilist(iloop)=-1
              iflag(iloop)=0
           enddo
           ni=0   !number of instable sites

           do j=1,Ly
              do i=1,Lx
                 track(i,j)=1
              enddo
           enddo
           
           
c******************************************************************
C           Initial configuration
c*******************************************************************
           do j=1,Ly
              do i=1,Lx
                 fth(i,j)=4000
                 f(i,j)=0
              enddo
           enddo
          
!           read(6002,'(a)')stringa
           do j=2,Ly-1
              do i=2,Lx-1
c$$$                 read(6002,*,end=96)fth(i,j),f(i,j),g(i,j)
c$$$                  if(fth(i,j).le.f(i,j))then
c$$$                    temp=f(i,j)
c$$$                    f(i,j)=fth(i,j)
c$$$                    fth(i,j)=temp
c$$$                 endif
                 rr1=ran2(seed1)
                 rr2=ran2(seed1)
                 fth(i,j)=1+sigmath*sqrt(-2*log(rr1))*cos(2*3.1415*rr2)
                 f(i,j)=rr1*fth(i,j)*u(i,j)
                 g(i,j)=0
              enddo
           enddo
 96        continue


c**************************************************************
          time=0
          kscrive=0
          do while(0.000001*time.le.final_time)
    
c******************************************************************
c              Evolution due to drive
c*******************************************************************
             ni=0
             ffmin=1E12
             do j=2,Ly-1
                do i=2,Lx-1
                   ff=fth(i,j)-f(i,j)
                   if(ff.lt.ffmin.and.ff.gt.0)then
                      iepi=i
                      jepi=j
                      ffmin=ff
                   endif
                enddo
             enddo
             !print*,ffmin
!pause
             
             time=time+ffmin/(k0*vp)
             !print*,time,ffmin
              do j=2,Ly-1
                 do i=2,Lx-1
                    f(i,j)=f(i,j)+ffmin
                 enddo
              enddo

              nseq=nseq+1

            ! print*,iepi,jepi,f(iepi,jepi),time
              !pause
              delta=1


c**********************************************************************
c              Avalanche
c**********************************************************************
 77           nio=-1
	      nif=1
	      ni=0
              iloop=0
              final_loop=2
              ilist(1)=iepi
              jlist(1)=jepi

	      do while (ni.gt.nio)
		 ni=ni+nif
		 nio=ni-nif

                 if(ni.ge.nimax)then
                    if(nif.ge.nimax)print*,'ATTENZIONE'
                    do iloop=nio+1,ni
                       ilist(iloop-nio)=ilist(iloop)
                       jlist(iloop-nio)=jlist(iloop)
                    enddo
                    ni=ni-nio
                    nio=0
                 endif
                 
		 nif=0
c		 print*,time,ni
		 do iloop=nio+1,ni
		    i=ilist(iloop)
		    j=jlist(iloop)
                    
c                    pause
		    if(f(i,j)+g(i,j)+1E-12.ge.fth(i,j))then
!!!!!!!!!!!!!!!!!!!!1E-12 added to avoid problems with roundings !!!!!!
                       !bor(i,j)=1
                       iflag(iloop)=1
 7                     continue
		       rr1=ran2(seed1)
		       rr2=ran2(seed1)
                       z=zetab  !-zetab*log(ran2(seed1))
                      !z=max(zetab,(fth(i,j)-f(i,j)-g(i,j))*1./(4*kk+k0))
		       zt(iloop)=zt(iloop)+z
                       !print*,fth(i+1,j),f(i+1,j),g(i+1,j)
                       !pause
		       g(i,j)=g(i,j)-4*KK*z*theta
		       f(i,j)=f(i,j)-(4*KK*(1-1.*theta)+K0)*z
                       fth(i,j)=1
     &                 +sigmath*sqrt(-2*log(rr1))*cos(2*3.1415*rr2)
		       mags=mags+z
		       if(f(i,j)+g(i,j).ge.fth(i,j))then
			  nif=nif+1
			  ilist(ni+nif)=i
			  jlist(ni+nif)=j
c			  goto 7
		       endif
                       
                       if(track(i,j).eq.1)then
                          mag=mag+1
                          track(i,j)=0
                          indx(mag)=i
                          indy(mag)=j
                       endif

		    endif
		 enddo
                 
		 do iloop=nio+1,ni
                    if(iflag(iloop).eq.1)then
                       i=ilist(iloop)
                       j=jlist(iloop)
                       do ii=1,4
                          i1=nn_i(i,j,ii)
                          j1=nn_j(i,j,ii)
                          if(zt(iloop).gt.0)then
                             !f(i1,j1)=f(i1,j1)+K1*zt(iloop)
                             f(i1,j1)=f(i1,j1)+KK*(1-1*theta)*z!*zt(iloop)
c     &                            +K0*zt(iloop)*(0.0025d0)
                             g(i1,j1)=g(i1,j1)+KK*theta*z!t(iloop)*gamma
                             if(f(i1,j1)+g(i1,j1).ge.fth(i1,j1))then
                                nif=nif+1
                                ilist(ni+nif)=i1
                                jlist(ni+nif)=j1
                             endif

                          endif

                       enddo
                    endif
                    iflag(iloop)=0
                    zt(iloop)=0
                 enddo

              enddo
              
              if(mags.ge.10)write(2101,*)time*0.000125,int(mag),iepi,
     &             jepi,mags,nseq
              call flush(2101)
              
              do k=1,mag
                 track(indx(k),indy(k))=1
              enddo
              
              !if(mags.gt.1E7)then
              !   print*,time,mag,iepi,jepi,ni
              !endif
              
              ni=0
              mag=0
              mags=0
             

c*************************************************************************
c              Relaxation 
c*************************************************************************

             deltatm=1E12
             flags=0
             do j=1,Ly-1
                do i=1,Lx-1
                   if(g(i,j).lt.0.and.fth(i,j)-f(i,j).lt.0)then
                      temp=(fth(i,j)-f(i,j))/g(i,j)
                      !print*,temp
                      !pause
                      deltat=exp((1-temp)/lambda)-1
                      if(deltat.gt.0.and.deltat.le.deltatm)then
                         deltatm=deltat
                         iepi=i
                         jepi=j
                         flags=1
                      endif
                      !print*,g(i,j)
                   endif
                enddo
             enddo
             !stop
      
             if(flags.eq.1)then

                !time=time+deltatm*(eta/(kk*theta))
                delta=1-lambda*log(1+deltatm)
                do j=1,Ly-1
                   do i=1,Lx-1
c                      g(i,j)=g(i,j)*(1-deltatm**(1/7))
                      g(i,j)=g(i,j)*(1-lambda*log(1+deltatm))                  
                   enddo
                enddo

                goto 77
             else
                do j=1,Ly-1
                   do i=1,Lx-1
                      g(i,j)=0
                   enddo
                enddo
             endif
             
             
          enddo                 !loop over time

       enddo                    !loop over realizations

                                                                          
c     ------------------------------------------------------------        
       stop                                                               
       end                                                                
	                                                                         
      FUNCTION ran2(idum)                                                 
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV         
      REAL ran2,AM,EPS,RNMX                                               
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,      
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,         
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)                    
      INTEGER idum2,j,k,iv(NTAB),iy                                       
      SAVE iv,iy,idum2                                                    
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/                            
      if (idum.le.0) then                                                 
        idum=max(-idum,1)                                                 
        idum2=idum                                                        
        do 11 j=NTAB+8,1,-1                                               
          k=idum/IQ1                                                      
          idum=IA1*(idum-k*IQ1)-k*IR1                                     
          if (idum.lt.0) idum=idum+IM1                                    
          if (j.le.NTAB) iv(j)=idum                                       
11      continue                                                          
        iy=iv(1)                                                          
      endif                                                               
      k=idum/IQ1                                                          
      idum=IA1*(idum-k*IQ1)-k*IR1                                         
      if (idum.lt.0) idum=idum+IM1                                        
      k=idum2/IQ2                                                         
      idum2=IA2*(idum2-k*IQ2)-k*IR2                                       
      if (idum2.lt.0) idum2=idum2+IM2                                     
      j=1+iy/NDIV                                                         
      iy=iv(j)-idum2                                                      
      iv(j)=idum                                                          
      if(iy.lt.1)iy=iy+IMM1                                               
      ran2=min(AM*iy,RNMX)                                                
      return                                                              
      END                                                                 
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          

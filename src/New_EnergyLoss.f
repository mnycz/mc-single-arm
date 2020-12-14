	subroutine New_EnergyLoss(len,dens,zeff,aeff,epart,mpart,typeflag,elossflag,Eloss)

	implicit none

	real*8 thick,len,dens,zeff,aeff,epart,mpart,Eloss
	real*8 x,chsi,lambda,gauss1,Eloss_mp,gamma
	real*8 denscorr,CO,hnup,log10bg,I,beta,Eloss_mp_new
	real*8 New_Eloss,lxi,lhbarw,lhbarwsq,Delta_p,num
	real*8  result,RANLAN
	integer typeflag          !1=normal eloss (picked from distribution)
                                  !2=min eloss
	                          !3=max eloss
                                  !4=most probable eloss

	integer elossflag      !1 = SIMC Eloss
	                       !2 Most probable energy loss
        integer numerr
        data numerr /0/

	real*8 me
        parameter(me=0.51099906)
c	real ld
	common /landau/ ld(0,1)
	real ld 
	save /landau/

	thick = len*dens
	gamma=epart/mpart
        beta = sqrt(1.-1./gamma**2)
c	write(6,*) len,zeff,aeff,epart,mpart,thick
	if (elossflag.eq.1) then 
	if(zeff.eq.1) then	!Ionization potential in MeV
	  I = 21.8e-06
	else
	  I = (16.*zeff**0.9)*1.0e-06
	endif

	hnup = 28.816e-06*sqrt(dens*zeff/aeff) !plasma frequency
	log10bg = log(beta*gamma)/log(10.)
	CO=log(hnup)-log(I)+0.5
c	write(6,*) CO
C DJG Get density effect correction (I got this from JV).

	if(log10bg.lt.0.) then
	  denscorr=0.
	elseif(log10bg.lt.3.) then
	  denscorr=CO+log(10.)*log10bg+abs(CO/27.)*(3.-log10bg)**3
	elseif(log10bg.lt.4.7) then
	  denscorr=CO+log(10.)*log10bg
	else
	  denscorr=CO+log(10.)*4.7
	endif
c	write(6,*) 'denscor',denscorr
	if (thick.le.0.) then
	  Eloss = 0.
	else
	  Eloss_mp = 0.1536e-03 * zeff/aeff * thick * ( 19.26 +
     &          log(thick/dens) )
	  Eloss_mp_new = 0.1536e-03 * zeff/aeff *thick/beta**2* (
     &          log(me/I**2) + 1.063 + 2.*log(gamma*beta) + 
     &		log(0.1536*zeff/aeff*thick/beta**2)-beta**2-denscorr)
c	  write(6,*) 'ELOSS',Eloss_mp,Eloss_mp_new 
! ........ convert to MeV, the unit of choice in THIS program
! ........ (cf. EVCOIN where GeV prevail)
	  Eloss_mp = Eloss_mp_new*1000.
	  chsi = 0.307075/2.*zeff/aeff*thick/beta**2
	  if(typeflag.eq.1)then
	    x=abs(gauss1(10.0e0))
	  elseif(typeflag.eq.2)then
	    x=3
	  elseif(typeflag.eq.3)then
	    x=0.0067
	  elseif(typeflag.eq.4)then
	    x=1
          endif
	  if(x.gt.0.0) then
	    lambda = -2.0*log(x)
	  else
	    lambda = 100000.
	  endif
	  Eloss = lambda*chsi+eloss_mp
c	  write(6,*) 'Eloss',Eloss
	endif
        if (eloss.gt.(epart-mpart)) then
	   eloss=(epart-mpart)-0.0000001
	   numerr=numerr+1
	   if (numerr.le.10) then
	      write(6,*) 'Eloss>Total KE; forcing Eloss=KE'
	      if (numerr.eq.10) write(6,*) '     FURTHER ELOSS ERRORS SUPPRESSED'
	   endif
        endif 
	
c swith to most probable energy loss: PDG 2019 equation 3.12 (or high energy approx 3.13)
	elseif(elossflag.eq.2) then
	   lxi =(0.307075/2.0)*(zeff/aeff)*(thick/beta)
c	   write (6,*) 'Xi',lxi
	   lhbarw =28.816*sqrt(dens*(zeff/aeff))*1e-6 ! MeV                                                                                                                            
           lhbarwsq = lhbarw*lhbarw    ! MeV^2                
	   num = 2.*me*lxi
	   if(lhbarw.gt.0)then
	      Delta_p =lxi*(log(num/lhbarwsq)+0.2) 
	      write(6,*) Delta_p,lxi
c	      result = ld(Delta_p,lxi)
c	      result=RANLAN(Delta_p,lxi)
	      result = RANLAN(1)
c	      Eloss=result
	      write(6,*) 'Eloss',result
	      endif	     
	elseif(elossflag.gt.2) then
	   write(6,*) 'Eloss routines must be 1 or 2' 	   
	endif
	return
	end

       subroutine sigmavdirect(relicname, sigmaV, sigmaV_ch)
        implicit none
      
#include "decl.h"



        double precision result_tmp(NCOMP)
        integer flags
        integer helicit
        external helicit
        integer i
        double precision sqrtS
        double precision costh
        double precision p12,p34
        character name*50
        double precision sigmaV
        double precision sigmaV_ch(25)
        character relicname*2
        character (len = 10), dimension (25) :: ch= [character(len=10) ::"ebare","nebarne","mbarm","nmbarnm", "nlbarnl", "lbarl", "sbars", "cbarc", "ubaru", "dbard", "tbart", "bbarb", "hh", "hhh","hh3", "hz", "hhhh", "hhh3", "hhz", "h3h3", "h3z", "hcbarhc", "hcbarw", "zz", "wbarw" ]
        integer, dimension(25) :: Sf34=[1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,2,1,1,2,1,1,1,2,1]
        double precision  particule3_masse(25), particule4_masse(25)
        integer h1, h2,h3,h4
        integer hel(25)
        double precision j,x
        if (relicname=='o1') then
           
            sqrtS=2.*relicmass
 
            costh=0.
            h4=helicit(FERMION,FERMION,VECTOR, VECTOR)
            h1=helicit(FERMION,FERMION,FERMION,FERMION)
            h2=helicit(FERMION,FERMION,SCALAR, SCALAR)
            h3=helicit(FERMION,FERMION,SCALAR, VECTOR)
            hel=[h1,h1,h1,h1,h1,h1,h1,h1,h1,h1,h1,h1,h2,h2,h2,h3,h2,h2,h3,h2,h3,h2,h3,h4,h4]
!            open (unit=11,file="/home/robbins/Packages/cross_section_stat/running_mass.txt",action="read",status="old")
!            
!           
!            read (unit=11,fmt=*) Ms, Mb, Mt, MW   ! lecture du fichier
!                           ! ecriture a l'ecran de ce
!                                    ! qu'on a lu
!             
!            
!            close(unit=11) 
            



            particule3_masse=[Me, 0.D0, Mm, 0.D0, 0.D0, Ml, Ms, Mc, Mu, Md, Mt, Mb, Mh, Mh, Mh, Mh, Mhh, Mhh, Mhh, Mh3, Mh3, Mhc, Mhc, MZ, MW]
            particule4_masse=[Me, 0.D0, Mm, 0.D0, 0.D0, Ml, Ms, Mc, Mu, Md, Mt, Mb, Mh, Mhh, Mh3, Mz, Mhh, Mh3, Mz, Mh3, Mz, Mhc, MW, MZ, MW]

          
    
        
            sigmaV=0
            flags=1
            do i = 1, NCOMP
                result_tmp(i)=0.D0
!         
            enddo
            call init_kin(relicmass,Mu,relicmass,Mu,sqrtS,costh,p12,p34)
            name="o1uo1u"
            call squared(trim(name),len_trim(name),result_tmp, h1, flags)
            print *,"sig()=", result_tmp(1)
            
!            do i=1, 25
!                if(sqrtS .gt. (particule3_masse(i)+particule4_masse(i))) then
!                name=trim(relicname)//trim(relicname)//trim(ch(i))
!                sigmaV_ch(i)=0.D0
!                call init_kin(relicmass,relicmass,particule3_masse(i),particule4_masse(i),sqrtS,costh,p12,p34)
!                call squared(trim(name),len_trim(name),result_tmp, hel(i), flags)
!                sigmaV_ch(i)=result_tmp(1)/(64.*pi*relicmass**2)*sqrt(1-(particule3_masse(i)**2+particule4_masse(i)**2)/(2.*relicmass**2) + ((particule3_masse(i)**2-particule4_masse(i)**2)**2)/(16*relicmass**4))*3.8937966e+8*2.997900e-26/2./Sf34(i)
!                print *,"sig(",i,")=", result_tmp(1)
!                sigmaV=sigmaV+sigmaV_ch(i)
!                endif
!            enddo
            
            
         
        
        else 
            print*, 'LSP is  not neutalino, cross-section cannot be computed'
        
        endif
  
        end

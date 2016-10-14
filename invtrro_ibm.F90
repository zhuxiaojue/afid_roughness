       forclo=1.d0 
       if(infig.gt.0) then
         usaldto = 1./aldto
         do n=1,npunte
           ic=indgeot(n,1)
           jc=indgeot(n,2)
           kc=indgeot(n,3)
           forclo(kc,jc,ic)=0.d0
           ie=indgeoet(n,1)
           je=indgeoet(n,2)
           ke=indgeoet(n,3)
           dense=((al*dt+aldto)*temp(ke,je,ie)-al*dt*densb(n))*usaldto
           rhs(kc,jc,ic) = -temp(kc,jc,ic) + dense*distbt(n)  &
                         +(1.-distbt(n))*temb(n)
           densb(n)= temp(ke,je,ie)
         end do
        end if
       
         
         call SolveImpEqnUpdate_Temp_ibm

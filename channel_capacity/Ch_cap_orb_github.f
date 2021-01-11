             Implicit Real*8(A-H,O-Z)
             Real*8  DM(2,2),  L_n(2,2)
             Real*8 LE, LE_old
c--------------            
c         Computes the channel capacity for selected orbits
c     versus  iteration number n for  corrugation eps fixed
c             
c     Input parameters: iteration number N_max, grid points number   N_grid  ( x only )
c                       corrugation amplitude eps  
c--------------
             open(10,file='Chan_cap.dat')
             pi=4*datan(1.d0)
             duepi=2*pi
            
c---
             N_grid=20
             N_max=500  
             eps=0.25d0
c---
             A_grid=N_grid
             
                     
             Ch_LE=0    ! channel capacity   <log(E_n^L)/n>
             Ch_RE=0    ! channel capacity   <log(E_n^R)/n>
             Ch_REM=0   ! channel capacity   <log(E_n^REM)/n>
            
             
             do Ix=1,N_grid               !       loop in x
                x0= pi*(-1+  2*Ix/A_grid) !      -pi < x < pi             
             write(*,*) Ix
            
                   call conv(x0,v0,xx0,y0)
                   write(15,'(I5,2g13.5)')   Ix,xx0,v0
                   write(20,*) '  '
                   write(20,*)  '####',  Ix,Iv, xx0,v0
                   write(20,*) '  '
                   
                   x=x0
                   v=v0

             L_n(1,1)=1
             L_n(1,2)=0
             L_n(2,1)=0
             L_n(2,2)=1                    !   Initialization  DM 
             

             LE_old=sqrt(2.d0)      ! Initialization for RE
             RE2=0

                
             

                   
                   do n=1,N_max   ! ----->   starts loop iterations n
                      x_old=x
                      v_old=v
                      call mappa(x,v,eps,tau,zer)

                      call conv(x,v,xx,y) !  xx=x/2 pi mod 1   y= acos(v)/duepi               
                      write(20,'(I5,12g13.5)') n,xx,v  ! writes orbits
                      
                call conv(x,v,xx,y) !  xx=x/2 pi mod 1   y= acos(v)/duepi               
                call tan_map(x_old,v_old,x,v,eps,tau,DM)                
                call mult(DM,L_n)

             Tr_n=L_n(1,1)**2 + L_n(1,2)**2 + L_n(2,1)**2 + L_n(2,2)**2    ! Tr( (L_n)^T L_n)
                             
                   LE= sqrt(Tr_n)                 ! error LE_n
                   RE2=RE2 + LE_old**2 + LE**2
                   RE=sqrt(RE2)                   ! error RE_n
                   LE_old=LE
                   an=n

                    write(20+Ix,'(I5,7g13.5)') n,LE,RE, dlog(LE)/an, 
     *              dlog(RE)/an 

         
                   enddo                     ! <---------- end loop n 
      
                   LE= sqrt(Tr_n)

c-------------    Errore REM 
                v=-v                             ! inversion
                do n=1,N_max
                call mappa(x,v,eps,tau,zer)      ! backawd iteration
                enddo
                v=-v
                call conv (x,v,xx,y)
                REM=sqrt((xx-xx0)**2+(v-v0)**2)/1E-16  ! error REM  for  n=N_max
c-------------
                   
             enddo   !  <----------------    end loop in  x 
       
 
             
             return
             end

           subroutine mappa  (x,v,eps,tau,zer)
             Implicit Real*8(A-H,O-Z)
c------------------
c             Reflection map: see formula  (9)  on the paper
c             input  x=x_{n-1}, v=v_{n-1}   output  x=x_n  v= v_n 
c------------------             
             duepi=8*datan(1.d0)
             precis=1.E-15    ! precision  of  the bisection method
             vz=sqrt(1-v*v)
c------------
c           bisection to compute  tau
c------------                
             tau_a=0.5/vz
             tau_b=1.5/vz
             G_a=1+eps*f(x+tau_a*v)-tau_a*vz
             G_b=1+eps*f(x+tau_b*v)-tau_b*vz
             if(G_a*G_b.gt.0)  write(*,*) ' ***** no zero ****'
             
             do iter=1,60
                tau_c=(tau_a+tau_b)/2
                G_c=1+eps*f(x+tau_c*v)-tau_c*vz
                if(G_a*G_c.lt.0) then
                   tau_b=tau_c
                   G_b=G_c
                else if(G_b*G_c.lt.0) then
                   tau_a=tau_c
                   G_a=G_c
                endif
                diff=(tau_b-tau_a)/(tau_b+tau_a)
                if(diff.le.precis) go to 1
             enddo
 1           tau=(tau_a+tau_b)/2
             zer=1+eps*f(x+tau*v)-tau*vz

c------------
c           eveluation of the map
c------------
             ffp=fp(x+tau*v)
             den= 1+(eps*ffp)**2
             v_new=v- 2*( (eps*ffp)**2*v -eps*ffp*vz)/den
             x_new=x+tau*v+ tau*vz* v_new/sqrt(1-v_new**2)
             x=x_new
             v=v_new             
             return
             end

              subroutine tan_map(x,v,x_new,v_new,eps,tau,DM)
              Implicit real*8(A-H,O-Z)
              real*8 DM(2,2)
c------------------
c              Tangent reflection map: see formulas on the paper
c              
c              input  x=x_{n-1}, v=v_{n-1}  and x_new=x_n  v_new= v_n
c              
c              output DM(x_{n-1}, v_{n-1})              
c------------------                           
              ffp=fp(x+tau*v)
              ffpp=fpp(x+tau*v)
              
              vz=sqrt(1-v**2)
              vz_new=sqrt(1-v_new**2)
              
              den_d= vz-eps*v*ffp
              d_tau_x= eps*ffp/den_d
              d_tau_v= (eps*tau*ffp+tau*v/vz)/den_d
              

              den= 1+(eps*ffp)**2
              coef_x= ( -2*eps*v*ffp + vz*(1-(eps*ffp)**2) )/den**2
              coef_x=coef_x * 2*eps*ffpp
              
c----------   tangent map  at  (x=x_{n-1},v=v_{n-1})
              
              
              DM(2,1)= coef_x*(1+v*d_tau_x)

              DM(2,2)= 1- (eps*ffp + v/vz)* 2*eps*ffp/den +
     *                     coef_x*(tau+v* d_tau_v)

              DM(1,1)= 1 + (v+v_new* vz/vz_new)*d_tau_x +
     *             tau*vz/vz_new**3 * DM(2,1)

              DM(1,2)= tau + (v+ v_new* vz/vz_new)*d_tau_v + 
     *             tau* vz/vz_new**3 * DM(2,2)- tau*v*v_new/(vz*vz_new)
              

c----------   Tangent map in the coordinates  (xx,v)   xx=x/(2*pi)
              
c                   duepi=8*datan(1.d0)
c                   DM(1,1)=  DM(1,1)
c                   DM(1,2)=  DM(1,2)/duepi
c                   DM(2,1)=  DM(2,1)*duepi
c                   DM(2,2)=  DM(2,2)
c
c                   Notice that ur coordinates are x, v we use xx= x/(2 pi) mod 1             
c                   only to plot the results 
c-------------              
              return
              end

               real*8 function f(x)
               Implicit real*8(A-H,O-Z)
c--------------
c            Corrugation function
c--------------               
               f=cos(x)
               return
               end

               real*8 function fp(x)
               Implicit real*8(A-H,O-Z)
c--------------
c            First derivative of orrugation function
c--------------                  
               fp=-sin(x)
               return
               end


               real*8 function fpp(x)
               Implicit real*8(A-H,O-Z)
               fpp=-cos(x)
c--------------
c             Second  derivative of orrugation function
c--------------                    
               return
               end

              subroutine conv(x,v,xx,y)
              Implicit real*8(A-H,O-Z)
              duepi=8*datan(1.d0)
c--------------
c             Conversion from x,v  to xx= x/(2 pi) mod 1   and  y= arcos(v)/(2 pi) 
c--------------                
              theta=acos(v)
              y=theta/duepi                        
              xx=x/duepi
              ix=xx
              xx=xx-ix
              if(xx.le.0) xx=xx+1
              if(xx.gt.0.5) xx=xx-1
              return
              end 


              Subroutine mult(A,B)
                Implicit real*8(A-H,O-Z)
                Real*8 A(2,2), B(2,2), C(2,2)
c------------
c                B = A*B
c------------
                do i=1,2
                do j=1,2
                   C(i,j)=0
                   do k=1,2
                      C(i,j)=C(i,j) + A(i,k)*b(k,j)
                   enddo
                enddo
                enddo

                do i=1,2
                   do j=1,2
                      B(i,j)= C(i,j)
                   enddo
                enddo

                return
                end


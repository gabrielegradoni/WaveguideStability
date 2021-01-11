             Implicit Real*8(A-H,O-Z)
             Real*8  DM(2,2),  L_n(2,2)
             Real*8 LE, LE_old
c------------------------------
c      Choice: 1     phase portrait  with initial conditions on a 10x10 grid
c      Choice: 2     errors LE, RE, REM and det on a grid    N_gridxN_grid    
c------------------------------       
                  Ich=1             ! orbits                   output fort.20
c                   Ich=2             ! LE color plot            output fort.30
c-------------------------------
c     Input parameters: iteration number, grid points number,   corrugation amplitude
c-------------------------------              
              
              N_max= 200           ! iterations number   
              N_grid=100            ! grid N_gridxN_grid
              
              if(Ich.eq.1) then    ! suggested values 
              N_max=2000
              N_grid=10
              endif

              eps=0.1              ! corrugation amplitude                    
              
c--------------
c     Compute the orbit  (n,v_n)= M(x_{n-1},v_{n-1})
c              
c     Compute the tangent map   by          L_n = DM_{n-1} *  L_{n-1}         L_0=I
c
c     Compute the  Lyapunov invariants are   I_1= Tr(L_n^T L_n)   and   I_2= det(L_n^T L_n)
c
c     Lyapunov error LE_n= sqrt(I_1)    and  det_n =sqrt(I_2) = det(L_n)        Notice that det should be 1             c 
c    
c     Reversibility error   RE defined by the recurrence 
c              
c     (RE_n)**2=  (RE_n-1)**2  +   (LE_n)**2 +  (LE_n-1)**2      (LE_0)**2=2=   RE_0=0
c
c     Two ways to compute det_n:
c
c     det_n=det(L_n)    round off error grows fast 
c
c     det_mod_n=  det_mod_{n-1} * det(DM_{n-1} )  with   det_mod_0   round off does grow excet neat v=1 or v=-1       
c--------------
              

             pi=4*datan(1.d0)
             duepi=2*pi        
             A_grid=N_grid                      
 
             do Ix=1,N_grid                                            ! iteration for x_0  on the grid
                x0= pi*(-1+  2*Ix/A_grid)      !   -pi < x < pi     
                write(*,*)  Ix
             do Iv=1, N_grid-1                                           ! iteration for v_0  on the grid
                   v0= -1+ 2* Iv/A_grid        ! -1 < v <1
                   call conv(x0,v0,xx0,y0)
                   x=x0
                   v=v0
                   
                   xx0=x0/duepi
                   write(15,'(I5,2g13.5)')   Ix,xx0,v0             
                   write(20,*) '  '
                   write(20,'(2I5,2g13.5)')  Ix,Iv, xx0,v0
                   write(20,*) '  '


             L_n(1,1)=1
             L_n(1,2)=0
             L_n(2,1)=0    ! Initialization tangent map
             L_n(2,2)=1
             
             det_n_mod=1   ! initialization of det(L_n)   

             LE_old=sqrt(2.d0)    ! nitialization for the recurrence of RE
             RE2=0                ! Initialization RE
                   
             do n=1,N_max       ! ----->   Map iteration loop
                
                      x_old=x                     ! x=x_old=x_{n-1}
                      v_old=v                     ! v=v_old=v(n-1}
                      
             call mappa(x,v,eps,tau,zer) !  orbit  recurrence x_{n-1}, v_{n-1}--->  x_n, v_n
                      
             call conv(x,v,xx,y)         !  xx=x/2 pi mod 1   y= acos(v)/duepi
c-----
             If(Ich.eq.1) write(20,'(I5,12g13.5)') n, xx, v, y, det, Tr    ! Phase portrait 
c-----                
             call tan_map(x_old,v_old,x,v,eps,tau,DM) ! DM  = DM_{n-1}= DM(x_old,v_old)
             
             call mult(DM,L_n)  ! L_n =DM_{n-1} * L_{n-1}
             
                det=DM(2,2)*DM(1,1)-DM(1,2)*DM(2,1)   ! det = det(DM(x_old,v_old))= det(L
                det_n_mod= det_n_mod*det              ! det_mod_n=  det_mod_{n-1} * det(DM_{n-1} )      
        det_n=L_n(2,2)*L_n(1,1)-L_n(1,2)*L_n(2,1)     ! det_n = det(L_n)     det_n**2 = det (L_n^T L_n)        
        
                       
        Tr_n=L_n(1,1)**2 + L_n(1,2)**2 + L_n(2,1)**2 + L_n(2,2)**2   !   Tr_n= Tr(L_n^T L_n)  first invariant I_1
c
                   LE= sqrt(Tr_n)                  ! error   LE_n 
                   RE2=RE2 + LE_old**2 + LE**2     ! 
                   RE=sqrt(RE2)                    ! error   RE_n
                   
                   LE_old=LE
                   enddo                  ! <-----   end iterations loop 
                    

c-------------     Error REM 
                v=-v                            ! velocity reversed 
                do n=1,N_max
                call mappa(x,v,eps,tau,zer)     ! backward iteration
                enddo
                v=-v
                call conv (x,v,xx,y)
                REM=sqrt((xx-xx0)**2+(v-v0)**2)/1E-16 !   error REM  for  n=N_max
c-------------

c                Restriction to be in the range of gnuplot: can be removed if   8 bytes reals
c                have eponent 10**(m) with -300 <m< 300. Gnuplot same exponent as Real*4  -30<m<30 
                
                 If(LE.ge.1.E30) LE=1.E30    ! For gnuplot
                 If(RE.ge.1.E30) RE=1.E30    ! For gnuplot
                 det_n=abs(det_n)
                 if(det_n.eq.0)    det_n=1.E-30 ! For gnuplot
                 if(det_n.ge.1E30) det_n=1.E30 ! For gnuplot
                 det_n_mod=abs(det_n_mod)
                 if(det_n_mod.eq.0)    det_n_mod=1.E-30  ! For gnuplot
                 if(det_n_mod.ge.1E30) det_n_mod=1.E30  ! For gnuplot
                 
             write(30,'(2I5,3x,8g13.5)') Ix, Iv, xx0, v0, y0, LE,
     *              RE, REM, det_n, det_n_mod 

c             Ch_LE= Ch_LE+ dlog(LE)/(N_max*N_grid**2)    ! Channel Capacity 
c             Ch_RE= Ch_RE+ dlog(RE)/(N_max*N_grid**2)     ! Channel Capacity 
c             Ch_REM=Ch_RE+ dlog(REM)/(N_max*N_grid**2)    ! Channel Capacity    
 

             enddo     ! end of the x_0 loop
             enddo     ! end of the v_0 loop 


c             write(*,'(5g13.5)')   eps,Ch_LE,Ch_RE,Ch_REM                         

             
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


          Implicit real*8(A-H,O-Z)
          real*8 nux,nuy,nuz  
          common/bb/kx,kv,n     ! rinuovere ???
          common/zero/zer
c---------
c         Output file of indicators on the grid:          fort.30
c          
c         Inacurrate values for bisection  written in:    fort.10 
c          
c---------
          duepi=8*datan(1.d0)
          pi=duepi/2
c----------------------------------------------
c          Corrugation function f(x,y) defined within  the function  g()
c
c          Choice of parameters
c----------------------------------------------          
         
          N_max= 200    ! Number of iterations 
          N_grid=100    ! Grid points  N_gridxN_grid 
          
          eps=0.1       ! corrugation amplitude 
          eta=1E-14     ! deviation amplitude for  shadow orbit to compute tangent map
                        ! derivatives replaced by finite diffrences
  
          y0=pi/4.d0    !  value of y0     kept fixed 
          phi0=0        !  value of  phi0  in [0,pi/2]  kept fixed 

                        !  x0, v0  vary on the grid   x0 in [-pi,pi]   v0  in [-1,1]

          A_grid=N_grid

           write(*,*)  '  Fixed values: y0, phi0 ', y0,phi0
           write(10,*) '  Inaccuracy  in bisection method'
           write(10,*) '  kx   kv    n   Iter      diff          zer'
c------------------------
c     Starts the conputation of the indicators LE, RE, REM on the grid
c------------------------          
                  
             do kx= 1,N_grid  ! 1, 10         ---------------------->   Loop on the grid points
              write(*,*) '----', kx
             do kv= 0,N_grid  !  0,10         --------------->  Loop on the grid points
         
                x0= pi*(-1+  2*kx/A_grid)
                v0=-0.98 + 2*0.98 *kv/A_grid    !  v0=sqrt(vx0**2+vy0**2)

                vx0=v0*cos(phi0)                !  vx0
                vy0=v0*sin(phi0)                !  vy0     
                
                call conv(x0,y0,xx0,yy0)
c---------------                
c         xx0= x0/(2 pi) mod 1      yy0= y0/(2 pi) mod 1
c---------------
              x=x0
              vx=vx0
              y=y0
              vy=vy0

c---------------------------------------------------------------------
c         Computation of tangent  L_n(x) = DM^n(x) with finite difference approximation
c--------
c         Letting M_1=M_x  M_2=M_v_x   M_3= m_y   M_4=M_v_y    x_1=x  x_2=v_x  x_3=y  x_4=v_y
c              
c         we approximate  (DM^n(x))_{i,j}(x_1,x_2,x_3,x_4) = d(M^n)_i)/dx_j   
c
c        (DM^n(x))_(i,1)= [ DM^n_i(x_1+eta,x_2,x_3,x_4)-DM^n_i(x_1,x_2,x_3,x_4) ]/eta
c        (DM^n(x))_(i,2)= [ DM^n_i(x_1,x_2+eta,x_3,x_4)-DM^n_i(x_1,x_2,x_3,x_4) ]/eta     etc. 
c----------------------------------------------------------------------              
             x_L1 =x0 + eta   ! initial displacemtent along x axis to compute   d M_i/dx 
             vx_L1=vx0         
             y_L1 =y0
             vy_L1=vy0

             x_L2 =x0         ! initial displacemtent along vx axis to compute   d M_i/dv_x 
             vx_L2=vx0 +eta
             y_L2 =y0
             vy_L2=vy0

             x_L3 =x0         ! initial displacemtent along y axis to compute   d M_i/dy 
             vx_L3=vx0
             y_L3 =y0 +eta
             vy_L3=vy0

             x_L4 =x0         ! initial displacemtent along y axis to compute   d M_i/dv_y 
             vx_L4=vx0
             y_L4 =y0
             vy_L4=vy0 +eta
                                  
             
             er2_RE=0        ! initialization of the squares of LE and LE
             er_LE_old=2

c------------
             
              do n=1, N_max     ! ---->  Iterations of the map
                
                 call mappa(x,y,vx,vy,tau, eps)
                 
                 call mappa(x_L1, y_L1, vx_L1, vy_L1,  tau, eps)
                 call mappa(x_L2, y_L2, vx_L2, vy_L2,  tau, eps)
                 call mappa(x_L3, y_L3, vx_L3, vy_L3,  tau, eps)
                 call mappa(x_L4, y_L4, vx_L4, vy_L4,  tau, eps)

           er2_1= (x_L1-x)**2+(y_L1-y)**2+ (vx_L1-vx)**2+ (vy_L1-vy)**2
           er2_2= (x_L2-x)**2+(y_L2-y)**2+ (vx_L2-vx)**2+ (vy_L2-vy)**2
           er2_3= (x_L3-x)**2+(y_L3-y)**2+ (vx_L3-vx)**2+ (vy_L3-vy)**2
           er2_4= (x_L4-x)**2+(y_L4-y)**2+ (vx_L4-vx)**2+ (vy_L4-vy)**2


c----------------------------------------------------------------------------
c          Tangent map elements  L_n =  DM^n  are given by
c           
c          (L_n)_11 = (x_L1-x)/\eta    (L_n)_21 = (vx_L1-vx)/eta
c          (L_n)_31 = (y_L1-y)/\eta    (L_n)_41 = (vy_L1-vy)/eta
c           
c          (L_n)_12 = (x_L2-x)/\eta    (L_n)_22 = (vx_L2-vx)/eta
c          (L_n)_32 = (y_L2-y)/\eta    (L_n)_42 = (vy_L2-vy)/eta
c
c          (L_n)_13 = (x_L3-x)/\eta    (L_n)_23 = (vx_L3-vx)/eta
c          (L_n)_33 = (y_L3-y)/\eta    (L_n)_43 = (vy_L3-vy)/eta
c           
c          (L_n)_14 = (x_L2-x)/\eta    (L_n)_24 = (vx_L2-vx)/eta
c          (L_n)_34 = (y_L2-y)/\eta    (L_n)_44 = (vy_L2-vy)/eta     
c-----------------------------------------------------------------------------
c           LE^2 = Tr( (L_n)^T L_n) =  sum_{ij}  L_n)_ij * L_n)_ij
c-----------------------------------------------------------------------------           
c
c-----------------------------           
c          Computation of  LE
c------------------------------           
           
             er2_LE= er2_1   + er2_2  + er2_3 + er2_4      ! er2_LE=  Tr( (L_n)^T L_n) *eta**2 

            
             er_LE= sqrt(er2_LE)/eta                       ! er_LE= sqrt[ Tr( (L_n)^T L_n) ]

c---------------------------
c      Computation of REM
c--------------------------             
                        
            er2_RE= er2_RE + er_LE**2 + er_LE_old**2      ! recurrence 
            er_RE= sqrt(er2_RE)
            
            er_LE_old= er_LE           

         enddo                  ! <-------- End of the iterations of the map
         
c--------------------
                                 
c------------------------------              
c    Reversibility error due to round off   REM           
c-----------------------------             
               vx=-vx                 !   inversion of velocities 
               vy=-vy                 !   inversion of velocities
               
                do n=1, N_max         ! backawrd interations of the map
                call mappa(x,y, vx,vy,tau, eps)              
                enddo                 ! emd of backward itereations 
                vx=-vx
                vy=-vy

             er2= (x-x0)**2+(y-y0)**2+ (vx-vx0)**2+ (vy-vy0)**2         
             er_REM= sqrt(er2)/1E-17
             
c---        errore calcolato con conversione da Conv
             
             call conv(x,y,xx,yy)
             err2= (xx-xx0)**2+(yy-yy0)**2+ (vx-vx0)**2+ (vy-vy0)**2         
             err_REM= sqrt(err2)/1E-17
c----------             
c           Output values 
c----------             
             write(30,'(2I5,10g15.5)') kx,kv,x0/duepi, v0, 
     *       y0/duepi, er_LE, er_RE, er_REM, zer        ! zer accuracy of bisection to compute tau 
                                          
          enddo    ! <-------------------    End of loop on grid points
          enddo    ! <-------------------------- End of loop on grid points
          end

   

          

      

          real*8 function g( x, y, vx, vy, tau, eps)
          Implicit real*8(A-H,O-Z)
c--------
c        To determine tau we have to sove   g=0
c--------          
          argx=x+tau*vx
          argy=y+tau*vy
          vz=sqrt(1-vx**2-vy**2)
c--------
c         corrugation function f(x,y) In the present case    f(x,y)= cos(x)cos(y)
c--------          
          f=cos(argx)*cos(argy)            
          g=tau*vz-1-eps*f       
          return
          end

c          real*8 function gp(x, y, vx,vy,tau, eps)
c          Implicit real*8(A-H,O-Z)
c          pi=4*datan(1.d0)          
c          argx=x+tau*vx
c          argy=y+tau*vy
c          fp=-sin(argx)*cos(argy)*vx -sin(argy)*cos(argx)*vy
c          vz=sqrt(1-vx**2-vy**2)
c          gp=vz - eps*fp
c          return
c          end
      

             subroutine mappa(x,y,vx,vy,tau,eps)
             Implicit real*8(A-H,O-Z)
             real*8 nux,nuy, nuz
             common/bb/kx,kv,n
             common/zero/zer
c--------------------------             
c            Reflections map: se formula (14) of the paper
c--------------------------
             
             pi=4*datan(1.d0)
             duepi=2*pi
             vz=sqrt(1-vx**2-vy**2)              !   vz= sin(theta) 
         

              call bisec(x,y,vx,vy,taub, eps)    ! computation of tau by bisection 
              tau=taub
              zer=g(x, y, vx, vy, tau,  eps)

           
              zer=g( x, y, vx, vy, tau,  eps)
              
              If(abs(zer).ge.1E-12) then
                 call conv(x,y,xx,yy)
                 write(20,'(3I6,3g12.4)') kx,kv,n,xx,vx,zer
              endif   
              
              argx=x+tau*vx
              argy=y+tau*vy

              fx=-cos(argy)*sin(argx)
              fy=-cos(argx)*sin(argy)

             sq=sqrt(1+ eps**2*(fx**2+fy**2) )
             nux=eps*fx/sq
             nuy=eps*fy/sq
             nuz=-1/sq

             vx_old= vx
             vy_old= vy
             vz_old= vz
            
             scal= vx_old*nux+ vy_old*nuy+ vz_old*nuz
             
             vx=vx_old-2*nux*scal
             vy=vy_old-2*nuy*scal
             vz=sqrt(1-vx**2-vy**2) 
             
             vz_coll=vz_old-2*nuz*scal ! Velocity  after reflection of  corrugared  plane  v_z* 
                                       ! Velocity  after reflection on  lower palne vz_coll= -vz*   

             x=x+ tau*( vx_old + vx* vz_old/vz)
             y=y+ tau*( vy_old + vy* vz_old/vz)
             
             return            
             end 

            subroutine conv(x,y,xx,yy)
            Implicit real*8(A-H,O-Z)
            duepi=8*datan(1.d0)
c------------------------------------------------
c      Conversion  xx= x/(2 pi)  mod 1   and   yy= y/(2 pi)  mod 1 
c------------------------------------------------
             xx=x/duepi
             ix=xx
             xx=xx-ix
             if(xx.le.0) xx=xx+1
             if(xx.gt.0.5) xx=xx-1

             yy=y/duepi
             iy=yy
             yy=yy-iy
             if(yy.le.0)   yy=yy+1
             if(yy.gt.0.5) yy=yy-1
             
             return            
             end 

             real*8 function phase(x,y)
             Implicit real*8(A-H,O-Z)
c------------
c      Phase in [0,2 pi] of a vector   (x,y) 
c------------             
             pi=4*datan(1.d0)
             duepi=2*pi
             phase =acos(x)
             if(y.lt.0)  phase=duepi-phase
             return
             end 


             

          subroutine bisec( x,y,vx,vy,tau, eps)
          implicit real*8(A-H,O-Z)
          common/bb/kx,kv,n   
c------------------        
c        Bisection method to compute tau by solving g=g(x,y,v_x,v_y, tau)=0
c-----------------          
             Precis=1.E-15             !  precision required by the bisection metod 
             vz=sqrt(1-vx**2-vy**2)
             tau_a=0.5d0/vz
             tau_b=1.5d0/vz
             G_a= g(x, y, vx, vy, tau_a,  eps) 
             G_b= g(x, y, vx, vy, tau_b,  eps)
             
             if(G_a*G_b.gt.0) write(*,*) ' ***** no zero ****'   ! no convergence

               do iter=1,60
                tau_c=(tau_a+tau_b)/2
                G_c= g(x, y, vx, vy, tau_c,  eps)
                if(G_a*G_c.lt.0) then
                   tau_b=tau_c
                   G_b=G_c
                else if(G_b*G_c.lt.0) then
                   tau_a=tau_c
                   G_a=G_c
                endif
                diff=(tau_b-tau_a)/(tau_b+tau_a)
              
                tau=(tau_a+tau_b)/2
                zer=g(x, y, vx, vy, tau,  eps)
               
                if(abs(diff).le.precis) go to 1   ! precisione raggiunta 
             enddo
 1         tau=(tau_a+tau_b)/2
           zer=g(x, y, vx, vy, tau,  eps)
c------------
           if(diff.gt.1E-12.or.abs(zer).ge.1.E-12)  then
           write(10,'(4I5,2g13.4)')  kx,kv,n,Iter,diff,zer
           endif
           
             return
           end

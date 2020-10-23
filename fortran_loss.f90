!------------------------------------------------------------------------
subroutine loop(Lz,Lr,nr,nmax,W, Rr, Zz, r, z, nx, nz, tz, T, q, f, phi, c0, angle,m,amp,C,P)

  implicit none
  
  integer, intent(in):: Lz,Lr,nr,nmax
  real*8, intent(in):: f, c0
  real*8, dimension(Lr, Lz), intent(in):: Rr, Zz
  real*8, dimension(nmax,nr), intent(in):: W, r, z, nx, nz, tz, T, q, phi, angle, m, amp, C
  complex*8, dimension(Lr, Lz), intent(out):: P

  integer:: k,l, zl, rl
  real, parameter:: pi = 3.14159265
  integer, parameter:: out_unit = 20
  complex:: i
  real*8:: A_, s_, n_, al, T_, q_, W_, r_, n_val

  i = (0,1)

  open (unit=out_unit,file="results.txt",action="write",status="replace")

  !number of rays
  do l = 1,nr
     !number of points per ray
     
     do k = 2,nmax

        if (W(k,l).ne.0) then
          !columns of R
           do rl = 1,Lr
              
              if ( (Rr(rl,1).gt.min(r(k,l),r(k-1,l))).and.(Rr(rl,1).lt.max(r(k,l),r(k-1,l)))) then
                
                !columns of Z
                 do zl = 1,Lz

                   n_ = abs((Rr(rl,zl) - r(k,l))*nx(k,l) + (Zz(rl,zl) - z(k,l))*nz(k,l) )
                   s_ = (Rr(rl,zl) - r(k-1,l))*nx(k-1,l) + (Zz(rl,zl) - z(k-1,l))*tz(k-1,l)
                   al = s_ / sqrt( (r(k-1,l) - r(k,l))**2 + (z(k-1,l) - z(k,l))**2 )

                   T_ = T(k-1,l) + al * (T(k,l) - T(k-1,l))
                   q_ = q(k-1,l) + al * (q(k,l) - q(k-1,l))
                   W_ = W(k-1,l) + al * (W(k,l) - W(k-1,l))
                   r_ = r(k-1,l) + al * (r(k,l) - r(k-1,l))

                   A_ = 1 / (4*pi) * amp(k,l) * i * m(k,l) * sqrt(abs( C(k,l) * cos( angle(k,l) / (r_ * c0 * q_))))

                   n_val = (W_ - n_) / W_
                   write (out_unit,*) n_val
                   if( n_val.gt.0) then
                      P(rl,zl) = P(rl,zl) + n_val * A_ * exp(i*2*pi*f*T_) * exp(i*phi(k,l))
                      
                   endif
                enddo
             endif
          enddo

        endif

     enddo
  enddo

  close (out_unit)

end subroutine loop
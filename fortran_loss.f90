!------------------------------------------------------------------------
subroutine loop(Lz,Lr,nr,nmax,W, Rr, Zz, r, z, nx, nz, tx, tz, T, q, f, phi, c0, angle,m,amp,C,P,n_bot)

  implicit none

  integer, intent(in):: Lz,Lr,nr,nmax
  real*8, intent(in):: f, c0
  real*8, dimension(Lr, Lz), intent(in):: Rr, Zz
  real*8, dimension(nmax,nr), intent(in):: W, r, z, nx, nz, tx, tz, T, q, phi, angle, m, amp, C,n_bot
  complex*16, dimension(Lr, Lz), intent(out):: P

  integer:: k,l, zl, rl
  real*8, parameter:: pi = 3.14159265
  integer, parameter:: out_unit = 20
  complex*16:: i
  complex*16::A_
  real*8::s_, n_, al, T_, q_, W_, r_, n_val

  i = (0,1)

  !open (unit=out_unit,file="results.txt",action="write",status="replace")

  !number of rays
  do l = 1,nr
    !l = 17
     !number of points per ray

     do k = 2,nmax-1

        if (W(k,l).ne.0) then
        !write (out_unit,*) W(k,l)
        !columns of R
         do rl = 1,Lr-1
            !write (out_unit,*) Rr(1,rl)
            if ((Rr(1,rl).ge.min(r(k,l),r(k-1,l))).and.(Rr(1,rl).lt.max(r(k,l),r(k-1,l))).and.(n_bot(k,l).lt.6)) then

              !columns of Z
               do zl = 1,Lz

                 n_ = abs((Rr(zl,rl) - r(k,l))*nx(k,l) + (Zz(zl,rl) - z(k,l))*nz(k,l) )
                 s_ = (Rr(zl,rl) - r(k-1,l))*tx(k-1,l) + (Zz(zl,rl) - z(k-1,l))*tz(k-1,l)
                 al = s_ / sqrt( (r(k-1,l) - r(k,l))**2 + (z(k-1,l) - z(k,l))**2 )

                 T_ = T(k-1,l) + al * (T(k,l) - T(k-1,l))
                 q_ = q(k-1,l) + al * (q(k,l) - q(k-1,l))
                 W_ = W(k-1,l) + al * (W(k,l) - W(k-1,l))
                 r_ = r(k-1,l) + al * (r(k,l) - r(k-1,l))

                 A_ = 1 / (4*pi) * amp(k,l) * (-i)**m(k,l) * sqrt(abs( C(k,l) * cos( angle(k,l)) / (r_ * c0 * q_)))
                 !A_ = (-i)**m(k,l)
                 n_val = (W_ - n_) / W_
                 !write (out_unit,*)  al

                 if (n_.le.W_) then
                    P(zl,rl) = P(zl,rl) + n_val * A_ * exp(i*2*pi*f*T_) * exp(i*phi(k,l))
                 endif
              enddo
             endif
          enddo

        endif

     enddo
  enddo

  close (out_unit)

end subroutine loop

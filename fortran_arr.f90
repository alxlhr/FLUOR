!------------------------------------------------------------------------
subroutine loop(nr,nmax,W, r_rcvr, z_rcvr, r, z, nx, nz, tx, tz, T, q, f, phi, c0, &
  angle,m,amp,C,Amp_rcvr,ray_num,Angle_rcvr,Delay_rcvr,eigen_ray)

  implicit none

  integer, intent(in):: nr,nmax
  real*8, intent(in):: f, c0, r_rcvr, z_rcvr
  real*8, dimension(nmax,nr), intent(in):: W, r, z, nx, nz, tx, tz, T, q, phi, angle, m, amp, C
  real*8, dimension(nr), intent(inout)::Amp_rcvr,Angle_rcvr,Delay_rcvr
  integer, intent(out):: eigen_ray
  integer, dimension(nr), intent(inout)::ray_num
  !complex*16, dimension(Lr, Lz), intent(out):: P

  integer:: k,l, zl, rl
  real*8, parameter:: pi = 3.14159265
  integer, parameter:: out_unit = 20
  complex*16:: i
  complex*16:: A_
  real*8::s_, n_, al, T_, q_, W_, r_, angle_

  i = (0,1)
  eigen_ray = 0 !remove 1 in the python code (fortran indexing starts in 1)
  !open (unit=out_unit,file="results_arr.txt",action="write",status="replace")
  !number of rays
  do l = 1,nr
    !number of points per ray

     do k = 2,nmax-1

        if (W(k,l).ne.0) then

            if ((r_rcvr.ge.min(r(k,l),r(k-1,l))).and.(r_rcvr.lt.max(r(k,l),r(k-1,l)))) then

                 n_ = abs((r_rcvr - r(k,l))*nx(k,l) + (z_rcvr - z(k,l))*nz(k,l) )

                 if (n_.le.2.0) then

                   s_ = (r_rcvr - r(k-1,l))*tx(k-1,l) + (z_rcvr - z(k-1,l))*tz(k-1,l)
                   al = s_ / sqrt( (r(k-1,l) - r(k,l))**2 + (z(k-1,l) - z(k,l))**2 )

                   T_ = T(k-1,l) + al * (T(k,l) - T(k-1,l))
                   q_ = q(k-1,l) + al * (q(k,l) - q(k-1,l))
                   W_ = W(k-1,l) + al * (W(k,l) - W(k-1,l))
                   r_ = r(k-1,l) + al * (r(k,l) - r(k-1,l))
                   angle_ = angle(k-1,l) + al * (angle(k,l) - angle(k-1,l))

                   A_ = 1 / (4*pi) * amp(k,l) * (-i)**m(k,l) * sqrt(abs( C(k,l) * cos( angle(k,l)) / (r(k,l) * c0 * q_)))
                   !A_ = (-i)**m(k,l)
                   eigen_ray = eigen_ray + 1
                   Angle_rcvr(eigen_ray) = angle_
                   Delay_rcvr(eigen_ray) = T_

                   Amp_rcvr(eigen_ray) = abs(A_*exp(i*phi(k,l)))

                   ray_num(eigen_ray) = l - 1 !to use in python indexing

               endif
             endif


        endif

     enddo
  enddo
  eigen_ray = eigen_ray - 1 !to use in python indexing
  !close (out_unit)
end subroutine loop

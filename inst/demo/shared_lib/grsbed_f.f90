subroutine mtgrsbed_core_f( &
    file, n, cls, af, m, nt, scale, &
    b_snp, grs_flat, nthreads, MG, JB, TB) bind(C, name="mtgrsbed_core_f")

  use iso_c_binding
  use omp_lib
  implicit none

  character(kind=c_char), dimension(*) :: file
  integer(c_int), value :: n, m, nt, nthreads, MG, JB, TB
  integer(c_int), dimension(*) :: cls
  real(c_double), dimension(*) :: af
  real(c_double), dimension(*) :: b_snp
  real(c_double), dimension(*) :: grs_flat
  logical(c_bool), value :: scale

  integer :: unit, ios, i_char
  integer(c_int) :: i0, imax, mlen, ii, i
  integer(c_int) :: j0, jmax, j, t0, tmax, t
  integer(c_int) :: byte0, byte1, kb, pos, jbase
  integer(c_int) :: tlen
  integer(c_int) :: nbytes
  integer(c_long_long) :: offset
  integer(c_int) :: gcode, ioff, joff

  integer(c_int8_t), allocatable :: block_buffer(:)
  real(c_double), allocatable :: map0(:), map1(:), map2(:), map3(:)

  real(c_double) :: p, denom, xj
  integer(c_int8_t) :: buf_k
  integer(c_int) :: error
  character(len=1024) :: filename

  error = 0
  filename = ""

  do i_char = 1, 1024
    if (file(i_char) == c_null_char) exit
    filename(i_char:i_char) = file(i_char)
  end do

  nbytes = (n + 3) / 4

  allocate(block_buffer(MG*nbytes))
  allocate(map0(MG), map1(MG), map2(MG), map3(MG))

  open(newunit=unit, file=trim(filename), access="stream", form="unformatted", &
       status="old", action="read", iostat=ios)

  if (ios /= 0) then
    stop "Could not open BED file"
  end if

  if (nthreads > 0) then
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(nthreads)
  end if

!$omp parallel private(i0, imax, mlen, ii, i, p, denom, &
!$omp& j0, jmax, byte0, byte1, kb, pos, jbase, t0, tmax, &
!$omp& tlen, t, xj, buf_k, offset, gcode, ioff, joff)

  do i0 = 1, m, MG
    imax = min(i0 + MG - 1, m)
    mlen = imax - i0 + 1

!$omp single
    do ii = 1, mlen
      if (error /= 0) cycle

      i = i0 + ii - 1
      offset = int((cls(i)-1)*nbytes + 3, kind=c_long_long)

      read(unit, pos=offset+1, iostat=ios) block_buffer((ii-1)*nbytes+1 : ii*nbytes)

      if (ios /= 0) then
!$omp atomic write
        error = 1
        cycle
      end if

      p = af(i)

      if (scale) then
        denom = sqrt(2.0d0*p*(1.0d0-p))
        if (denom <= 0.0d0) then
          map0(ii)=0.0d0
          map1(ii)=0.0d0
          map2(ii)=0.0d0
          map3(ii)=0.0d0
        else
          map0(ii)=(2.0d0-2.0d0*p)/denom
          map1(ii)=0.0d0
          map2(ii)=(1.0d0-2.0d0*p)/denom
          map3(ii)=(-2.0d0*p)/denom
        end if
      else
        map0(ii)=2.0d0
        map1(ii)=-2.0d0*p
        map2(ii)=1.0d0
        map3(ii)=0.0d0
      end if
    end do
!$omp end single

!$omp do schedule(static)
    do j0 = 1, n, JB
      if (error /= 0) cycle

      jmax = min(j0 + JB - 1, n)
      byte0 = (j0-1) / 4 + 1
      byte1 = (jmax + 3) / 4

      do t0 = 1, nt, TB
        tmax = min(t0 + TB - 1, nt)
        tlen = tmax - t0 + 1

        do ii = 1, mlen
          i = i0 + ii - 1
          ioff = (i-1)*nt + t0

          do kb = byte0, byte1
            buf_k = block_buffer((ii-1)*nbytes + kb)
            jbase = (kb-1)*4 + 1

            do pos = 0, 3
              j = jbase + pos

              if (j < j0 .or. j > jmax .or. j > n) then
                buf_k = ishft(buf_k, -2)
                cycle
              end if

              gcode = iand(int(buf_k, kind=c_int), 3_c_int)

              select case (gcode)
              case (0)
                xj = map0(ii)
              case (1)
                xj = map1(ii)
              case (2)
                xj = map2(ii)
              case default
                xj = map3(ii)
              end select

              buf_k = ishft(buf_k, -2)
              joff = (j-1)*nt + t0

!$omp simd
              do t = 0, tlen - 1
                grs_flat(joff+t) = grs_flat(joff+t) + b_snp(ioff+t) * xj
              end do

            end do
          end do
        end do
      end do
    end do
!$omp end do

  end do
!$omp end parallel

  close(unit)
  deallocate(block_buffer, map0, map1, map2, map3)

end subroutine mtgrsbed_core_f

program unit_test_solar_geometry_mo

  use, intrinsic :: iso_fortran_env, only: real64
  use solar_geometry_mo

  implicit none

  integer      :: doy_summer, doy_winter, nfail
  real(real64) :: elev, az, i0
  real(real64), parameter :: TOKYO_LAT = 35.69_real64, TOKYO_LON = 139.69_real64

  nfail = 0
  print *, '=== fortran-solar : solar_geometry_mo unit tests ==='

  ! ---- day_of_year (閏年対応) ----
  call check_i('day_of_year 2026-01-01',           day_of_year(2026, 1, 1),  1,   nfail)
  call check_i('day_of_year 2026-12-31 (non-leap)', day_of_year(2026,12,31),  365, nfail)
  call check_i('day_of_year 2024-12-31 (leap)',     day_of_year(2024,12,31),  366, nfail)
  call check_i('day_of_year 2024-03-01 (leap)',     day_of_year(2024, 3, 1),  61,  nfail)
  call check_i('day_of_year 2026-03-01 (non-leap)', day_of_year(2026, 3, 1),  60,  nfail)

  doy_summer = day_of_year(2026, 6, 21)   ! 夏至
  doy_winter = day_of_year(2026,12, 22)   ! 冬至

  ! ---- 太陽高度角 (Tokyo) ----
  elev = solar_elevation_deg(TOKYO_LAT, TOKYO_LON, doy_summer, 12.0_real64)
  call check_range('summer noon elevation high', elev, 70.0_real64, 90.0_real64, nfail)
  elev = solar_elevation_deg(TOKYO_LAT, TOKYO_LON, doy_winter, 12.0_real64)
  call check_range('winter noon elevation low',  elev, 25.0_real64, 40.0_real64, nfail)
  elev = solar_elevation_deg(TOKYO_LAT, TOKYO_LON, doy_summer, 0.0_real64)
  call check_range('midnight elevation negative', elev, -90.0_real64, 0.0_real64, nfail)

  ! ---- 太陽方位角 (正午は概ね南 ≈ 180°) ----
  az = solar_azimuth_deg(TOKYO_LAT, TOKYO_LON, doy_summer, 12.0_real64)
  call check_range('summer noon azimuth ~ south', az, 150.0_real64, 210.0_real64, nfail)

  ! ---- 大気外水平面日射量 ----
  i0 = extraterrestrial_radiation_wm2(TOKYO_LAT, TOKYO_LON, doy_summer, 12.0_real64)
  call check_range('noon I0 plausible', i0, 800.0_real64, 1400.0_real64, nfail)
  i0 = extraterrestrial_radiation_wm2(TOKYO_LAT, TOKYO_LON, doy_summer, 0.0_real64)
  call check_range('night I0 == 0', i0, -0.001_real64, 0.001_real64, nfail)

  print *, '----------------------------------------'
  if (nfail == 0) then
    print *, 'ALL TESTS PASSED'
  else
    print '(1x,i0,a)', nfail, ' TEST(S) FAILED'
    error stop 1
  end if

contains

  subroutine check_i(name, got, want, nf)
    character(*), intent(in)    :: name
    integer,      intent(in)    :: got, want
    integer,      intent(inout) :: nf
    if (got == want) then
      print '(1x,a,a,a,i0)', 'PASS  ', name, ' = ', got
    else
      print '(1x,a,a,a,i0,a,i0)', 'FAIL  ', name, ' got ', got, ' want ', want
      nf = nf + 1
    end if
  end subroutine

  subroutine check_range(name, got, lo, hi, nf)
    character(*), intent(in)    :: name
    real(real64), intent(in)    :: got, lo, hi
    integer,      intent(inout) :: nf
    if (got >= lo .and. got <= hi) then
      print '(1x,a,a,a,f9.3)', 'PASS  ', name, ' = ', got
    else
      print '(1x,a,a,a,f9.3,a,f9.3,a,f9.3,a)', 'FAIL  ', name, ' = ', got, ' not in [', lo, ',', hi, ']'
      nf = nf + 1
    end if
  end subroutine

end program unit_test_solar_geometry_mo

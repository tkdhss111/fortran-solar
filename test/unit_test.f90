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
  call check_i( 'day_of_year 2026-01-01',           day_of_year( 2026, 1, 1 ),  1,   nfail )
  call check_i( 'day_of_year 2026-12-31 (non-leap)', day_of_year( 2026,12,31 ),  365, nfail )
  call check_i( 'day_of_year 2024-12-31 (leap)',     day_of_year( 2024,12,31 ),  366, nfail )
  call check_i( 'day_of_year 2024-03-01 (leap)',     day_of_year( 2024, 3, 1 ),  61,  nfail )
  call check_i( 'day_of_year 2026-03-01 (non-leap)', day_of_year( 2026, 3, 1 ),  60,  nfail )

  doy_summer = day_of_year( 2026, 6, 21 )   ! 夏至
  doy_winter = day_of_year( 2026,12, 22 )   ! 冬至

  ! ---- 太陽高度角 (Tokyo) ----
  elev = solar_elevation_deg( TOKYO_LAT, TOKYO_LON, doy_summer, 12.0_real64 )
  call check_range( 'summer noon elevation high', elev, 70.0_real64, 90.0_real64, nfail )
  elev = solar_elevation_deg( TOKYO_LAT, TOKYO_LON, doy_winter, 12.0_real64 )
  call check_range( 'winter noon elevation low',  elev, 25.0_real64, 40.0_real64, nfail )
  elev = solar_elevation_deg( TOKYO_LAT, TOKYO_LON, doy_summer, 0.0_real64 )
  call check_range( 'midnight elevation negative', elev, -90.0_real64, 0.0_real64, nfail )

  ! ---- 太陽方位角 (正午は概ね南 ≈ 180°) ----
  az = solar_azimuth_deg( TOKYO_LAT, TOKYO_LON, doy_summer, 12.0_real64 )
  call check_range( 'summer noon azimuth ~ south', az, 150.0_real64, 210.0_real64, nfail )

  ! ---- 大気外水平面日射量 ----
  i0 = extraterrestrial_radiation_wm2( TOKYO_LAT, TOKYO_LON, doy_summer, 12.0_real64 )
  call check_range( 'noon I0 plausible', i0, 800.0_real64, 1400.0_real64, nfail )
  i0 = extraterrestrial_radiation_wm2( TOKYO_LAT, TOKYO_LON, doy_summer, 0.0_real64 )
  call check_range( 'night I0 == 0', i0, -0.001_real64, 0.001_real64, nfail )

  ! ---- 入射角（任意方位・傾斜）----
  ! 春分の太陽南中で、傾斜=緯度の南向き面はほぼ正対（入射角≈0）
  block
    integer      :: doy_eq
    real(real64) :: inc_noon, solar_noon
    doy_eq     = day_of_year( 2026, 3, 20 )
    solar_noon = 12.0_real64 + (135.0_real64 - TOKYO_LON) * 4.0_real64 / 60.0_real64
    inc_noon   = incidence_angle_deg( TOKYO_LAT, TOKYO_LON, doy_eq, solar_noon, TOKYO_LAT, 180.0_real64 )
    call check_range( 'equinox noon incidence (tilt=lat,south) ~0', inc_noon, 0.0_real64, 8.0_real64, nfail )
  end block

  ! ---- 最適傾斜角（任意方位）----
  block
    real(real64) :: t_south, t_naha, t_wakkanai, t_flat, t_se, t_east, t_west
    t_south    = optimal_tilt_deg( TOKYO_LAT, TOKYO_LON )                       ! 真南（既定）
    t_naha     = optimal_tilt_deg( 26.20_real64, 127.69_real64 )               ! 那覇   lat 26.2
    t_wakkanai = optimal_tilt_deg( 45.42_real64, 141.68_real64 )               ! 稚内   lat 45.4
    call check_range( 'optimal tilt Tokyo south (deg)', t_south, 18.0_real64, 36.0_real64, nfail )
    call check_range( 'lower lat (Naha) shallower',  t_naha,     0.0_real64, t_south - 1.0_real64, nfail )
    call check_range( 'higher lat (Wakkanai) steeper', t_wakkanai, t_south + 1.0_real64, 90.0_real64, nfail )
    ! 散乱比率↑ → 浅く（keyword 引数）
    t_flat = optimal_tilt_deg( TOKYO_LAT, TOKYO_LON, diffuse_fraction=0.8_real64 )
    call check_range( 'higher diffuse => flatter', t_flat, 0.0_real64, t_south - 1.0_real64, nfail )
    ! 任意方位: 南東(135)は南より浅い、東(90)と西(270)はほぼ対称
    t_se   = optimal_tilt_deg( TOKYO_LAT, TOKYO_LON, panel_azimuth_deg=135.0_real64 )
    t_east = optimal_tilt_deg( TOKYO_LAT, TOKYO_LON, panel_azimuth_deg=90.0_real64 )
    t_west = optimal_tilt_deg( TOKYO_LAT, TOKYO_LON, panel_azimuth_deg=270.0_real64 )
    call check_range( 'SE (135) tilt <= south', t_se, 0.0_real64, t_south, nfail )
    call check_range( 'E/W near symmetric', t_east - t_west, -3.0_real64, 3.0_real64, nfail )
  end block

  ! ---- 直散分離（Engerer2/Erbs）+ Perez 傾斜面変換 ----
  block
    real(real64), parameter :: D2R = 3.14159265358979323846_real64 / 180.0_real64
    real(real64) :: noon, ghics, dni, dhi, sinh, recon, poa_h
    integer      :: doy_s
    doy_s = day_of_year( 2026, 6, 21 )
    noon  = 12.0_real64 + (135.0_real64 - TOKYO_LON) * 4.0_real64 / 60.0_real64  ! 太陽南中
    sinh  = sin( solar_elevation_deg( TOKYO_LAT, TOKYO_LON, doy_s, noon ) * D2R )

    ghics = clear_sky_ghi_wm2( TOKYO_LAT, TOKYO_LON, doy_s, noon )
    call check_range( 'clear-sky GHI plausible (W/m2)', ghics, 700.0_real64, 1050.0_real64, nfail )

    ! 雲量フォールバック（Kasten–Czeplak）: cc=0 → 晴天一致、cc=100 → 0.25×晴天、
    ! 中間は単調に減衰、夜間は 0。
    call check_range( 'cloudy(cc=0) == clear-sky', &
                      cloudy_sky_ghi_wm2( TOKYO_LAT, TOKYO_LON, doy_s, noon, 0.0_real64 ) - ghics, &
                      -0.5_real64, 0.5_real64, nfail )
    call check_range( 'cloudy(cc=100) == 0.25*clear-sky', &
                      cloudy_sky_ghi_wm2( TOKYO_LAT, TOKYO_LON, doy_s, noon, 100.0_real64 ) - 0.25_real64*ghics, &
                      -0.5_real64, 0.5_real64, nfail )
    call check_range( 'cloudy(cc=50) within (overcast, clear)', &
                      cloudy_sky_ghi_wm2( TOKYO_LAT, TOKYO_LON, doy_s, noon, 50.0_real64 ), &
                      0.25_real64*ghics, ghics, nfail )
    call check_range( 'cloudy at night == 0', &
                      cloudy_sky_ghi_wm2( TOKYO_LAT, TOKYO_LON, doy_s, 0.0_real64, 50.0_real64 ), &
                      -0.001_real64, 0.001_real64, nfail )

    ! 晴天相当 GHI → kt 高 → Engerer2 Kd 小・DNI 大、エネルギー保存
    call decompose_engerer2( ghics, TOKYO_LAT, TOKYO_LON, doy_s, noon, dni, dhi, interval_min = 0.0_real64 )
    call check_range( 'clear-day diffuse fraction small', dhi/ghics, 0.0_real64, 0.45_real64, nfail )
    call check_range( 'clear-day DNI substantial', dni, 400.0_real64, 1100.0_real64, nfail )
    recon = dhi + dni * sinh
    call check_range( 'engerer2 energy conservation', recon - ghics, -0.5_real64, 0.5_real64, nfail )

    ! Perez 傾斜面: 水平面 (tilt=0) は GHI に厳密一致
    poa_h = poa_perez( dni, dhi, ghics, TOKYO_LAT, TOKYO_LON, doy_s, noon, &
                       0.0_real64, 180.0_real64, interval_min = 0.0_real64 )
    call check_range( 'Perez horizontal == GHI', poa_h - ghics, -0.5_real64, 0.5_real64, nfail )

    ! 曇天 (低 GHI) → ほぼ全散乱・DNI≈0
    call decompose_engerer2( 80.0_real64, TOKYO_LAT, TOKYO_LON, doy_s, noon, dni, dhi, interval_min = 0.0_real64 )
    call check_range( 'overcast diffuse fraction high', dhi/80.0_real64, 0.70_real64, 1.0_real64, nfail )
    call check_range( 'overcast DNI small', dni, 0.0_real64, 120.0_real64, nfail )

    ! Erbs も水平面で GHI 一致（比較経路の健全性）
    call decompose_erbs( ghics, TOKYO_LAT, TOKYO_LON, doy_s, noon, dni, dhi, interval_min = 0.0_real64 )
    poa_h = poa_perez( dni, dhi, ghics, TOKYO_LAT, TOKYO_LON, doy_s, noon, &
                       0.0_real64, 180.0_real64, interval_min = 0.0_real64 )
    call check_range( 'Erbs+Perez horizontal == GHI', poa_h - ghics, -0.5_real64, 0.5_real64, nfail )
  end block

  print *, '----------------------------------------'
  if ( nfail == 0 ) then
    print *, 'ALL TESTS PASSED'
  else
    print '(1x,i0,a)', nfail, ' TEST(S) FAILED'
    error stop 1
  end if

contains

  subroutine check_i( name, got, want, nf )
    character(*), intent(in)    :: name
    integer,      intent(in)    :: got, want
    integer,      intent(inout) :: nf
    if ( got == want ) then
      print '(1x,a,a,a,i0)', 'PASS  ', name, ' = ', got
    else
      print '(1x,a,a,a,i0,a,i0)', 'FAIL  ', name, ' got ', got, ' want ', want
      nf = nf + 1
    end if
  end subroutine

  subroutine check_range( name, got, lo, hi, nf )
    character(*), intent(in)    :: name
    real(real64), intent(in)    :: got, lo, hi
    integer,      intent(inout) :: nf
    if ( got >= lo .and. got <= hi ) then
      print '(1x,a,a,a,f9.3)', 'PASS  ', name, ' = ', got
    else
      print '(1x,a,a,a,f9.3,a,f9.3,a,f9.3,a)', 'FAIL  ', name, ' = ', got, ' not in [', lo, ',', hi, ']'
      nf = nf + 1
    end if
  end subroutine

end program unit_test_solar_geometry_mo

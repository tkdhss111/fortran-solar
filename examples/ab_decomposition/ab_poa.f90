program ab_poa
  !! GHI A/B ハーネス: Engerer2 vs Erbs 直散分離 → Perez 傾斜面変換 で POA を比較。
  !!
  !! 標準入力 CSV（ヘッダ: site,jst,lat,lon,ghi）を読み、各行で両手法の POA を計算し
  !! CSV を標準出力へ、要約を標準エラーへ出力する。GHI は LFM の 30分平均を想定。
  !!
  !! 引数: [tilt_deg=30] [panel_azimuth_deg=180] [albedo=0.2]
  use, intrinsic :: iso_fortran_env, only: real64, input_unit, output_unit, error_unit
  use solar_geometry_mo
  implicit none

  character(1024) :: line, arg
  character(64)   :: site
  character(32)   :: jst
  real(real64) :: lat, lon, ghi, hour_jst, tilt, azi, alb
  real(real64) :: dni, dhi, poa_e2, poa_erbs, diff
  integer      :: doy, yr, mo, dy, hh, mi, ios, n, nargs
  real(real64) :: s_ghi, s_e2, s_erbs, s_absdiff, s_diff, s_sqdiff, mx_abs

  tilt = 30.0_real64 ; azi = 180.0_real64 ; alb = 0.2_real64
  nargs = command_argument_count()
  if ( nargs >= 1 ) then ; call get_command_argument(1, arg) ; read(arg,*,iostat=ios) tilt ; end if
  if ( nargs >= 2 ) then ; call get_command_argument(2, arg) ; read(arg,*,iostat=ios) azi  ; end if
  if ( nargs >= 3 ) then ; call get_command_argument(3, arg) ; read(arg,*,iostat=ios) alb  ; end if

  write(output_unit,'(a)') 'site,jst,ghi,dni_e2,dhi_e2,poa_engerer2,poa_erbs,diff'
  read(input_unit,'(a)',iostat=ios) line   ! skip header

  n = 0
  s_ghi = 0 ; s_e2 = 0 ; s_erbs = 0 ; s_absdiff = 0 ; s_diff = 0 ; s_sqdiff = 0 ; mx_abs = 0
  do
    read(input_unit,'(a)',iostat=ios) line
    if ( ios /= 0 ) exit
    if ( len_trim(line) == 0 ) cycle
    if ( .not. parse_row( line, site, jst, lat, lon, ghi ) ) cycle
    if ( ghi <= 0.0_real64 ) cycle
    if ( .not. parse_jst( jst, yr, mo, dy, hh, mi ) ) cycle
    doy      = day_of_year( yr, mo, dy )
    hour_jst = real(hh, real64) + real(mi, real64) / 60.0_real64

    ! 既定 Engerer2-30min + Perez（interval_min=30 既定）
    call decompose_engerer2( ghi, lat, lon, doy, hour_jst, dni, dhi )
    poa_e2   = poa_perez( dni, dhi, ghi, lat, lon, doy, hour_jst, tilt, azi, albedo = alb )
    ! 比較: Erbs + Perez
    poa_erbs = poa_from_ghi( ghi, lat, lon, doy, hour_jst, tilt, azi, albedo = alb, use_erbs = .true. )
    diff = poa_e2 - poa_erbs

    write(output_unit,'(a,",",a,",",f0.1,5(",",f0.2))') &
      trim(site), trim(jst), ghi, dni, dhi, poa_e2, poa_erbs, diff

    n = n + 1
    s_ghi = s_ghi + ghi ; s_e2 = s_e2 + poa_e2 ; s_erbs = s_erbs + poa_erbs
    s_diff = s_diff + diff ; s_absdiff = s_absdiff + abs(diff) ; s_sqdiff = s_sqdiff + diff*diff
    if ( abs(diff) > mx_abs ) mx_abs = abs(diff)
  end do

  if ( n > 0 ) then
    write(error_unit,'(a)')        '=== Engerer2 vs Erbs (-> Perez POA) A/B summary ==='
    write(error_unit,'(a,i0)')     '  daytime samples (GHI>0) : ', n
    write(error_unit,'(a,f0.1)')   '  mean GHI         [W/m2] : ', s_ghi / n
    write(error_unit,'(a,f0.1)')   '  mean POA Engerer2[W/m2] : ', s_e2 / n
    write(error_unit,'(a,f0.1)')   '  mean POA Erbs    [W/m2] : ', s_erbs / n
    write(error_unit,'(a,f0.2)')   '  mean diff (E2-Erbs)     : ', s_diff / n
    write(error_unit,'(a,f0.2)')   '  mean |diff|             : ', s_absdiff / n
    write(error_unit,'(a,f0.2)')   '  RMS diff                : ', sqrt(s_sqdiff / n)
    write(error_unit,'(a,f0.2)')   '  max |diff|              : ', mx_abs
    write(error_unit,'(a,f0.2,a)') '  mean |diff| as % of POA : ', 100.0_real64*s_absdiff/max(1.0e-9_real64,s_e2), ' %'
  else
    write(error_unit,'(a)') 'no daytime samples'
  end if

contains

  ! CSV 1 行を site,jst,lat,lon,ghi に分解（jst は空白を含むがカンマは無い前提）
  logical function parse_row( ln, site, jst, lat, lon, ghi )
    character(*), intent(in)  :: ln
    character(*), intent(out) :: site, jst
    real(real64), intent(out) :: lat, lon, ghi
    integer :: p(4), i, st, ios
    parse_row = .false.
    st = 1
    do i = 1, 4
      p(i) = index( ln(st:), ',' )
      if ( p(i) <= 0 ) return
      p(i) = p(i) + st - 1
      st = p(i) + 1
    end do
    site = ln(1:p(1)-1)
    jst  = ln(p(1)+1:p(2)-1)
    read( ln(p(2)+1:p(3)-1), *, iostat=ios ) lat ; if ( ios /= 0 ) return
    read( ln(p(3)+1:p(4)-1), *, iostat=ios ) lon ; if ( ios /= 0 ) return
    read( ln(p(4)+1:),       *, iostat=ios ) ghi ; if ( ios /= 0 ) return
    parse_row = .true.
  end function

  ! "YYYY-MM-DD HH:MM:SS" を分解
  logical function parse_jst( s, yr, mo, dy, hh, mi )
    character(*), intent(in)  :: s
    integer,      intent(out) :: yr, mo, dy, hh, mi
    integer :: ios
    parse_jst = .false.
    if ( len_trim(s) < 16 ) return
    read( s(1:4),   *, iostat=ios ) yr ; if ( ios /= 0 ) return
    read( s(6:7),   *, iostat=ios ) mo ; if ( ios /= 0 ) return
    read( s(9:10),  *, iostat=ios ) dy ; if ( ios /= 0 ) return
    read( s(12:13), *, iostat=ios ) hh ; if ( ios /= 0 ) return
    read( s(15:16), *, iostat=ios ) mi ; if ( ios /= 0 ) return
    parse_jst = .true.
  end function

end program ab_poa

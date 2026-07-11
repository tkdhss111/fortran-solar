
module solar_geometry_mo
  !! Spencer (1971) 式による太陽位置計算モジュール。
  !! 角度は特に記載がない限り度（degrees）。_rad 接尾辞はラジアン。
  !! コード規約: ユーザ定義手続きの呼出し/定義は括弧内に空白 ( a1, a2 ) を入れ、
  !!            配列添字 arr(i,j) ・型種別 real(real64) ・組込関数 sin(x) と区別する。
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: solar_elevation_deg, solar_azimuth_deg, &
            extraterrestrial_radiation_wm2, day_of_year, &
            incidence_angle_deg, optimal_tilt_deg, &
            apparent_solar_time_h, clear_sky_ghi_wm2, cloudy_sky_ghi_wm2, &
            air_mass_kastenyoung, &
            decompose_erbs, decompose_engerer2, poa_perez, poa_from_ghi

  real(real64), parameter :: PI      = 3.14159265358979323846_real64
  real(real64), parameter :: DEG2RAD = PI / 180.0_real64
  real(real64), parameter :: RAD2DEG = 180.0_real64 / PI
  real(real64), parameter :: SOLAR_CONSTANT = 1361.0_real64  ! W/m2 (太陽定数)
  ! 年間平均の散乱日射比率（散乱/全天）。日本の平年値は概ね 0.45–0.55。
  ! 最適傾斜角の既定値に使用（引数で上書き可）。
  real(real64), parameter :: DIFFUSE_FRACTION_DEFAULT = 0.5_real64

contains

  !-----------------------------------------------------------------
  ! 年間通日（1-366）。閏年対応。
  !-----------------------------------------------------------------
  pure integer function day_of_year( year, month, day )
    integer, intent(in) :: year, month, day
    integer :: m
    integer, parameter :: mdays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    day_of_year = day
    do m = 1, month - 1
      day_of_year = day_of_year + mdays(m)
    end do
    if ( month > 2 .and. is_leap( year ) ) day_of_year = day_of_year + 1
  end function

  !-----------------------------------------------------------------
  ! 閏年判定
  !-----------------------------------------------------------------
  pure logical function is_leap( year )
    integer, intent(in) :: year
    is_leap = ( mod(year,4) == 0 .and. mod(year,100) /= 0 ) .or. mod(year,400) == 0
  end function

  !-----------------------------------------------------------------
  ! 日角（ラジアン）— Spencer (1971)
  !-----------------------------------------------------------------
  pure real(real64) function day_angle_rad( doy )
    integer, intent(in) :: doy
    day_angle_rad = 2.0_real64 * PI * real(doy - 1, real64) / 365.0_real64
  end function

  !-----------------------------------------------------------------
  ! 太陽赤緯（ラジアン）— Spencer (1971)
  !-----------------------------------------------------------------
  pure real(real64) function solar_declination_rad( doy )
    integer, intent(in) :: doy
    real(real64) :: b
    b = day_angle_rad( doy )
    solar_declination_rad = 0.006918_real64 &
      - 0.399912_real64 * cos(b) + 0.070257_real64 * sin(b) &
      - 0.006758_real64 * cos(2.0_real64*b) + 0.000907_real64 * sin(2.0_real64*b) &
      - 0.002697_real64 * cos(3.0_real64*b) + 0.001480_real64 * sin(3.0_real64*b)
  end function

  !-----------------------------------------------------------------
  ! 均時差（分）— Spencer (1971)
  !-----------------------------------------------------------------
  pure real(real64) function equation_of_time_min( doy )
    integer, intent(in) :: doy
    real(real64) :: b
    b = day_angle_rad( doy )
    equation_of_time_min = 229.18_real64 * ( &
        0.000075_real64 &
      + 0.001868_real64 * cos(b) - 0.032077_real64 * sin(b) &
      - 0.014615_real64 * cos(2.0_real64*b) - 0.040849_real64 * sin(2.0_real64*b) )
  end function

  !-----------------------------------------------------------------
  ! 時角（ラジアン）
  !   hour_jst: JST 時刻（0-23、小数可）
  !   lon_deg:  経度（東経、度）
  !-----------------------------------------------------------------
  pure real(real64) function hour_angle_rad( hour_jst, lon_deg, doy )
    real(real64), intent(in) :: hour_jst, lon_deg
    integer,      intent(in) :: doy
    real(real64) :: solar_time_h
    ! JST 標準子午線は東経 135°
    ! 太陽時 = 地方標準時 + 均時差/60 + (経度 - 標準子午線) × 4/60
    solar_time_h = hour_jst &
      + equation_of_time_min( doy ) / 60.0_real64 &
      + (lon_deg - 135.0_real64) * 4.0_real64 / 60.0_real64
    hour_angle_rad = (solar_time_h - 12.0_real64) * 15.0_real64 * DEG2RAD
  end function

  !-----------------------------------------------------------------
  ! 太陽高度角（度）
  !-----------------------------------------------------------------
  pure real(real64) function solar_elevation_deg( lat_deg, lon_deg, doy, hour_jst )
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst
    integer,      intent(in) :: doy
    real(real64) :: lat_r, dec_r, ha_r, sin_elev
    lat_r = lat_deg * DEG2RAD
    dec_r = solar_declination_rad( doy )
    ha_r  = hour_angle_rad( hour_jst, lon_deg, doy )
    sin_elev = sin(lat_r) * sin(dec_r) + cos(lat_r) * cos(dec_r) * cos(ha_r)
    solar_elevation_deg = asin(max(-1.0_real64, min(1.0_real64, sin_elev))) * RAD2DEG
  end function

  !-----------------------------------------------------------------
  ! 太陽方位角（度、北から時計回り）
  !-----------------------------------------------------------------
  pure real(real64) function solar_azimuth_deg( lat_deg, lon_deg, doy, hour_jst )
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst
    integer,      intent(in) :: doy
    real(real64) :: lat_r, dec_r, ha_r, sin_elev, cos_elev
    real(real64) :: cos_az, sin_az, az
    lat_r = lat_deg * DEG2RAD
    dec_r = solar_declination_rad( doy )
    ha_r  = hour_angle_rad( hour_jst, lon_deg, doy )
    sin_elev = sin(lat_r) * sin(dec_r) + cos(lat_r) * cos(dec_r) * cos(ha_r)
    cos_elev = sqrt(max(0.0_real64, 1.0_real64 - sin_elev**2))
    if ( cos_elev < 1.0e-10_real64 ) then
      solar_azimuth_deg = 0.0_real64
      return
    end if
    cos_az = (sin(dec_r) - sin_elev * sin(lat_r)) / (cos_elev * cos(lat_r))
    cos_az = max(-1.0_real64, min(1.0_real64, cos_az))
    sin_az = -cos(dec_r) * sin(ha_r) / cos_elev
    az = acos(cos_az) * RAD2DEG
    if ( sin_az < 0.0_real64 ) az = 360.0_real64 - az
    solar_azimuth_deg = az
  end function

  !-----------------------------------------------------------------
  ! 大気外水平面日射量（W/m2）
  !   太陽が地平線以下のときは 0 を返す。
  !-----------------------------------------------------------------
  pure real(real64) function extraterrestrial_radiation_wm2( lat_deg, lon_deg, doy, hour_jst )
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst
    integer,      intent(in) :: doy
    real(real64) :: b, eccentricity, elev
    b = day_angle_rad( doy )
    ! 地球-太陽距離補正係数 — Spencer (1971)
    eccentricity = 1.000110_real64 &
      + 0.034221_real64 * cos(b) + 0.001280_real64 * sin(b) &
      + 0.000719_real64 * cos(2.0_real64*b) + 0.000077_real64 * sin(2.0_real64*b)
    elev = solar_elevation_deg( lat_deg, lon_deg, doy, hour_jst )
    if ( elev <= 0.0_real64 ) then
      extraterrestrial_radiation_wm2 = 0.0_real64
    else
      extraterrestrial_radiation_wm2 = SOLAR_CONSTANT * eccentricity * sin(elev * DEG2RAD)
    end if
  end function

  !-----------------------------------------------------------------
  ! 任意方位・傾斜面への太陽光線の入射角余弦 cos(θ)。
  !   tilt_deg          : 傾斜角（水平=0, 鉛直=90, 度）
  !   panel_azimuth_deg : パネル方位（北=0, 東=90, 南=180, 西=270;
  !                       solar_azimuth_deg と同基準）
  !   cosθ = sin(h)·cos(β) + cos(h)·sin(β)·cos(γs - γp)
  !          h=太陽高度角, β=傾斜角, γs=太陽方位, γp=パネル方位
  !-----------------------------------------------------------------
  pure real(real64) function cos_incidence( lat_deg, lon_deg, doy, hour_jst, tilt_deg, panel_azimuth_deg )
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst, tilt_deg, panel_azimuth_deg
    integer,      intent(in) :: doy
    real(real64) :: elev_r, saz_deg, b_r
    elev_r  = solar_elevation_deg( lat_deg, lon_deg, doy, hour_jst ) * DEG2RAD
    saz_deg = solar_azimuth_deg( lat_deg, lon_deg, doy, hour_jst )
    b_r     = tilt_deg * DEG2RAD
    cos_incidence = sin(elev_r) * cos(b_r) &
                  + cos(elev_r) * sin(b_r) * cos( (saz_deg - panel_azimuth_deg) * DEG2RAD )
  end function

  !-----------------------------------------------------------------
  ! 入射角（度）— 任意方位・傾斜面と太陽光線のなす角。
  !   傾斜面日射（transposition: 等方/Hay-Davies/Perez 等）の基本量。
  !   背面側（太陽がパネル裏）では 90° 超を返す。
  !-----------------------------------------------------------------
  pure real(real64) function incidence_angle_deg( lat_deg, lon_deg, doy, hour_jst, tilt_deg, panel_azimuth_deg )
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst, tilt_deg, panel_azimuth_deg
    integer,      intent(in) :: doy
    real(real64) :: ci
    ci = cos_incidence( lat_deg, lon_deg, doy, hour_jst, tilt_deg, panel_azimuth_deg )
    incidence_angle_deg = acos( max(-1.0_real64, min(1.0_real64, ci)) ) * RAD2DEG
  end function

  !-----------------------------------------------------------------
  ! 任意方位・傾斜面の年間傾斜面日射量（相対値, isotropic sky モデル）。
  !   最適傾斜角探索用の設計時近似。大気外水平面日射を全天日射の代理に、
  !   散乱比率 fd で直達/散乱に配分する。
  !   ※ 予報精度向上の transposition（Perez 等）は別途実装予定。
  !-----------------------------------------------------------------
  pure real(real64) function annual_poa( lat_deg, lon_deg, tilt_deg, panel_azimuth_deg, fd )
    real(real64), intent(in) :: lat_deg, lon_deg, tilt_deg, panel_azimuth_deg, fd
    integer      :: doy, ih
    real(real64) :: hour, elev, sinelev, g0h, cosi, beam, diff, view_tilt, total
    view_tilt = (1.0_real64 + cos(tilt_deg * DEG2RAD)) * 0.5_real64   ! 天空率 (1+cosβ)/2
    total = 0.0_real64
    do doy = 1, 365
      do ih = 0, 23
        hour = real(ih, real64)
        elev = solar_elevation_deg( lat_deg, lon_deg, doy, hour )
        if ( elev <= 0.0_real64 ) cycle                ! 夜間はスキップ
        sinelev = sin(elev * DEG2RAD)
        g0h     = extraterrestrial_radiation_wm2( lat_deg, lon_deg, doy, hour )  ! 水平面・大気外
        cosi    = cos_incidence( lat_deg, lon_deg, doy, hour, tilt_deg, panel_azimuth_deg )
        ! 直達(傾斜面) = GHI(1-fd) * Rb,  Rb = cosθ/sinh (背面側 cosθ<0 は 0)
        beam = (1.0_real64 - fd) * g0h * max(0.0_real64, cosi) / sinelev
        ! 散乱(等方天空) = GHI*fd * (1+cosβ)/2
        diff = fd * g0h * view_tilt
        total = total + beam + diff
      end do
    end do
    annual_poa = total
  end function

  !-----------------------------------------------------------------
  ! 最適傾斜角（度）— 指定方位の固定設置で年間日射量を最大化する傾斜角。
  !   lat_deg / lon_deg : 設置地点（北緯・東経, 度）
  !   panel_azimuth_deg : 任意。パネル方位（既定 180 = 真南）。
  !                       北=0, 東=90, 南=180, 西=270。
  !   diffuse_fraction  : 任意。年間散乱日射比率（既定 0.5）。
  !   返り値            : 0–90 度（1 度刻みで探索）。
  !   ※ 等方天空・大気外日射近似による設計時推定。予報用の高精度
  !     transposition（Perez 等）は別途。
  !-----------------------------------------------------------------
  pure real(real64) function optimal_tilt_deg( lat_deg, lon_deg, panel_azimuth_deg, diffuse_fraction )
    real(real64), intent(in)           :: lat_deg, lon_deg
    real(real64), intent(in), optional :: panel_azimuth_deg, diffuse_fraction
    real(real64) :: az, fd, poa, best_poa
    integer      :: i, best_i
    az = 180.0_real64                                  ! 既定: 真南
    if ( present(panel_azimuth_deg) ) az = panel_azimuth_deg
    fd = DIFFUSE_FRACTION_DEFAULT
    if ( present(diffuse_fraction) ) fd = diffuse_fraction
    best_poa = -1.0_real64
    best_i   = 0
    do i = 0, 90
      poa = annual_poa( lat_deg, lon_deg, real(i, real64), az, fd )
      if ( poa > best_poa ) then
        best_poa = poa
        best_i   = i
      end if
    end do
    optimal_tilt_deg = real(best_i, real64)
  end function

  !=================================================================
  ! PV 傾斜面日射: GHI → DNI/DHI 分離 → Perez 傾斜面変換
  !
  ! 入力 GHI は LFM GPV の 30 分平均（区間 [hour-Δ, hour]）を想定する。
  ! 幾何量は区間代表点（中点 hour - Δ/2）で評価する（平滑関数の区間平均
  ! ≒ 中点値）。interval_min はその Δ（既定 30 分; 0 で瞬時値）。
  !=================================================================

  !-----------------------------------------------------------------
  ! 視太陽時（hours, 0–24）。AST = 12 + 時角(度)/15。
  !-----------------------------------------------------------------
  pure real(real64) function apparent_solar_time_h( lon_deg, doy, hour_jst )
    real(real64), intent(in) :: lon_deg, hour_jst
    integer,      intent(in) :: doy
    apparent_solar_time_h = 12.0_real64 &
      + hour_angle_rad( hour_jst, lon_deg, doy ) * RAD2DEG / 15.0_real64
  end function

  !-----------------------------------------------------------------
  ! 晴天時全天日射（W/m2, ASHRAE / Threlkeld–Jordan "TJ" モデル）。
  ! Engerer2 30分版の元適合と同一の晴天モデル（dktc/kde の整合のため）。
  ! 係数 A,k,C は doy の正弦関数（角度は度）:
  !   A = 1160 + 75·sin(360(doy-275)/365)   [W/m2]
  !   k = 0.174 + 0.035·sin(360(doy-100)/365)
  !   C = 0.095 + 0.040·sin(360(doy-100)/365)
  ! DNIcs = A·exp(-k/cosZ),  GHIcs = DNIcs·(cosZ + C),  cosZ = sin(高度角)。
  !-----------------------------------------------------------------
  pure real(real64) function clear_sky_ghi_wm2( lat_deg, lon_deg, doy, hour_jst )
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst
    integer,      intent(in) :: doy
    real(real64) :: elev, cosz, a, k, c, dnics, dd
    elev = solar_elevation_deg( lat_deg, lon_deg, doy, hour_jst )
    if ( elev <= 0.0_real64 ) then
      clear_sky_ghi_wm2 = 0.0_real64
      return
    end if
    cosz = sin(elev * DEG2RAD)                          ! cos(天頂角) = sin(高度角)
    dd = real(doy, real64)
    a = 1160.0_real64 + 75.0_real64  * sin( 360.0_real64*(dd-275.0_real64)/365.0_real64 * DEG2RAD )
    k = 0.174_real64  + 0.035_real64 * sin( 360.0_real64*(dd-100.0_real64)/365.0_real64 * DEG2RAD )
    c = 0.095_real64  + 0.040_real64 * sin( 360.0_real64*(dd-100.0_real64)/365.0_real64 * DEG2RAD )
    dnics = a * exp( -k / cosz )
    clear_sky_ghi_wm2 = dnics * ( cosz + c )
  end function

  !-----------------------------------------------------------------
  ! 雲量を考慮した全天日射 GHI [W/m2]。Kasten–Czeplak (1980):
  !   GHI = GHIcs · ( 1 − 0.75·N^3.4 ),  N = 全雲量（0–1）= cloud_cover_pct/100
  ! GHIcs は clear_sky_ghi_wm2。夜間は GHIcs=0 のため自動的に 0 を返す。
  ! cloud_cover_pct は 0–100 [%]（範囲外は 0–100 に clamp）。
  ! 全雲量のみから GHI を推定する汎用手法（例: 簡易日射推定、日射欠測の物理補間）。
  !-----------------------------------------------------------------
  pure real(real64) function cloudy_sky_ghi_wm2( lat_deg, lon_deg, doy, hour_jst, cloud_cover_pct )
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst, cloud_cover_pct
    integer,      intent(in) :: doy
    real(real64) :: ghics, n
    ghics = clear_sky_ghi_wm2( lat_deg, lon_deg, doy, hour_jst )
    n = max( 0.0_real64, min( 100.0_real64, cloud_cover_pct ) ) / 100.0_real64
    cloudy_sky_ghi_wm2 = ghics * ( 1.0_real64 - 0.75_real64 * n**3.4_real64 )
  end function

  !-----------------------------------------------------------------
  ! 相対大気質量（Kasten–Young 1989）。zenith_deg = 天頂角（度）。
  !-----------------------------------------------------------------
  pure real(real64) function air_mass_kastenyoung( zenith_deg )
    real(real64), intent(in) :: zenith_deg
    real(real64) :: z
    z = min( zenith_deg, 90.0_real64 )
    air_mass_kastenyoung = 1.0_real64 / &
      ( cos(z * DEG2RAD) + 0.50572_real64 * (96.07995_real64 - z)**(-1.6364_real64) )
  end function

  !-----------------------------------------------------------------
  ! Perez(1990) 天空清明度 ε のビン番号（1–8）。
  !-----------------------------------------------------------------
  pure integer function perez_bin( eps )
    real(real64), intent(in) :: eps
    if      ( eps < 1.065_real64 ) then ; perez_bin = 1
    else if ( eps < 1.230_real64 ) then ; perez_bin = 2
    else if ( eps < 1.500_real64 ) then ; perez_bin = 3
    else if ( eps < 1.950_real64 ) then ; perez_bin = 4
    else if ( eps < 2.800_real64 ) then ; perez_bin = 5
    else if ( eps < 4.500_real64 ) then ; perez_bin = 6
    else if ( eps < 6.200_real64 ) then ; perez_bin = 7
    else                                ; perez_bin = 8
    end if
  end function

  !-----------------------------------------------------------------
  ! Erbs(1982) 直散分離。GHI[W/m2] → DNI, DHI[W/m2]。
  !-----------------------------------------------------------------
  pure subroutine decompose_erbs( ghi, lat_deg, lon_deg, doy, hour_jst, dni, dhi, interval_min )
    real(real64), intent(in)            :: ghi, lat_deg, lon_deg, hour_jst
    integer,      intent(in)            :: doy
    real(real64), intent(out)           :: dni, dhi
    real(real64), intent(in), optional  :: interval_min
    real(real64) :: iv, t_rep, elev, sinelev, g0h, kt, df
    dni = 0.0_real64 ; dhi = 0.0_real64
    iv = 30.0_real64 ; if ( present(interval_min) ) iv = interval_min
    t_rep = hour_jst - iv / 120.0_real64
    elev = solar_elevation_deg( lat_deg, lon_deg, doy, t_rep )
    if ( elev <= 0.0_real64 .or. ghi <= 0.0_real64 ) return
    sinelev = sin(elev * DEG2RAD)
    g0h = extraterrestrial_radiation_wm2( lat_deg, lon_deg, doy, t_rep )
    if ( g0h <= 0.0_real64 ) return
    kt = max( 0.0_real64, min( 1.0_real64, ghi / g0h ) )
    if      ( kt <= 0.22_real64 ) then
      df = 1.0_real64 - 0.09_real64 * kt
    else if ( kt <= 0.80_real64 ) then
      df = 0.9511_real64 - 0.1604_real64*kt + 4.388_real64*kt**2 &
         - 16.638_real64*kt**3 + 12.336_real64*kt**4
    else
      df = 0.165_real64
    end if
    dhi = df * ghi
    dni = (ghi - dhi) / sinelev
  end subroutine

  !-----------------------------------------------------------------
  ! Engerer2 直散分離（30 分版係数, Bright & Engerer 2019）。
  !   Kd = C + (1-C)/(1+exp(B0+B1·kt+B2·ast+B3·zen+B4·dktc)) + B5·kde
  !   kt=GHI/G0h, zen=天頂角[rad], ast=視太陽時[h], dktc=ktc-kt,
  !   kde=max(0,(GHI-GHIcs)/GHI)。DHI=GHI·Kd, DNI=(GHI-DHI)/cos(zen)。
  !-----------------------------------------------------------------
  pure subroutine decompose_engerer2( ghi, lat_deg, lon_deg, doy, hour_jst, dni, dhi, interval_min )
    real(real64), intent(in)            :: ghi, lat_deg, lon_deg, hour_jst
    integer,      intent(in)            :: doy
    real(real64), intent(out)           :: dni, dhi
    real(real64), intent(in), optional  :: interval_min
    real(real64), parameter :: C  = 0.0326750_real64, B0 = -4.8681_real64, &
                               B1 = 8.1867_real64,    B2 = 0.015829_real64, &
                               B3 = 0.0059922_real64, B4 = -4.0304_real64,  &
                               B5 = 0.47371_real64
    real(real64) :: iv, t_rep, elev, sinelev, zen_r, g0h, ghics, kt, ktc, dktc, kde, ast, kd
    dni = 0.0_real64 ; dhi = 0.0_real64
    iv = 30.0_real64 ; if ( present(interval_min) ) iv = interval_min
    t_rep = hour_jst - iv / 120.0_real64
    elev = solar_elevation_deg( lat_deg, lon_deg, doy, t_rep )
    if ( elev <= 0.0_real64 .or. ghi <= 0.0_real64 ) return
    sinelev = sin(elev * DEG2RAD)
    zen_r   = (90.0_real64 - elev) * DEG2RAD
    g0h = extraterrestrial_radiation_wm2( lat_deg, lon_deg, doy, t_rep )
    if ( g0h <= 0.0_real64 ) return
    kt    = ghi / g0h
    ghics = clear_sky_ghi_wm2( lat_deg, lon_deg, doy, t_rep )
    ktc   = ghics / g0h
    dktc  = ktc - kt
    kde   = 0.0_real64
    if ( ghi > ghics ) kde = (ghi - ghics) / ghi
    ast = apparent_solar_time_h( lon_deg, doy, t_rep )
    kd  = C + (1.0_real64 - C) / &
          ( 1.0_real64 + exp( B0 + B1*kt + B2*ast + B3*zen_r + B4*dktc ) ) + B5*kde
    kd  = max( 0.0_real64, min( 1.0_real64, kd ) )
    dhi = kd * ghi
    dni = (ghi - dhi) / sinelev
  end subroutine

  !-----------------------------------------------------------------
  ! Perez(1990) 傾斜面全天日射（W/m2, allsitescomposite1990 係数）。
  !   POA = DNI·cosθ + DHI·R_d(Perez) + GHI·ρ·(1-cosβ)/2
  !   albedo ρ 既定 0.2。tilt_deg/panel_azimuth_deg は incidence と同基準。
  !-----------------------------------------------------------------
  pure real(real64) function poa_perez( dni, dhi, ghi, lat_deg, lon_deg, doy, hour_jst, &
                                        tilt_deg, panel_azimuth_deg, albedo, interval_min )
    real(real64), intent(in)           :: dni, dhi, ghi, lat_deg, lon_deg, hour_jst, &
                                          tilt_deg, panel_azimuth_deg
    integer,      intent(in)           :: doy
    real(real64), intent(in), optional :: albedo, interval_min
    ! Perez allsitescomposite1990: 各 bin の [F11,F12,F13,F21,F22,F23]
    real(real64), parameter :: PF(6,8) = reshape([ &
      -0.0083117_real64,  0.5877285_real64, -0.0620636_real64, -0.0596012_real64,  0.0721249_real64, -0.0220216_real64, &
       0.1299457_real64,  0.6825954_real64, -0.1513752_real64, -0.0189325_real64,  0.0659650_real64, -0.0288748_real64, &
       0.3296958_real64,  0.4868735_real64, -0.2210958_real64,  0.0554140_real64, -0.0639588_real64, -0.0260542_real64, &
       0.5682053_real64,  0.1874525_real64, -0.2951290_real64,  0.1088631_real64, -0.1519229_real64, -0.0139754_real64, &
       0.8730280_real64, -0.3920403_real64, -0.3616149_real64,  0.2255647_real64, -0.4620442_real64,  0.0012448_real64, &
       1.1326077_real64, -1.2367284_real64, -0.4118494_real64,  0.2877813_real64, -0.8230357_real64,  0.0558651_real64, &
       1.0601591_real64, -1.5999137_real64, -0.3589221_real64,  0.2642124_real64, -1.1272340_real64,  0.1310694_real64, &
       0.6777470_real64, -0.3272588_real64, -0.2504286_real64,  0.1561313_real64, -1.3765031_real64,  0.2506212_real64  &
      ], [6,8])
    real(real64) :: rho, iv, t_rep, elev, b_r, ground, zen_deg, zen_r, cosi, &
                    g0n, eps, delta, am, f1, f2, a_term, b_term, sky, beam, kappa
    integer :: ib
    rho = 0.2_real64  ; if ( present(albedo) )       rho = albedo
    iv  = 30.0_real64 ; if ( present(interval_min) ) iv  = interval_min
    t_rep = hour_jst - iv / 120.0_real64
    b_r   = tilt_deg * DEG2RAD
    ground = max( 0.0_real64, ghi ) * rho * (1.0_real64 - cos(b_r)) * 0.5_real64
    elev = solar_elevation_deg( lat_deg, lon_deg, doy, t_rep )
    if ( elev <= 0.0_real64 ) then                     ! 夜間: 地面反射のみ
      poa_perez = ground
      return
    end if
    zen_deg = 90.0_real64 - elev
    zen_r   = zen_deg * DEG2RAD
    cosi = cos( incidence_angle_deg( lat_deg, lon_deg, doy, t_rep, tilt_deg, panel_azimuth_deg ) * DEG2RAD )
    beam = dni * max( 0.0_real64, cosi )
    sky  = 0.0_real64
    if ( dhi > 0.0_real64 ) then
      kappa = 1.041_real64
      eps   = ( (dhi + dni)/dhi + kappa*zen_r**3 ) / ( 1.0_real64 + kappa*zen_r**3 )
      am    = air_mass_kastenyoung( zen_deg )
      g0n   = extraterrestrial_radiation_wm2( lat_deg, lon_deg, doy, t_rep ) / max( 1.0e-6_real64, sin(elev*DEG2RAD) )
      delta = dhi * am / g0n
      ib    = perez_bin( eps )
      f1    = max( 0.0_real64, PF(1,ib) + PF(2,ib)*delta + PF(3,ib)*zen_r )
      f2    =                  PF(4,ib) + PF(5,ib)*delta + PF(6,ib)*zen_r
      a_term = max( 0.0_real64, cosi )
      b_term = max( cos(85.0_real64*DEG2RAD), cos(zen_r) )
      sky = dhi * ( 0.5_real64*(1.0_real64 - f1)*(1.0_real64 + cos(b_r)) &
                  + f1*a_term/b_term + f2*sin(b_r) )
      sky = max( 0.0_real64, sky )
    end if
    poa_perez = beam + sky + ground
  end function

  !-----------------------------------------------------------------
  ! GHI → 傾斜面全天日射 POA（既定: Engerer2-30min 分離 + Perez 変換）。
  !   use_erbs=.true. で分離を Erbs に切替（比較用）。
  !-----------------------------------------------------------------
  pure real(real64) function poa_from_ghi( ghi, lat_deg, lon_deg, doy, hour_jst, &
                                           tilt_deg, panel_azimuth_deg, albedo, interval_min, use_erbs )
    real(real64), intent(in)           :: ghi, lat_deg, lon_deg, hour_jst, tilt_deg, panel_azimuth_deg
    integer,      intent(in)           :: doy
    real(real64), intent(in), optional :: albedo, interval_min
    logical,      intent(in), optional :: use_erbs
    real(real64) :: dni, dhi, iv, alb
    logical      :: erbs
    iv   = 30.0_real64  ; if ( present(interval_min) ) iv  = interval_min
    alb  = 0.2_real64   ; if ( present(albedo) )       alb = albedo
    erbs = .false.      ; if ( present(use_erbs) )     erbs = use_erbs
    if ( erbs ) then
      call decompose_erbs( ghi, lat_deg, lon_deg, doy, hour_jst, dni, dhi, interval_min = iv )
    else
      call decompose_engerer2( ghi, lat_deg, lon_deg, doy, hour_jst, dni, dhi, interval_min = iv )
    end if
    poa_from_ghi = poa_perez( dni, dhi, ghi, lat_deg, lon_deg, doy, hour_jst, &
                              tilt_deg, panel_azimuth_deg, albedo = alb, interval_min = iv )
  end function

end module solar_geometry_mo

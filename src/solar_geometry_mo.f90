
module solar_geometry_mo
  !! Spencer (1971) 式による太陽位置計算モジュール。
  !! 角度は特に記載がない限り度（degrees）。_rad 接尾辞はラジアン。
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: solar_elevation_deg, solar_azimuth_deg, &
            extraterrestrial_radiation_wm2, day_of_year

  real(real64), parameter :: PI      = 3.14159265358979323846_real64
  real(real64), parameter :: DEG2RAD = PI / 180.0_real64
  real(real64), parameter :: RAD2DEG = 180.0_real64 / PI
  real(real64), parameter :: SOLAR_CONSTANT = 1361.0_real64  ! W/m2 (太陽定数)

contains

  !-----------------------------------------------------------------
  ! 年間通日（1-366）。閏年対応。
  !-----------------------------------------------------------------
  pure integer function day_of_year(year, month, day)
    integer, intent(in) :: year, month, day
    integer :: m
    integer, parameter :: mdays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    day_of_year = day
    do m = 1, month - 1
      day_of_year = day_of_year + mdays(m)
    end do
    if (month > 2 .and. is_leap(year)) day_of_year = day_of_year + 1
  end function

  !-----------------------------------------------------------------
  ! 閏年判定
  !-----------------------------------------------------------------
  pure logical function is_leap(year)
    integer, intent(in) :: year
    is_leap = (mod(year,4) == 0 .and. mod(year,100) /= 0) .or. mod(year,400) == 0
  end function

  !-----------------------------------------------------------------
  ! 日角（ラジアン）— Spencer (1971)
  !-----------------------------------------------------------------
  pure real(real64) function day_angle_rad(doy)
    integer, intent(in) :: doy
    day_angle_rad = 2.0_real64 * PI * real(doy - 1, real64) / 365.0_real64
  end function

  !-----------------------------------------------------------------
  ! 太陽赤緯（ラジアン）— Spencer (1971)
  !-----------------------------------------------------------------
  pure real(real64) function solar_declination_rad(doy)
    integer, intent(in) :: doy
    real(real64) :: b
    b = day_angle_rad(doy)
    solar_declination_rad = 0.006918_real64 &
      - 0.399912_real64 * cos(b) + 0.070257_real64 * sin(b) &
      - 0.006758_real64 * cos(2.0_real64*b) + 0.000907_real64 * sin(2.0_real64*b) &
      - 0.002697_real64 * cos(3.0_real64*b) + 0.001480_real64 * sin(3.0_real64*b)
  end function

  !-----------------------------------------------------------------
  ! 均時差（分）— Spencer (1971)
  !-----------------------------------------------------------------
  pure real(real64) function equation_of_time_min(doy)
    integer, intent(in) :: doy
    real(real64) :: b
    b = day_angle_rad(doy)
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
  pure real(real64) function hour_angle_rad(hour_jst, lon_deg, doy)
    real(real64), intent(in) :: hour_jst, lon_deg
    integer,      intent(in) :: doy
    real(real64) :: solar_time_h
    ! JST 標準子午線は東経 135°
    ! 太陽時 = 地方標準時 + 均時差/60 + (経度 - 標準子午線) × 4/60
    solar_time_h = hour_jst &
      + equation_of_time_min(doy) / 60.0_real64 &
      + (lon_deg - 135.0_real64) * 4.0_real64 / 60.0_real64
    hour_angle_rad = (solar_time_h - 12.0_real64) * 15.0_real64 * DEG2RAD
  end function

  !-----------------------------------------------------------------
  ! 太陽高度角（度）
  !-----------------------------------------------------------------
  pure real(real64) function solar_elevation_deg(lat_deg, lon_deg, doy, hour_jst)
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst
    integer,      intent(in) :: doy
    real(real64) :: lat_r, dec_r, ha_r, sin_elev
    lat_r = lat_deg * DEG2RAD
    dec_r = solar_declination_rad(doy)
    ha_r  = hour_angle_rad(hour_jst, lon_deg, doy)
    sin_elev = sin(lat_r) * sin(dec_r) + cos(lat_r) * cos(dec_r) * cos(ha_r)
    solar_elevation_deg = asin(max(-1.0_real64, min(1.0_real64, sin_elev))) * RAD2DEG
  end function

  !-----------------------------------------------------------------
  ! 太陽方位角（度、北から時計回り）
  !-----------------------------------------------------------------
  pure real(real64) function solar_azimuth_deg(lat_deg, lon_deg, doy, hour_jst)
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst
    integer,      intent(in) :: doy
    real(real64) :: lat_r, dec_r, ha_r, sin_elev, cos_elev
    real(real64) :: cos_az, sin_az, az
    lat_r = lat_deg * DEG2RAD
    dec_r = solar_declination_rad(doy)
    ha_r  = hour_angle_rad(hour_jst, lon_deg, doy)
    sin_elev = sin(lat_r) * sin(dec_r) + cos(lat_r) * cos(dec_r) * cos(ha_r)
    cos_elev = sqrt(max(0.0_real64, 1.0_real64 - sin_elev**2))
    if (cos_elev < 1.0e-10_real64) then
      solar_azimuth_deg = 0.0_real64
      return
    end if
    cos_az = (sin(dec_r) - sin_elev * sin(lat_r)) / (cos_elev * cos(lat_r))
    cos_az = max(-1.0_real64, min(1.0_real64, cos_az))
    sin_az = -cos(dec_r) * sin(ha_r) / cos_elev
    az = acos(cos_az) * RAD2DEG
    if (sin_az < 0.0_real64) az = 360.0_real64 - az
    solar_azimuth_deg = az
  end function

  !-----------------------------------------------------------------
  ! 大気外水平面日射量（W/m2）
  !   太陽が地平線以下のときは 0 を返す。
  !-----------------------------------------------------------------
  pure real(real64) function extraterrestrial_radiation_wm2(lat_deg, lon_deg, doy, hour_jst)
    real(real64), intent(in) :: lat_deg, lon_deg, hour_jst
    integer,      intent(in) :: doy
    real(real64) :: b, eccentricity, elev
    b = day_angle_rad(doy)
    ! 地球-太陽距離補正係数 — Spencer (1971)
    eccentricity = 1.000110_real64 &
      + 0.034221_real64 * cos(b) + 0.001280_real64 * sin(b) &
      + 0.000719_real64 * cos(2.0_real64*b) + 0.000077_real64 * sin(2.0_real64*b)
    elev = solar_elevation_deg(lat_deg, lon_deg, doy, hour_jst)
    if (elev <= 0.0_real64) then
      extraterrestrial_radiation_wm2 = 0.0_real64
    else
      extraterrestrial_radiation_wm2 = SOLAR_CONSTANT * eccentricity * sin(elev * DEG2RAD)
    end if
  end function

end module solar_geometry_mo

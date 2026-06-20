# fortran-solar

Physics-based solar-position functions for Fortran (Spencer, 1971).

- **Zero dependencies** — only `iso_fortran_env` (intrinsic). No external modules, no database/IO.
- Pure / elemental-friendly `real(real64)` functions.
- JST-based (standard meridian 135°E); angles in **degrees** unless suffixed `_rad`.

Shared by `tkd-05-solar` (日射 fill/forecast) and `tkd-pv-forecast` (太陽光発電 forecast).

## Public API

`module solar_geometry_mo`

| function | returns | description |
|---|---|---|
| `solar_elevation_deg(lat_deg, lon_deg, doy, hour_jst)` | `real(real64)` | 太陽高度角 (degrees) |
| `solar_azimuth_deg(lat_deg, lon_deg, doy, hour_jst)` | `real(real64)` | 太陽方位角 (degrees, 北から時計回り) |
| `extraterrestrial_radiation_wm2(lat_deg, lon_deg, doy, hour_jst)` | `real(real64)` | 大気外水平面日射量 (W/m²); 太陽が地平線以下なら 0 |
| `day_of_year(year, month, day)` | `integer` | 年間通日 (1–366, 閏年対応) |

Arguments: `lat_deg`/`lon_deg` = 緯度・経度 (東経, degrees), `doy` = 年間通日,
`hour_jst` = JST 時刻 (0–23, 小数可).

## Usage

```fortran
use solar_geometry_mo
real(real64) :: elev, az, i0
integer :: doy

doy  = day_of_year(2026, 6, 21)                 ! 夏至
elev = solar_elevation_deg(35.69_8, 139.69_8, doy, 12.0_8)   ! Tokyo, 正午
az   = solar_azimuth_deg  (35.69_8, 139.69_8, doy, 12.0_8)
i0   = extraterrestrial_radiation_wm2(35.69_8, 139.69_8, doy, 12.0_8)
```

## Consume it

Symlink the module into your project's `src/` (as the tkd-* servers do with the
other `fortran-*` tools), or `curl` it in a Dockerfile:

```
https://raw.githubusercontent.com/tkdhss111/fortran-solar/main/src/solar_geometry_mo.f90
```

## Test

```
make test          # cmake + ninja + ctest
```

## Reference

Spencer, J. W. (1971). *Fourier series representation of the position of the sun.*

## License

MIT © Hisashi Takeda

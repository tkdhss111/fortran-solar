# fortran-solar

太陽位置・大気外日射量を計算する **依存ゼロ** の Fortran モジュール（Spencer, 1971）。

`tkd-05-solar`（日射 `rj` の補完・予測）と、今後構築する `tkd-pv-forecast`
（太陽光発電予測）の **共有ライブラリ**。物理ベースの純粋関数のみを提供する。

## 特徴

- **依存ゼロ** — 組込モジュール `iso_fortran_env` のみを使用。外部モジュール・
  データベース・ファイル I/O は一切なし。`#include` も不要（`PI` はモジュール内に定義）。
- すべて `real(real64)` の `pure` 関数。副作用なし・スレッド安全。
- **JST 基準**（標準子午線 東経 135°）。角度は特記なき限り **度（degree）**、
  接尾辞 `_rad` はラジアン。
- 太陽赤緯・均時差・地球–太陽距離補正は Spencer (1971) のフーリエ級数近似。

## 公開 API（`module solar_geometry_mo`）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `solar_elevation_deg(lat_deg, lon_deg, doy, hour_jst)` | `real(real64)` | 太陽高度角（度） |
| `solar_azimuth_deg(lat_deg, lon_deg, doy, hour_jst)` | `real(real64)` | 太陽方位角（度、**北から時計回り**。東=90, 南=180, 西=270） |
| `extraterrestrial_radiation_wm2(lat_deg, lon_deg, doy, hour_jst)` | `real(real64)` | 大気外水平面日射量（W/m²）。太陽が地平線以下なら 0 |
| `day_of_year(year, month, day)` | `integer` | 年間通日（1–366、閏年対応） |

### 引数

| 引数 | 型 | 説明 |
|---|---|---|
| `lat_deg`  | `real(real64)` | 緯度（北緯、度） |
| `lon_deg`  | `real(real64)` | 経度（東経、度） |
| `doy`      | `integer`      | 年間通日（`day_of_year` で算出） |
| `hour_jst` | `real(real64)` | JST 時刻（0–23、小数可。例：`12.5` = 12:30） |

## 使い方

```fortran
use solar_geometry_mo
real(real64) :: elev, az, i0
integer      :: doy

doy  = day_of_year(2026, 6, 21)                                ! 夏至
elev = solar_elevation_deg(35.69_8, 139.69_8, doy, 12.0_8)     ! 東京・正午の高度角
az   = solar_azimuth_deg  (35.69_8, 139.69_8, doy, 12.0_8)     ! 同・方位角
i0   = extraterrestrial_radiation_wm2(35.69_8, 139.69_8, doy, 12.0_8)
```

## 物理モデル（Spencer, 1971）

- **日角** `B = 2π(doy-1)/365`
- **太陽赤緯 δ**・**均時差 EoT**：B のフーリエ級数近似
- **太陽時** `= JST + EoT/60 + (経度 - 135)×4/60`（135° は JST 標準子午線）
- **時角** `ω = (太陽時 - 12)×15°`
- **高度角** `sin(h) = sinφ·sinδ + cosφ·cosδ·cosω`
- **方位角**：高度角・赤緯・緯度から算出（北基準・時計回り）
- **大気外日射** `I0 = 1361 × 離心率補正 × sin(h)`（太陽定数 1361 W/m²、`h ≤ 0` で 0）

## 「正午」と方位角についての注意

`hour_jst = 12.0` は **太陽南中時刻と一致しない**。南中時刻は経度と均時差で決まり、
例えば東京（東経 139.69°）では概ね **11:41 JST** 頃。したがって 12:00 JST の太陽は
南よりやや **西** に位置し、方位角は約 **198°**（180° ではない）。これは仕様どおりの
正しい挙動である（「正午＝真南」を前提にしたテストは誤り）。真南が必要なら経度 135°
かつ均時差 0 の条件で評価すること。

## 取り込み方

各 `tkd-*` サーバと同様、本モジュールを `src/` にシンボリックリンクするか、
Dockerfile で GitHub から取得する。

```bash
# ローカル開発（シンボリックリンク）
ln -s /home/hss/2_tools/fortran-solar/src/solar_geometry_mo.f90 \
      <your-project>/src/solar_geometry_mo.f90
```

```dockerfile
# コンテナビルド（他の fortran-* 依存と同じパターン）
ARG FORTRAN_SOLAR_REF=main
RUN test -e src/solar_geometry_mo.f90 || \
    curl -fsSL https://raw.githubusercontent.com/tkdhss111/fortran-solar/${FORTRAN_SOLAR_REF}/src/solar_geometry_mo.f90 \
      -o src/solar_geometry_mo.f90
```

## テスト

```bash
make test          # cmake + ninja + ctest（太陽位置のサニティチェック）
```

## 参考文献

Spencer, J. W. (1971). *Fourier series representation of the position of the sun.*
Search 2(5), 172.

## ライセンス

MIT © Hisashi Takeda

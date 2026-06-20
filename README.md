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
| `incidence_angle_deg(lat_deg, lon_deg, doy, hour_jst, tilt_deg, panel_azimuth_deg)` | `real(real64)` | **入射角**（度）。任意方位・傾斜面と太陽光線のなす角。傾斜面日射（transposition）の基本量。下記参照 |
| `optimal_tilt_deg(lat_deg, lon_deg [, panel_azimuth_deg] [, diffuse_fraction])` | `real(real64)` | **最適傾斜角**（度）。指定方位の固定パネルで年間日射量を最大化する傾斜角。下記参照 |

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
real(real64) :: elev, az, i0, tilt
integer      :: doy

doy  = day_of_year( 2026, 6, 21 )                                 ! 夏至
elev = solar_elevation_deg( 35.69_8, 139.69_8, doy, 12.0_8 )      ! 東京・正午の高度角
az   = solar_azimuth_deg( 35.69_8, 139.69_8, doy, 12.0_8 )        ! 同・方位角
i0   = extraterrestrial_radiation_wm2( 35.69_8, 139.69_8, doy, 12.0_8 )
tilt = optimal_tilt_deg( 35.69_8, 139.69_8 )                      ! 東京の最適傾斜角（度）
```

> コード規約: ユーザ定義手続きの呼出し/定義は括弧内に空白 `f( a1, a2 )` を入れ、
> 配列添字 `arr(i,j)`・型種別 `real(real64)`・組込関数 `sin(x)` と区別する。

## 物理モデル（Spencer, 1971）

- **日角** `B = 2π(doy-1)/365`
- **太陽赤緯 δ**・**均時差 EoT**：B のフーリエ級数近似
- **太陽時** `= JST + EoT/60 + (経度 - 135)×4/60`（135° は JST 標準子午線）
- **時角** `ω = (太陽時 - 12)×15°`
- **高度角** `sin(h) = sinφ·sinδ + cosφ·cosδ·cosω`
- **方位角**：高度角・赤緯・緯度から算出（北基準・時計回り）
- **大気外日射** `I0 = 1361 × 離心率補正 × sin(h)`（太陽定数 1361 W/m²、`h ≤ 0` で 0）

## 最適傾斜角（`optimal_tilt_deg`）

南向きに固定設置した太陽光パネルが **年間日射量を最大化** する傾斜角（度）を返す。
PV 発電量予測（tkd-pv-forecast）で、設置地点ごとのパネル角を決めるために用いる。

### 計算手順
1. 傾斜角 β を 0°〜90° まで 1° 刻みで走査する。
2. 各 β について 1 年間（365 日 × 毎時）の **傾斜面日射量 POA** を積算する。
3. POA が最大となる β を最適傾斜角として返す。

### 日射モデル（等方天空 / isotropic sky）

傾斜面日射 = **直達** + **散乱** の和:

```
POA = (1 - fd)·GHI·Rb           直達（傾斜面）   Rb = cosθ_tilt / sin(h)
    +      fd ·GHI·(1+cosβ)/2   散乱（等方天空）  (1+cosβ)/2 = 天空率
```

- `GHI`（全天日射）の代理として **大気外水平面日射** を用いる（晴天指数 Kt は一様と
  仮定し、最大化では相殺するため省略）。
- `cosθ_tilt`：南向き傾斜面への直達入射角余弦（等価緯度 φ−β で評価。背面側は 0）。
- `fd`：散乱日射比率（散乱/全天）。`diffuse_fraction` 引数（既定 0.5）。

### 引数

| 引数 | 既定 | 説明 |
|---|---|---|
| `lat_deg` | — | 緯度（北緯, 度） |
| `lon_deg` | — | 経度（東経, 度） |
| `panel_azimuth_deg` | 180（真南） | 任意。パネル方位（北=0, 東=90, 南=180, 西=270） |
| `diffuse_fraction` | 0.5 | 任意。年間散乱日射比率。大きいほど最適傾斜角は浅くなる |

> 任意引数は keyword で指定する: `optimal_tilt_deg( lat, lon, panel_azimuth_deg=135.0_8 )`。
> 例（東京）: 真南→24°、南東 135°→21°（南より浅い）、東 90°≈西 270°（対称）。

### 結果の目安（既定 fd = 0.5）

| 地点 | 緯度 | 最適傾斜角 |
|---|---|---|
| 那覇 | 26.2° | 17° |
| 東京 | 35.7° | 24° |
| 稚内 | 45.4° | 31° |

- **緯度が高いほど急**（冬の低い太陽を捉えるため）。
- **散乱比率 fd を上げると浅く**（散乱光は水平面で多く受かる）。例: 東京 fd=0.8 → 12°。

### 注意・限界

- **等方天空モデルの近似**：実際の空は太陽周辺光・地平輝度（異方性）があり、最適傾斜角は
  本モデルよりやや急になる傾向。経験則「傾斜角 ≒ 緯度」や異方性モデル（Perez 等）は
  東京で ~30° を与えるのに対し、本モデルは既定 fd=0.5 で ~24°（数度浅め）。
- **fd で校正可能**：地域の散乱日射比率に合わせて `diffuse_fraction` を調整する。fd を
  小さく（~0.3）すると緯度に近い急めの値になる。
- **南向き固定のみ**（方位角は真南を仮定）。東西向き屋根や追尾架台は対象外。
- 計算量は約 91（傾斜角）× 8760（時間）回の太陽位置評価。`pure` 関数で副作用なく
  数十 ms で完了する。

## 入射角と傾斜面日射（transposition）

`incidence_angle_deg( lat_deg, lon_deg, doy, hour_jst, tilt_deg, panel_azimuth_deg )`
は任意方位・傾斜面と太陽光線のなす **入射角 θ**（度）を返す。PV の傾斜面日射
（transposition: 水平面日射 → 傾斜面日射）モデルの基本量である。

```
cosθ = sin(h)·cos(β) + cos(h)·sin(β)·cos(γs - γp)
       h=太陽高度角, β=傾斜角, γs=太陽方位, γp=パネル方位
```

### GHI → POA パイプライン（実装済み・既定 Engerer2-30min + Perez）

予報精度（＝当社の商品価値）のため、傾斜面全天日射 POA を次式で構成する:

```
POA = DNI·cosθ              直達（incidence_angle_deg が基盤）
    + DHI·R_d(model)        天空散乱
    + GHI·ρ·(1 - cosβ)/2    地面反射（アルベド ρ）
```

天空散乱係数 `R_d` の世界標準モデル（精度の低い順）:

| モデル | 内容 | 精度 |
|---|---|---|
| 等方（Liu–Jordan） | `(1+cosβ)/2` のみ | 基準 |
| Hay–Davies | 太陽周辺光（circumsolar）を異方的に考慮 | 良 |
| Reindl | Hay–Davies + 地平輝度 | 良＋ |
| **Perez** | circumsolar + 地平輝度を係数ビンで精緻化 | **最良（PV 業界標準）** |

LFM GPV は **GHI のみ**（`solar_radiation_wm2`, 30分平均）なので、GHI を分離してから
変換する。実装関数（`solar_geometry_mo`、依存ゼロ）:

| 関数 | 説明 |
|---|---|
| `decompose_engerer2( ghi, lat, lon, doy, hour_jst, dni, dhi [, interval_min] )` | **Engerer2(30分版係数)** 直散分離。GHI→DNI,DHI |
| `decompose_erbs( … )` | Erbs(1982) 直散分離（A/B 比較用） |
| `poa_perez( dni, dhi, ghi, lat, lon, doy, hour_jst, tilt, panel_azimuth [, albedo] [, interval_min] )` | **Perez(1990)** 傾斜面 POA |
| `poa_from_ghi( ghi, lat, lon, doy, hour_jst, tilt, panel_azimuth [, albedo] [, interval_min] [, use_erbs] )` | GHI→POA 一括（既定 Engerer2+Perez。`use_erbs=.true.` で Erbs） |
| `clear_sky_ghi_wm2` / `apparent_solar_time_h` / `air_mass_kastenyoung` | 補助量 |

```fortran
! LFM: solar_radiation_wm2 = GHI[W/m2]。jst から doy, hour_jst を算出して:
poa = poa_from_ghi( ghi, lat, lon, doy, hour_jst, tilt_deg=30.0_8, panel_azimuth_deg=180.0_8 )
```

**精度に関する注意**:
- `interval_min`（既定 **30**）で 30 分平均 GHI に整合（幾何量は区間中点で評価）。瞬時データは 0。
- 晴天モデルは **TJ (ASHRAE/Threlkeld–Jordan)** を内蔵（Engerer2 30分版の元適合と同一。
  A,k,C は doy の正弦関数）。`ktc`/`kde` の整合が取れる。
- A/B 比較ハーネス: `examples/ab_decomposition/`（LFM GHI に対し Engerer2 vs Erbs の
  POA 差を集計）。実データ 51,862 点では平均 |差| ≈ 0.85%、最大の乖離は中清明度
  （GHI 500–800）帯に現れる。
- Perez 係数は **allsitescomposite1990**（ソースに全値明示。要時 pvlib と照合可）。
- 検証済み: 水平面で `poa_perez ≡ GHI`（厳密）、分離は `DHI + DNI·sin h ≡ GHI`（エネルギー保存）。

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

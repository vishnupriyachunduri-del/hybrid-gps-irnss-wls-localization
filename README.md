# Hybrid GPS–IRNSS Positioning and Accuracy Enhancement using Weighted Least Squares
### Advanced Communication Laboratory Project | ECE Department

---

## Overview

This MATLAB simulation implements a **Hybrid GPS–IRNSS positioning system** that fuses signals from both GPS (US) and IRNSS/NavIC (India) satellite constellations to compute receiver position with high accuracy using **Weighted Least Squares (WLS)**.

The upgraded version adds **ionospheric correction (Klobuchar model)**, **tropospheric correction (Saastamoinen model)**, **RAIM integrity monitoring**, **PDOP/HDOP/VDOP analysis**, and ECEF-to-geodetic coordinate conversion — making it a near-realistic simulation of a dual-constellation GNSS receiver.

---

## Why Hybrid GPS + IRNSS?

| Feature | GPS Only | IRNSS Only | Hybrid GPS+IRNSS |
|---|---|---|---|
| Coverage | Global | South Asia | Global + enhanced India |
| No. of Satellites | ~31 | 7 | 11+ combined |
| Positioning Accuracy | ~5–10 m | ~10–20 m | **~3–6 m (improved)** |
| PDOP (typical) | ~2.5 | ~3.5 | **~1.5–2.0** |
| Integrity (RAIM) | Possible | Limited | **Robust** |

Combining both constellations increases the number of visible satellites, which improves **DOP geometry**, reduces **positioning error**, and enables **RAIM fault detection**.

---

## Project Structure

```
Hybrid_GPS_IRNSS_WLS.m       ← Main MATLAB simulation file
README.md                     ← This file
outputs/
  fig1_position_error.png     ← 3D position error per epoch
  fig2_dop_values.png         ← PDOP / HDOP / VDOP per epoch
  fig3_accuracy_summary.png   ← Bar chart of accuracy metrics
  fig4_cdf_error.png          ← Cumulative distribution of error
  fig5_ecef_components.png    ← X / Y / Z error components
  fig6_atmospheric_corrections.png  ← Iono & tropo delays per satellite
  fig7_command_window_output.png    ← Expected terminal output
```

---

## Upgrades Over Base Version

| Feature | Base Version | Upgraded Version |
|---|---|---|
| Satellites | 8 random | 6 GPS + 5 IRNSS (realistic geometry) |
| Epochs | 10 | 20 |
| Receiver position | Random ECEF | Hyderabad, India (real coordinates) |
| Ionospheric correction | ❌ | ✅ Klobuchar model |
| Tropospheric correction | ❌ | ✅ Saastamoinen model |
| Weighting | SNR only | SNR + Elevation-based |
| DOP computation | ❌ | ✅ PDOP, HDOP, VDOP per epoch |
| RAIM | ❌ | ✅ Chi-squared residual test |
| Coordinate output | ECEF only | ✅ ECEF + Lat/Lon/Alt |
| Accuracy metrics | Mean, RMS | ✅ + Max, Std, CEP50, CEP95 |
| Plots | 2 | 6 (error, DOP, bar, CDF, components, atmospheric) |
| Reproducibility | Random seed | ✅ Fixed seed (rng(42)) |

---

## How It Works — Step by Step

### Step 1: Satellite Geometry Setup
Satellite positions are placed in ECEF coordinates at realistic orbital altitudes:
- **GPS**: ~26,200 km (Medium Earth Orbit)
- **IRNSS**: ~36,000 km (Geostationary/Geosynchronous Orbit)

The true receiver position is set to **Hyderabad, India** (Lat: 17.385°N, Lon: 78.487°E, Alt: 531 m), converted from geodetic to ECEF using the WGS-84 ellipsoid.

### Step 2: Atmospheric Delay Correction

**Klobuchar Ionospheric Model (GPS L1)**  
Models the delay experienced by the GPS signal travelling through the ionosphere. Uses broadcast coefficients (alpha, beta) and accounts for local time and receiver latitude.

```
Iono delay = F × [5ns + A × (1 − x²/2 + x⁴/24)]   [for |x| < π/2]
```
where F is the obliquity factor and A is the amplitude determined from alpha coefficients.

**Saastamoinen Tropospheric Model**  
Models signal delay through the neutral atmosphere (troposphere). Uses standard pressure (1013.25 hPa), temperature (288.15 K), and relative humidity:

```
Trop delay (zenith) = 0.002277 × (P + (1255/T + 0.05) × RH × e_s)
Trop delay (slant)  = Zenith delay / sin(elevation)
```

### Step 3: Pseudorange Generation
For each satellite–epoch pair:
```
PR = true_range + c × clock_error + iono_delay + trop_delay + noise
```
Noise is Gaussian with σ = 3 m (reduced from 5 m in the base version).

### Step 4: Weighted Least Squares (WLS) Positioning
The WLS algorithm iteratively solves for 4 unknowns: receiver X, Y, Z, and clock bias (dt).

**Observation model:**
```
PR_i = ρ_i + c·dt + ε_i
```

**Design matrix H (linearised):**
```
H_i = [-(sat_x - x)/ρ,  -(sat_y - y)/ρ,  -(sat_z - z)/ρ,  1]
```

**WLS update:**
```
Δx = (Hᵀ W H)⁻¹ Hᵀ W y
```

**Weight matrix W** combines SNR and elevation:
```
w_i = SNR_weight_i × sin²(elevation_i)
```
Higher elevation and higher SNR = higher weight = less noise influence.

Iteration continues until position update `||Δx|| < 0.001 m`.

### Step 5: RAIM (Receiver Autonomous Integrity Monitoring)
After solving, residuals are computed and a test statistic is calculated:
```
test_stat = Σ(residuals²) / (m − 4)
```
If `test_stat > 30`, a RAIM alarm is raised for that epoch, flagging possible satellite fault.

### Step 6: DOP Computation
From the geometry matrix Q = (HᵀH)⁻¹:
- **PDOP** = √(Q₁₁ + Q₂₂ + Q₃₃)
- **HDOP** = √(Q_ENU₁₁ + Q_ENU₂₂)  ← horizontal (ENU frame)
- **VDOP** = √(Q_ENU₃₃)              ← vertical

Lower DOP = better satellite geometry = more accurate position.

---

## Expected Output (Command Window)

```
=============================================================
    Hybrid GPS-IRNSS Positioning Simulation (Upgraded)
=============================================================

[1/6] Initialising simulation parameters...
      GPS satellites  : 6
      IRNSS satellites: 5
      Total epochs    : 20

[3/6] Applying atmospheric correction models...
      Ionospheric corrections  (Klobuchar): applied.
      Tropospheric corrections (Saastamoinen): applied.

[5/6] Running WLS positioning with RAIM integrity monitoring...
      WLS positioning complete.

=============================================================
                    ACCURACY RESULTS
=============================================================
  Mean 3D Position Error  :   X.XXXX  m
  RMS  3D Position Error  :   X.XXXX  m
  Max  3D Position Error  :   X.XXXX  m
  CEP 50%  (2D)           :   X.XXXX  m
  CEP 95%  (2D)           :   X.XXXX  m
  Mean PDOP               :   X.XXXX
  Mean HDOP               :   X.XXXX
  Mean VDOP               :   X.XXXX
  RAIM Alarms             :   0 / 20 epochs
=============================================================

=== Simulation Completed Successfully ===
```

---

## Software Requirements

- **MATLAB R2019b or later** (or GNU Octave 6.0+)
- No additional toolboxes required
- All code uses base MATLAB functions only

---

## How to Run

1. Open MATLAB
2. Navigate to the folder containing `Hybrid_GPS_IRNSS_WLS.m`
3. Type in the Command Window:
   ```matlab
   run('Hybrid_GPS_IRNSS_WLS.m')
   ```
   or simply open the file and press **Run (F5)**
4. Six figure windows will open automatically

---

## Key Concepts Covered

- GNSS signal propagation and pseudorange measurement
- WGS-84 coordinate system and ECEF-to-geodetic conversion
- Weighted Least Squares estimation
- Klobuchar ionospheric delay model
- Saastamoinen tropospheric delay model
- Dilution of Precision (PDOP, HDOP, VDOP)
- RAIM integrity monitoring
- Satellite constellation geometry analysis
- Statistical accuracy metrics (CEP, RMS, CDF)

---

## Authors

**Team:** [Chunduri Vishnu Priya] | ECE Department, MGIT Hyderabad  
**Course:** Advanced Communication Laboratory  
**Institution:** Mahatma Gandhi Institute of Technology, Hyderabad

---

*This is a simulation. Real GNSS receivers additionally handle signal acquisition, carrier phase, multipath, and hardware biases.*

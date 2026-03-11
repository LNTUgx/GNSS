# EGZTD

This repository provides the EGZTD model, a high-precision empirical Zenith Tropospheric Delay (ZTD) model derived from ERA5 data. The repository includes necessary coefficient matrices files and computation scripts (MATLAB) that allow users to instantly resolve ZTD for any given spatiotemporal information.

## 📌 Model Overview

- **Spatial Resolution**:  
  - Longitude: 1.0°  
  - Latitude: 1.0°  

## 🔧 Files and Descriptions

| File/Folder            | Description                                                  |
| ---------------------- | ------------------------------------------------------------ |
| `latlon_dh.txt`      | Latitude, longitude, and geoid undulation for the 1°×1° sliding grids of the EGZTD model.      |
| `surface_height.mat`         | Ellipsoidal heights of all grid points. |
| `ZTD_7pre`     | Annual mean (Col. 1), annual amplitudes (Col. 2-3), semi-annual amplitudes (Col. 4-5), and diurnal amplitudes (Col. 6-7) of ZTD. |
| `ctx1` | Annual mean (Col. 1), annual amplitudes (Col. 2-3), and semi-annual amplitudes (Col. 4-5) of first-order coefficient for the exponential function. |
| `ctx2` | Annual mean (Col. 1), annual amplitudes (Col. 2-3), and semi-annual amplitudes (Col. 4-5) of second-order coefficient for the exponential function. |
| `ptx1` | Annual mean (Col. 1), annual amplitudes (Col. 2-3), and semi-annual amplitudes (Col. 4-5) of first-order coefficient for the quadratic polynomial. |
| `ptx2` | Annual mean (Col. 1), annual amplitudes (Col. 2-3), and semi-annual amplitudes (Col. 4-5) of second-order coefficient for the quadratic polynomial. |
| `EGZTD.m` | The core function file of EGZTD that invokes the coefficient matrix for ZTD calculation. |
| `main.m` | The main function for ZTD calculation using user-input spatiotemporal information. |

## 📖 Reference

The following literature is available for users’ reference:

- Jiang C, Gao X, Zhu H, et al. An improved global pressure and zenith wet delay model with optimized vertical correction considering the spatiotemporal variability in multiple height-scale factors. Geoscientific Model Development, 2024, 17(15): 5939-5959.
- Gao X, Li P, Chen S, et al. A stacking ensemble learning approach for enhancing global real-time zenith total delay modeling using global forecast system forecasts. Measurement Science and Technology, 2025, 36(11): 116311.

## 📬 Contact

For questions or collaboration, please contact:

**Xiang Gao**  
📧 xianggao@chd.edu.cn

---

Thank you for using EGZTD!


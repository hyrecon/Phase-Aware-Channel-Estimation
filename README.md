# 📡 Phase-Aware Channel Estimation for Monostatic Backscatter Communications

> Lifted LS/LMMSE estimators with closed-form phase correction for ultra-low-power monostatic backscatter links.

## 📌 About

Conventional channel estimators for backscatter systems assume that the residual carrier frequency offset (CFO) after coarse compensation is negligible. In practice, however, even a small residual CFO introduces a per-sample phase drift that accumulates across the frame — causing the standard LS/LMMSE estimators to exhibit irreducible error floors regardless of SNR.

This project implements **phase-aware lifted channel estimators** that explicitly model the residual phase slope as an additional unknown parameter. By constructing an augmented observation model with pilot-based sufficient statistics, both the composite channel gain and the residual phase drift are jointly resolved in closed form — without iterative optimization or grid search.

## 🎯 Key Results

- NMSE suppression: up to **27.8 dB** improvement over conventional LS at moderate residual CFO
- Complexity: **O(N)** per frame — no matrix inversion required for LS mode
- Modulation support: OOK and BPSK
- BER: < 10⁻⁴ at SNR = 20 dB with residual phase drift
- Estimator modes: **LS** (no prior) and **LMMSE** (with Bayesian prior)
- Validated via over-the-air image transmission on a monostatic backscatter platform

## 🔧 System Architecture

**Algorithm 1 — Phase-Aware Lifted Channel Estimation:**
Accumulates sufficient statistics (s₀, s₁, s₂, t₀, t₁) over the pilot set, then computes the lifted parameter vector [θ̂₀, θ̂₁] via closed-form LS or LMMSE solution. The composite channel gain ĥ and residual phase slope φ̂₁ are extracted from the ratio of the estimated parameters.

**Algorithm 2 — Monostatic Backscatter Receiver Processing:**
Full baseband receiver pipeline: DC offset removal → CFO estimation (autocorrelation-based) → CFO compensation → FFT-based low-pass filtering → per-frame phase-aware channel estimation (Algorithm 1) → single-tap equalization → OOK/BPSK detection.

**Platform:**
Semi-passive tag (SPDT switch) + USRP N200 SDR reader + GNU Radio-based transceiver + MATLAB post-processing.

## 🛠️ Built With

`MATLAB` · `Python` · `SDR`

USRP N200 · semi-passive backscatter tag · planar Yagi-Uda antenna array · microwave absorber

## 📄 Publication

**H. Ryu** and S. Kim, "Ultra-low-power Monostatic Backscatter Platform with Phase-Aware Channel Estimation and System-Level Validation," [IEEE Internet of Things Journal]([https://arxiv.org/abs/2601.02227](https://ieeexplore.ieee.org/document/11535284)), May. 2026.

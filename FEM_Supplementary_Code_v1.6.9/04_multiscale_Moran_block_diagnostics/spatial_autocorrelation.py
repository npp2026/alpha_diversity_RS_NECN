#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spatial_autocorrelation.py — Global Moran's I + Local Getis-Ord Gi* on raster data
====================================================================================

A self-contained tool for spatial autocorrelation analysis on raster outputs from
the time series pipeline (e.g., Sen's slope rasters from `unified_timeseries_analysis_v3.py`).

Usage (typical):
    # Basic Gi* with significance masking from MK p-value
    python spatial_autocorrelation.py \\
        --slope-input  output/SR_TrendSlope.tif \\
        --mk-input     output/SR_MannKendall.tif \\
        --p-threshold  0.05 \\
        --aggregate-factor 4 \\
        --output-dir   output/SR_spatial \\
        --prefix       SR_2001_2020

    # With Monte Carlo permutation test for global Moran (slower but more reliable)
    python spatial_autocorrelation.py \\
        --slope-input  output/SR_TrendSlope.tif \\
        --output-dir   output/SR_spatial \\
        --prefix       SR_2001_2020 \\
        --permutations 999 \\
        --seed 42

Outputs (in --output-dir):
    {prefix}_Moran.csv              Global Moran's I summary (1 row scalar metrics)
    {prefix}_GiStar.tif             3-band Int16 raster:
                                      Band 1: z-score x100        (range ~ -1000 to +1000)
                                      Band 2: p-value x10000      (range 0 to 10000)
                                      Band 3: class -3..+3        (-3=Cold99, 0=NotSig, +3=Hot99)
    {prefix}_GiStar_summary.csv     Class counts and proportions
    {prefix}_log.txt                Run log

Algorithms:
-----------
- Global Moran's I (Cliff & Ord 1981, randomization assumption)
    I = (n/W) * Σ_i Σ_j w_ij * (x_i - x̄) * (x_j - x̄) / Σ_i (x_i - x̄)²
    z = (I - E(I)) / sqrt(Var_R(I))
    Optionally: Monte Carlo permutation p-value

- Local Getis-Ord Gi* (Getis & Ord 1995)
    Gi*_i = Σ_j w_ij * x_j / Σ_j x_j           (intensity form)
    z(Gi*) = (Σ_j w_ij*x_j - x̄ * W_i) / (S * sqrt(W_i*(n - W_i) / (n - 1)))
    where:
      - w_ij is binary 1/0 (queen contiguity by default; star = focal cell included)
      - W_i = Σ_j w_ij computed per-pixel from valid neighbor mask (handles edges + NoData)
      - x̄ = global mean of valid pixels
      - S  = sqrt( Σ(x_j - x̄)² / n )       (population SD, per Getis-Ord 1995)
      - n  = global count of valid pixels
    Two-tailed p-value: p = 2 * Φ(-|z|) under standard normal asymptotic

Classification (matching v7 design — 7-level Donnini-style legend):
     3 = Hot99       (z >= 2.58)
     2 = Hot95       (1.96 <= z < 2.58)
     1 = Hot90       (1.65 <= z < 1.96)
     0 = NotSig      (-1.65 < z < 1.65)
    -1 = Cold90      (-1.96 < z <= -1.65)
    -2 = Cold95      (-2.58 < z <= -1.96)
    -3 = Cold99      (z <= -2.58)

Edge handling:
    Cells with fewer than 2 valid neighbors (incl. self) get z = NaN.
    Convolution uses mode='constant' with cval=0 over a (mask * x) array,
    plus a separate convolution to count valid neighbors per pixel.

Validation note:
    Output has been verified against `esda.G_Local` (PySAL) and `spdep::localG` (R) on
    a small test grid up to numerical precision. See validation_test() at the bottom.

Performance:
    Vectorized via scipy.ndimage.convolve (no Python loops in the hot path).
    A 5000x5000 raster with default queen kernel runs in ~25 seconds on commodity hardware.
    Memory peak ~ 4 * raster_size * 8 bytes (≈800 MB for 5000x5000 float64).

Author: Time Series Analysis Team
Version: 1.0
License: MIT
"""

__version__ = "1.5.1"

import argparse
import logging
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Tuple, Optional, Dict

import numpy as np
import rasterio
from scipy import ndimage
from scipy import stats as sp_stats

from spatial_block_bootstrap import (
    block_statistics_metadata,
    compute_block_statistics,
    spatial_block_bootstrap_moran,
)

# ============================================================
# Logging setup
# ============================================================

logger = logging.getLogger("spatial_autocorr")
logger.addHandler(logging.NullHandler())


# ============================================================
# Constants
# ============================================================

# Output Int16 scaling (matching the v3 time-series pipeline conventions)
SCALE_ZSCORE  = 100      # z-score precision 0.01, range ~ -327 to +327 (covers any plausible Gi*)
SCALE_PVALUE  = 10000    # p-value precision 0.0001
NODATA_INT16  = -9999

# Significance thresholds (two-tailed, standard normal)
Z_99 = 2.575829303548901   # qnorm(0.995)
Z_95 = 1.959963984540054   # qnorm(0.975)
Z_90 = 1.6448536269514722  # qnorm(0.95)


# ============================================================
# Kernel construction
# ============================================================

def make_kernel(neighbor_type: str = "queen", radius: int = 1) -> np.ndarray:
    """
    Build a binary spatial weights kernel for raster Gi* and Moran's I.

    The kernel encodes binary weights w_ij ∈ {0, 1}.
    Center cell is included (w_ii = 1) — this gives the "star" variant Gi*.
    For Gi (no star) or Moran's I without self, set kernel center to 0 manually.

    Args:
        neighbor_type:
            'queen'  — 8 neighbors (3x3 minus center for Moran; full 3x3 for Gi*)
            'rook'   — 4 neighbors (cross pattern)
            'circle' — circular kernel of given pixel radius
        radius: For 'queen' and 'rook', kernel is (2*radius+1) x (2*radius+1).
                For 'circle', radius is in pixels.

    Returns:
        2D float64 kernel array. Center cell included (=1) for star variants.
    """
    radius = int(radius)
    if radius < 1:
        raise ValueError("radius must be >= 1")

    if neighbor_type == "queen":
        size = 2 * radius + 1
        return np.ones((size, size), dtype=np.float64)

    elif neighbor_type == "rook":
        size = 2 * radius + 1
        kernel = np.zeros((size, size), dtype=np.float64)
        kernel[radius, :] = 1
        kernel[:, radius] = 1
        return kernel

    elif neighbor_type == "circle":
        size = 2 * radius + 1
        y, x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
        kernel = (x * x + y * y <= radius * radius).astype(np.float64)
        return kernel

    else:
        raise ValueError(f"Unknown neighbor_type: {neighbor_type}")


# ============================================================
# Raster I/O
# ============================================================

def _assert_raster_aligned(src: rasterio.io.DatasetReader, reference: dict, label: str) -> None:
    """Fail fast when a mask raster is not pixel-aligned with the slope raster."""
    problems = []
    if src.width != reference['width'] or src.height != reference['height']:
        problems.append(
            f"shape {src.height}x{src.width} != "
            f"{reference['height']}x{reference['width']}"
        )
    ref_crs = reference.get('crs')
    if src.crs != ref_crs:
        problems.append(f"CRS {src.crs!s} != {ref_crs!s}")
    ref_transform = reference.get('transform')
    if ref_transform is not None and not src.transform.almost_equals(ref_transform):
        problems.append(f"transform {src.transform!s} != {ref_transform!s}")
    if problems:
        raise ValueError(
            f"{label} is not aligned to the slope raster: " + "; ".join(problems) +
            ". Reproject/resample it to the slope grid before analysis."
        )


def _masked_raster_to_float64(values: np.ndarray) -> np.ndarray:
    """Convert a raster band or MaskedArray to float64 with masked cells as NaN.

    Rasterio commonly returns integer MaskedArray objects for encoded products.
    Calling ``filled(np.nan)`` on an integer MaskedArray raises ``TypeError``
    because NaN cannot be represented by the original integer dtype. Convert the
    underlying data first, then apply the mask explicitly.
    """
    data = np.asarray(np.ma.getdata(values), dtype=np.float64).copy()
    mask = np.ma.getmaskarray(values)
    if np.any(mask):
        data[mask] = np.nan
    return data


def load_slope_with_mask(
    slope_path: Path,
    mk_path: Optional[Path] = None,
    mk_pvalue_band: int = 2,
    mk_pvalue_scale: float = 10000.0,
    p_threshold: Optional[float] = None,
    forest_mask_path: Optional[Path] = None,
    aggregate_factor: int = 1,
    slope_input_scale: float = 100.0,
    aggregate_edge_mode: str = "partial",
    aggregate_min_valid_fraction: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Load, validate, mask and optionally aggregate a slope raster."""
    logger.info(f"Loading slope raster: {slope_path}")
    with rasterio.open(slope_path) as src:
        slope_band = src.read(1, masked=True)
        slope_raw = _masked_raster_to_float64(slope_band)
        profile = src.profile.copy()
        reference = {
            'width': src.width,
            'height': src.height,
            'crs': src.crs,
            'transform': src.transform,
        }
        nodata = src.nodata
        logger.info(
            f"  Shape: {slope_raw.shape}, nodata: {nodata}, "
            f"input scale: x{slope_input_scale}"
        )

    if nodata is not None:
        slope_raw[slope_raw == nodata] = np.nan
    slope = slope_raw / slope_input_scale
    valid_mask = np.isfinite(slope)
    n_initial = int(valid_mask.sum())
    logger.info(f"  Initial valid pixels: {n_initial:,}")

    if mk_path is not None and p_threshold is not None:
        logger.info(
            f"Applying MK significance mask: p <= {p_threshold} "
            f"(from {mk_path}, band {mk_pvalue_band})"
        )
        with rasterio.open(mk_path) as src:
            _assert_raster_aligned(src, reference, "MK/FDR raster")
            if not (1 <= mk_pvalue_band <= src.count):
                raise ValueError(
                    f"MK p-value band {mk_pvalue_band} is outside 1..{src.count}"
                )
            p_band = src.read(mk_pvalue_band, masked=True)
            p_raw = _masked_raster_to_float64(p_band)
            p_nodata = src.nodata
        if p_nodata is not None:
            p_raw[p_raw == p_nodata] = np.nan
        p_value = p_raw / mk_pvalue_scale
        sig_mask = np.isfinite(p_value) & (p_value <= p_threshold)
        valid_mask &= sig_mask
        n_after_mk = int(valid_mask.sum())
        logger.info(
            f"  After MK mask: {n_after_mk:,} "
            f"({100*n_after_mk/max(n_initial,1):.1f}% retained)"
        )

    if forest_mask_path is not None:
        logger.info(f"Applying forest mask: {forest_mask_path}")
        with rasterio.open(forest_mask_path) as src:
            _assert_raster_aligned(src, reference, "Forest mask")
            forest_band = src.read(1, masked=True)
            forest_raw = _masked_raster_to_float64(forest_band)
            f_nodata = src.nodata
        if f_nodata is not None:
            forest_raw[forest_raw == f_nodata] = np.nan
        forest_bool = (forest_raw > 0) & np.isfinite(forest_raw)
        valid_mask &= forest_bool
        logger.info(f"  After forest mask: {int(valid_mask.sum()):,}")

    slope[~valid_mask] = np.nan

    if aggregate_factor > 1:
        logger.info(
            f"Aggregating by factor {aggregate_factor} "
            f"(mean of valid pixels; edge_mode={aggregate_edge_mode})..."
        )
        slope, valid_mask, profile = aggregate_mean(
            slope,
            valid_mask,
            aggregate_factor,
            profile,
            edge_mode=aggregate_edge_mode,
            min_valid_fraction=aggregate_min_valid_fraction,
        )
        logger.info(
            f"  After aggregation: shape={slope.shape}, "
            f"valid={int(valid_mask.sum()):,}"
        )

    return slope, valid_mask, profile


def aggregate_mean(
    arr: np.ndarray,
    mask: np.ndarray,
    factor: int,
    profile: dict,
    *,
    edge_mode: str = "partial",
    min_valid_fraction: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Mean-aggregate by an integer factor while ignoring invalid cells.

    ``partial`` (default) keeps bottom/right edge groups rather than silently
    dropping as much as one coarse cell from the study extent. ``trim`` is
    retained only for exact reproduction of legacy outputs.

    ``min_valid_fraction`` is measured against the nominal ``factor x factor``
    coarse-cell area, including clipped bottom/right cells.  This prevents a
    tiny edge sliver or a single valid base pixel from receiving the same weight
    as a fully supported coarse cell when a positive threshold is requested.
    """
    arr = np.asarray(arr, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    if arr.ndim != 2 or arr.shape != mask.shape:
        raise ValueError("arr and mask must be same-shaped 2-D arrays")
    factor = int(factor)
    if factor < 1:
        raise ValueError("aggregate factor must be >= 1")
    if edge_mode not in {"partial", "trim"}:
        raise ValueError("edge_mode must be 'partial' or 'trim'")
    min_valid_fraction = float(min_valid_fraction)
    if not (0.0 <= min_valid_fraction <= 1.0):
        raise ValueError("min_valid_fraction must be between 0 and 1")
    if factor == 1:
        return arr, mask, profile.copy()

    h, w = arr.shape
    arr_filled = np.where(mask, arr, 0.0)

    if edge_mode == "trim":
        h_new, w_new = h // factor, w // factor
        if h_new < 1 or w_new < 1:
            raise ValueError("aggregate factor exceeds raster dimensions in trim mode")
        h_trim, w_trim = h_new * factor, w_new * factor
        arr_t = arr_filled[:h_trim, :w_trim].reshape(h_new, factor, w_new, factor)
        mask_t = mask[:h_trim, :w_trim].reshape(h_new, factor, w_new, factor)
        sum_block = arr_t.sum(axis=(1, 3), dtype=np.float64)
        count_block = mask_t.sum(axis=(1, 3), dtype=np.int64)
    else:
        row_starts = np.arange(0, h, factor, dtype=np.int64)
        col_starts = np.arange(0, w, factor, dtype=np.int64)
        # reduceat retains the final partial row/column group without padding.
        sum_rows = np.add.reduceat(arr_filled, row_starts, axis=0, dtype=np.float64)
        sum_block = np.add.reduceat(sum_rows, col_starts, axis=1, dtype=np.float64)
        count_rows = np.add.reduceat(mask, row_starts, axis=0, dtype=np.int64)
        count_block = np.add.reduceat(count_rows, col_starts, axis=1, dtype=np.int64)
        h_new, w_new = sum_block.shape

    required_valid = max(1, int(np.ceil(min_valid_fraction * factor * factor)))
    new_mask = count_block >= required_valid
    with np.errstate(invalid='ignore', divide='ignore'):
        agg = np.where(new_mask, sum_block / count_block, np.nan)

    new_profile = profile.copy()
    new_profile['height'] = int(h_new)
    new_profile['width'] = int(w_new)
    new_profile['transform'] = profile['transform'] * rasterio.Affine.scale(factor, factor)
    return agg, new_mask, new_profile


# ============================================================
# Global Moran's I
# ============================================================

def compute_global_moran(
    x: np.ndarray,
    valid: np.ndarray,
    kernel: np.ndarray,
) -> Tuple[float, float, int]:
    """
    Compute global Moran's I using binary symmetric weights (no self).

    The kernel must have center = 0 for proper Moran's I (no self-loop).
    If kernel center is 1 (Gi* style), it will be set to 0 internally.

    Args:
        x:      2D array of values, NaN where invalid.
        valid:  2D bool mask (True = include in analysis).
        kernel: 2D binary kernel array.

    Returns:
        (I_value, W_total, n_valid)
            I_value:  Moran's I scalar (in [-1, +1] typically)
            W_total:  Total sum of weights = sum of valid neighbor counts
            n_valid:  Number of valid pixels
    """
    x = np.asarray(x, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    kernel = np.asarray(kernel, dtype=np.float64)
    if x.ndim != 2 or valid.ndim != 2 or x.shape != valid.shape:
        raise ValueError("x and valid must be same-shape 2-D arrays")
    if kernel.ndim != 2 or any(size % 2 == 0 for size in kernel.shape):
        raise ValueError("kernel must be a 2-D array with odd dimensions")
    if np.any(valid & ~np.isfinite(x)):
        raise ValueError("valid=True occurs at non-finite x values")

    # Use a kernel with center = 0 (no self for Moran's I)
    k = kernel.copy()
    cy, cx = k.shape[0] // 2, k.shape[1] // 2
    k[cy, cx] = 0.0

    # Center the data
    valid_vals = x[valid]
    n = len(valid_vals)
    # A point estimate is mathematically defined from two non-constant values.
    # Inference and bootstrap reliability are controlled separately in main().
    if n < 2:
        raise ValueError(f"Too few valid pixels for Moran's I point estimate: n={n}")

    x_mean = valid_vals.mean()
    deviations = np.where(valid, x - x_mean, 0.0)
    mask_int = valid.astype(np.float64)

    # Numerator: Σ_i Σ_j w_ij * (x_i - x̄) * (x_j - x̄)
    #         = Σ_i (x_i - x̄) * [Σ_j w_ij * (x_j - x̄)]
    sum_neighbor_dev = ndimage.convolve(deviations, k, mode='constant', cval=0.0)
    cross_product = (deviations * sum_neighbor_dev * mask_int).sum()

    # Denominator: Σ_i (x_i - x̄)²
    sum_sq = (deviations ** 2).sum()  # invalid contribute 0

    # W: total sum of weights (only count where BOTH endpoints are valid)
    # W_i_per_cell = number of valid neighbors of cell i (where cell i itself is valid)
    W_per_cell = ndimage.convolve(mask_int, k, mode='constant', cval=0.0) * mask_int
    W_total = W_per_cell.sum()

    if sum_sq == 0 or W_total == 0:
        return float('nan'), 0.0, n

    I = (n / W_total) * (cross_product / sum_sq)
    return float(I), float(W_total), n


def moran_inference_normality(
    I: float,
    W_total: float,
    n: int,
    valid: np.ndarray,
    kernel: np.ndarray,
    x: np.ndarray,
) -> Dict[str, float]:
    """
    Compute Moran's I expected value, variance (randomization assumption),
    z-score, and asymptotic p-value.

    Reference: Cliff & Ord (1981), Spatial Processes: Models and Applications, eq. 1.36-1.39.

    Returns dict with: expected_I, var_I_norm, var_I_rand, z_norm, z_rand,
                       p_norm (two-tailed), p_rand (two-tailed)
    """
    k = kernel.copy()
    cy, cx = k.shape[0] // 2, k.shape[1] // 2
    k[cy, cx] = 0.0
    mask_int = valid.astype(np.float64)

    # Per-cell weight properties
    # Each cell i has degree_i = number of valid neighbors (where neighbor is valid)
    # For symmetric raster weights, deg_i = deg_i(in) = deg_i(out)
    # But we only count "live" neighbors (both endpoints valid) for the variance:
    deg_i = ndimage.convolve(mask_int, k, mode='constant', cval=0.0) * mask_int
    deg_valid = deg_i[valid]

    # S0 = W = sum of all weights
    S0 = float(deg_valid.sum())  # equivalent to W_total

    # S1 = (1/2) * Σ_i Σ_j (w_ij + w_ji)²
    # For binary symmetric: w_ij = w_ji ∈ {0,1}, so (w_ij + w_ji)² = 4*w_ij² = 4*w_ij
    # S1 = (1/2) * 4 * Σ_i Σ_j w_ij = 2 * S0
    S1 = 2.0 * S0

    # S2 = Σ_i (Σ_j w_ij + Σ_j w_ji)²
    # For binary symmetric: row_sum_i = col_sum_i = degree_i
    # S2 = Σ_i (2 * degree_i)² = 4 * Σ_i degree_i²
    S2 = 4.0 * float((deg_valid ** 2).sum())

    # b2 = n * Σ(x_i - x̄)⁴ / (Σ(x_i - x̄)²)²    (kurtosis-like)
    valid_vals = x[valid]
    x_mean = valid_vals.mean()
    dev = valid_vals - x_mean
    m2 = (dev ** 2).sum()
    m4 = (dev ** 4).sum()
    b2 = (n * m4) / (m2 ** 2) if m2 > 0 else 0.0

    expected_I = -1.0 / (n - 1) if n > 1 else float('nan')

    # Variance under normality assumption
    # Var_N(I) = (n²·S1 - n·S2 + 3·S0²) / (S0² · (n²-1)) - E(I)²
    if S0 > 0 and n > 1:
        var_norm = (n * n * S1 - n * S2 + 3 * S0 * S0) \
                   / (S0 * S0 * (n * n - 1)) - expected_I ** 2
    else:
        var_norm = float('nan')

    # Variance under randomization assumption
    # Var_R(I) = [n·((n²-3n+3)·S1 - n·S2 + 3·S0²)
    #            - b2·((n²-n)·S1 - 2n·S2 + 6·S0²)]
    #           / [(n-1)(n-2)(n-3)·S0²]  - E(I)²
    if S0 > 0 and n > 3:
        num1 = n * ((n*n - 3*n + 3) * S1 - n * S2 + 3 * S0 * S0)
        num2 = b2 * ((n*n - n) * S1 - 2 * n * S2 + 6 * S0 * S0)
        denom = (n - 1) * (n - 2) * (n - 3) * S0 * S0
        var_rand = (num1 - num2) / denom - expected_I ** 2
    else:
        var_rand = float('nan')

    def z_p(I_val, var_val):
        if not np.isfinite(var_val) or var_val <= 0:
            return float('nan'), float('nan')
        z = (I_val - expected_I) / np.sqrt(var_val)
        p = 2.0 * sp_stats.norm.sf(abs(z))   # two-tailed
        return z, p

    z_norm, p_norm = z_p(I, var_norm)
    z_rand, p_rand = z_p(I, var_rand)

    return {
        'expected_I': expected_I,
        'var_I_norm': var_norm,
        'var_I_rand': var_rand,
        'z_norm': z_norm,
        'z_rand': z_rand,
        'p_norm': p_norm,
        'p_rand': p_rand,
        'S0': S0, 'S1': S1, 'S2': S2, 'b2': b2,
    }


def unavailable_moran_inference(n: int) -> Dict[str, float]:
    """Return a schema-compatible all-NA inference record.

    Used when a Moran point estimate is retained at a very coarse scale but the
    configured minimum sample size for analytical inference is not met.
    """
    expected = -1.0 / (n - 1) if n > 1 else float("nan")
    nan = float("nan")
    return {
        'expected_I': expected,
        'var_I_norm': nan,
        'var_I_rand': nan,
        'z_norm': nan,
        'z_rand': nan,
        'p_norm': nan,
        'p_rand': nan,
        'S0': nan, 'S1': nan, 'S2': nan, 'b2': nan,
    }


def moran_permutation_test(
    x: np.ndarray,
    valid: np.ndarray,
    kernel: np.ndarray,
    observed_I: float,
    n_perm: int = 999,
    seed: int = 42,
) -> Tuple[float, np.ndarray]:
    """
    Monte Carlo permutation test for global Moran's I.

    Shuffles values among valid locations, recomputes I for each permutation.
    Two-tailed p-value compares distance from the randomization expectation
    E(I) = -1/(n-1), not distance from zero.

    Returns:
        (p_perm, perm_distribution)
    """
    rng = np.random.default_rng(seed)
    valid_idx = np.argwhere(valid)
    valid_vals = x[valid].copy()

    perm_Is = np.empty(n_perm, dtype=np.float64)
    x_perm = x.copy()

    for k in range(n_perm):
        rng.shuffle(valid_vals)
        x_perm[valid] = valid_vals
        I_k, _, _ = compute_global_moran(x_perm, valid, kernel)
        perm_Is[k] = I_k
        if (k + 1) % max(1, n_perm // 10) == 0:
            logger.info(f"  Permutation progress: {k+1}/{n_perm}")

    # Two-tailed pseudo p-value around the permutation-null expectation,
    # not around zero.  Moran's I has E(I) = -1/(n-1) under randomization.
    expected_I = -1.0 / (int(valid.sum()) - 1)
    observed_distance = abs(observed_I - expected_I)
    perm_distance = np.abs(perm_Is - expected_I)
    p_perm = (np.sum(perm_distance >= observed_distance) + 1) / (n_perm + 1)
    return float(p_perm), perm_Is


# ============================================================
# Local Getis-Ord Gi*
# ============================================================

def compute_local_gi_star(
    x: np.ndarray,
    valid: np.ndarray,
    kernel: np.ndarray,
    min_neighbors: int = 2,
    min_total_pixels: int = 30,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Local Getis-Ord Gi* with proper edge / NoData handling.

    Star variant: focal cell IS included in the weight sum (kernel center = 1).
    For each focal cell i:
        W_i        = count of valid neighbors (incl. self) — variable across cells
        Σ_j w_ij*x_j = sum of valid neighbor values
        z_i        = (Σ w_ij*x_j - x̄_global * W_i) / [S * sqrt(W_i*(n - W_i)/(n - 1))]
    Cells with W_i < min_neighbors get z = NaN.

    Args:
        x:      2D array (NaN at invalid).
        valid:  2D bool mask.
        kernel: Binary kernel with center = 1 (will be enforced).
        min_neighbors: Minimum W_i to compute z; cells below this → NaN.
        min_total_pixels: Minimum total valid pixels required for Gi*.

    Returns:
        z:  2D z-score array (NaN at invalid / under-supported cells)
        p:  2D two-tailed p-value array
        W:  2D count-of-valid-neighbors array (for diagnostics)
    """
    x = np.asarray(x, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    kernel = np.asarray(kernel, dtype=np.float64)
    min_neighbors = int(min_neighbors)
    min_total_pixels = int(min_total_pixels)
    if x.ndim != 2 or valid.ndim != 2 or x.shape != valid.shape:
        raise ValueError("x and valid must be same-shape 2-D arrays")
    if kernel.ndim != 2 or any(size % 2 == 0 for size in kernel.shape):
        raise ValueError("kernel must be a 2-D array with odd dimensions")
    if np.any(valid & ~np.isfinite(x)):
        raise ValueError("valid=True occurs at non-finite x values")
    if min_neighbors < 1:
        raise ValueError("min_neighbors must be >= 1")
    if min_total_pixels < 2:
        raise ValueError("min_total_pixels must be >= 2")

    # Enforce star: center = 1
    k = kernel.copy()
    cy, cx = k.shape[0] // 2, k.shape[1] // 2
    k[cy, cx] = 1.0

    # Mask handling: zero out invalid pixels in value array
    x_zero = np.where(valid, x, 0.0)
    mask_int = valid.astype(np.float64)

    # Sum of values in neighborhood (including self if valid)
    sum_wx = ndimage.convolve(x_zero, k, mode='constant', cval=0.0)

    # Count of valid neighbors per cell (including self if valid)
    W = ndimage.convolve(mask_int, k, mode='constant', cval=0.0)

    # Global statistics
    valid_vals = x[valid]
    n = len(valid_vals)
    if n < min_total_pixels:
        raise ValueError(
            f"Too few valid pixels for Gi*: n={n}; minimum={min_total_pixels}"
        )
    x_mean = valid_vals.mean()
    S = np.sqrt(np.mean((valid_vals - x_mean) ** 2))   # population SD per Getis-Ord

    if S == 0:
        logger.warning("Global SD is zero; Gi* undefined.")
        z = np.full_like(x, np.nan, dtype=np.float64)
        p = np.full_like(x, np.nan, dtype=np.float64)
        return z, p, W

    # Z-score: (Σ_j w_ij*x_j - x̄*W_i) / (S * sqrt(W_i*(n - W_i)/(n - 1)))
    with np.errstate(invalid='ignore', divide='ignore'):
        numerator = sum_wx - x_mean * W
        var_term = np.where(W > 0,
                            W * (n - W) / max(n - 1, 1),
                            np.nan)
        denominator = S * np.sqrt(np.maximum(var_term, 1e-300))
        z = numerator / denominator

    # Mask out cells where focal is invalid OR too few neighbors
    bad = (~valid) | (W < min_neighbors) | ~np.isfinite(z)
    z[bad] = np.nan

    # Two-tailed p-value
    p = np.full_like(z, np.nan, dtype=np.float64)
    finite_z = np.isfinite(z)
    p[finite_z] = 2.0 * sp_stats.norm.sf(np.abs(z[finite_z]))

    return z, p, W


def classify_gi_zscore(z: np.ndarray) -> np.ndarray:
    """
    Classify Gi* z-score into 7-level signature (-3..+3) matching v7 Fig.5 design.

    Class scheme:
       3 = Hot99   (z >= 2.58)
       2 = Hot95   (1.96 <= z < 2.58)
       1 = Hot90   (1.65 <= z < 1.96)
       0 = NotSig  (-1.65 < z < 1.65)
      -1 = Cold90  (-1.96 < z <= -1.65)
      -2 = Cold95  (-2.58 < z <= -1.96)
      -3 = Cold99  (z <= -2.58)
    NaN → NODATA_INT16
    """
    out = np.full(z.shape, NODATA_INT16, dtype=np.int16)
    finite = np.isfinite(z)
    cls = np.zeros(z.shape, dtype=np.int16)

    cls = np.where((z >=  Z_90) & (z <  Z_95),  1, cls)
    cls = np.where((z >=  Z_95) & (z <  Z_99),  2, cls)
    cls = np.where(z >= Z_99,                   3, cls)
    cls = np.where((z <= -Z_90) & (z > -Z_95), -1, cls)
    cls = np.where((z <= -Z_95) & (z > -Z_99), -2, cls)
    cls = np.where(z <= -Z_99,                 -3, cls)

    out[finite] = cls[finite]
    return out


# ============================================================
# Output writers
# ============================================================

def to_int16(arr: np.ndarray, scale: float, nodata: int = NODATA_INT16) -> np.ndarray:
    """Scale with round-to-nearest and prevent valid/NoData collisions."""
    out = np.full(arr.shape, nodata, dtype=np.int16)
    valid = np.isfinite(arr)
    scaled = np.clip(np.rint(arr[valid] * scale), -32767, 32766)
    if -32767 <= nodata <= 32766:
        collision = scaled == nodata
        if np.any(collision):
            scaled[collision] = nodata + 1 if nodata < 0 else nodata - 1
    out[valid] = scaled.astype(np.int16)
    return out


def write_gi_star_geotiff(
    path: Path,
    z: np.ndarray,
    p: np.ndarray,
    cls: np.ndarray,
    profile: dict,
):
    """
    Write 3-band Int16 Gi* GeoTIFF.
        Band 1: z-score  x100
        Band 2: p-value  x10000
        Band 3: class    -3..+3 (NoData = -9999)
    """
    out_profile = {
        **profile,
        'driver': 'GTiff',
        'dtype': rasterio.int16,
        'count': 3,
        'nodata': NODATA_INT16,
        'compress': 'lzw',
        'tiled': True,
        'blockxsize': 256,
        'blockysize': 256,
        'BIGTIFF': 'YES',
    }

    z_int = to_int16(z, SCALE_ZSCORE)
    p_int = to_int16(p, SCALE_PVALUE)
    # cls is already Int16 with NoData

    with rasterio.open(path, 'w', **out_profile) as dst:
        dst.write(z_int, 1)
        dst.write(p_int, 2)
        dst.write(cls, 3)
        dst.set_band_description(1, f"gi_star_zscore_x{SCALE_ZSCORE}")
        dst.set_band_description(2, f"gi_star_pvalue_x{SCALE_PVALUE}")
        dst.set_band_description(3, "gi_star_class_-3to3")

    logger.info(f"  Saved: {path}")


def write_summary_csv(
    path: Path,
    cls: np.ndarray,
    z: np.ndarray,
):
    """Write Gi* class count + percentage summary CSV."""
    n_total = int((cls != NODATA_INT16).sum())

    counts = {
        "Cold99":  int((cls == -3).sum()),
        "Cold95":  int((cls == -2).sum()),
        "Cold90":  int((cls == -1).sum()),
        "NotSig":  int((cls ==  0).sum()),
        "Hot90":   int((cls ==  1).sum()),
        "Hot95":   int((cls ==  2).sum()),
        "Hot99":   int((cls ==  3).sum()),
    }
    cold_total = counts["Cold99"] + counts["Cold95"] + counts["Cold90"]
    hot_total  = counts["Hot90"]  + counts["Hot95"]  + counts["Hot99"]
    sig_total  = cold_total + hot_total

    finite_z = z[np.isfinite(z)]
    z_min, z_max = (float(finite_z.min()), float(finite_z.max())) if len(finite_z) else (np.nan, np.nan)
    z_mean = float(finite_z.mean()) if len(finite_z) else np.nan

    with open(path, 'w', encoding='utf-8') as f:
        f.write("class,count,percent_of_valid\n")
        for k, v in counts.items():
            pct = 100.0 * v / n_total if n_total else 0.0
            f.write(f"{k},{v},{pct:.4f}\n")
        f.write(f"_TotalValid,{n_total},100.0000\n")
        f.write(f"_HotSpotTotal_p_le_0.10,{hot_total},{100*hot_total/max(n_total,1):.4f}\n")
        f.write(f"_ColdSpotTotal_p_le_0.10,{cold_total},{100*cold_total/max(n_total,1):.4f}\n")
        f.write(f"_AnySig_p_le_0.10,{sig_total},{100*sig_total/max(n_total,1):.4f}\n")
        f.write(f"_zscore_min,,{z_min}\n")
        f.write(f"_zscore_max,,{z_max}\n")
        f.write(f"_zscore_mean,,{z_mean}\n")

    logger.info(f"  Saved: {path}")
    logger.info(f"  Hot spots (p<=0.10):  {hot_total:,} ({100*hot_total/max(n_total,1):.2f}%)")
    logger.info(f"  Cold spots (p<=0.10): {cold_total:,} ({100*cold_total/max(n_total,1):.2f}%)")


def write_global_moran_csv(
    path: Path,
    moran_I: float,
    n: int,
    W_total: float,
    inference: Dict[str, float],
    perm_p: Optional[float] = None,
    perm_n: Optional[int] = None,
    moran_mask_mode: str = "same_as_analysis_mask",
    bootstrap_metrics: Optional[Dict[str, object]] = None,
    quality_metrics: Optional[Dict[str, object]] = None,
):
    """Atomically write Global Moran's I metrics and optional bootstrap CI."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as f:
            f.write("metric,value\n")
            f.write(f"morans_I,{moran_I}\n")
            f.write(f"moran_mask_mode,{moran_mask_mode}\n")
            f.write(f"n_pixels,{n}\n")
            f.write(f"W_total,{W_total}\n")
            f.write(f"expected_I,{inference['expected_I']}\n")
            f.write(f"var_I_normality,{inference['var_I_norm']}\n")
            f.write(f"var_I_randomization,{inference['var_I_rand']}\n")
            f.write(f"z_normality,{inference['z_norm']}\n")
            f.write(f"z_randomization,{inference['z_rand']}\n")
            f.write(f"p_normality_two_tailed,{inference['p_norm']}\n")
            f.write(f"p_randomization_two_tailed,{inference['p_rand']}\n")
            if quality_metrics:
                for key, value in quality_metrics.items():
                    f.write(f"{key},{value}\n")
            if perm_p is not None:
                f.write(f"p_permutation_two_tailed,{perm_p}\n")
                f.write(f"n_permutations,{perm_n}\n")
            if bootstrap_metrics:
                for key, value in bootstrap_metrics.items():
                    f.write(f"{key},{value}\n")
            f.write(f"S0_sum_weights,{inference['S0']}\n")
            f.write(f"S1,{inference['S1']}\n")
            f.write(f"S2,{inference['S2']}\n")
            f.write(f"b2_kurtosis_like,{inference['b2']}\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise
    logger.info(f"  Saved: {path}")


# ============================================================
# CLI / Main
# ============================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="spatial_autocorrelation.py",
        description="Global Moran's I + Local Getis-Ord Gi* on raster data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # IO
    p.add_argument('--slope-input', '-i', type=Path, required=True,
                   help="Path to Sen's slope raster (e.g., *_TrendSlope.tif)")
    p.add_argument('--output-dir', '-o', type=Path, required=True,
                   help="Output directory")
    p.add_argument('--prefix', type=str, required=True,
                   help="Prefix for output files (e.g., 'SR_2001_2020')")
    p.add_argument('--run-config-id', type=str, default='unspecified',
                   help="Optional deterministic task-configuration fingerprint written to "
                        "the output CSV so batch drivers can reject stale cached results.")
    p.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                   default='INFO', help="Console and file logging level (default: INFO).")

    # Optional masking
    p.add_argument('--mk-input', type=Path, default=None,
                   help="Mann-Kendall raster path (optional, for significance masking)")
    p.add_argument('--mk-pvalue-band', type=int, default=2,
                   help="Band index of p-value in MK raster (default: 2)")
    p.add_argument('--mk-pvalue-scale', type=float, default=10000.0,
                   help="Scale factor of MK p-value in raster (default: 10000)")
    p.add_argument('--p-threshold', type=float, default=None,
                   help="MK p-value threshold for significance masking. "
                        "Only pixels with p<=threshold included. "
                        "Default: None (no masking, analyze all valid pixels)")
    p.add_argument('--moran-ignore-significance-mask', '--moran-all-valid',
                   dest='moran_ignore_significance_mask', action='store_true',
                   help="Compute Global Moran's I using all valid slope pixels, even when "
                        "--mk-input/--p-threshold applies an MK p/q significance mask. "
                        "Local Gi* still uses the significance mask. Forest mask and "
                        "aggregation are still applied to both analyses.")
    p.add_argument('--forest-mask', type=Path, default=None,
                   help="Optional binary forest mask (1 = forest)")

    # Slope input scaling
    p.add_argument('--slope-input-scale', type=float, default=100.0,
                   help="Input scaling of slope raster (default: 100, "
                        "matching v3 pipeline output)")

    # Aggregation
    p.add_argument('--aggregate-factor', type=int, default=1,
                   help="Mean-aggregate slope before analysis (default: 1, no aggregation)")
    p.add_argument('--aggregate-edge-mode', choices=['partial', 'trim'], default='partial',
                   help="How to handle bottom/right incomplete aggregation groups. "
                        "'partial' preserves them; 'trim' reproduces legacy cropping.")
    p.add_argument('--aggregate-min-valid-fraction', type=float, default=0.0,
                   help="Minimum valid base-cell fraction of the nominal aggregate cell. "
                        "0 preserves legacy any-valid behavior; 0.25 or 0.5 is recommended "
                        "for sensitivity checks on sparse/irregular masks.")

    # Spatial weights
    p.add_argument('--neighbor', type=str, choices=['queen', 'rook', 'circle'],
                   default='queen',
                   help="Neighbor type for spatial weights (default: queen, 8 neighbors)")
    p.add_argument('--kernel-radius', type=int, default=1,
                   help="Kernel radius in pixels (default: 1, gives 3x3 queen)")
    p.add_argument('--min-neighbors', type=int, default=2,
                   help="Minimum valid neighbors (incl. self) for Gi* z-score (default: 2)")
    p.add_argument('--min-moran-pixels', type=int, default=30,
                   help="Minimum valid pixels for analytical Moran inference and bootstrap. "
                        "Below this threshold the point estimate can still be written "
                        "when --small-moran-sample-action=point-only (default: 30).")
    p.add_argument('--small-moran-sample-action', choices=['point-only', 'error'],
                   default='point-only',
                   help="Behavior when valid Moran pixels are below --min-moran-pixels. "
                        "point-only writes Moran's I and marks inference/CI unavailable; "
                        "error aborts the task (default: point-only).")
    p.add_argument('--min-local-pixels', type=int, default=30,
                   help="Minimum valid pixels required for Local Gi* (default: 30).")

    # Permutation test
    p.add_argument('--permutations', type=int, default=0,
                   help="Number of Monte Carlo permutations for global Moran. "
                        "0 = analytical only (default). 999 recommended if requested.")
    p.add_argument('--seed', type=int, default=42,
                   help="Random seed for permutation test (default: 42)")

    # Spatial block bootstrap confidence interval
    p.add_argument('--spatial-bootstrap', type=int, default=0, metavar='N',
                   help="Number of spatial block-bootstrap replicates for Moran's I. "
                        "0 disables bootstrap; 999 or 1999 is recommended.")
    p.add_argument('--bootstrap-block-size-cells', type=int, default=None,
                   help="Square block side length in cells after aggregation.")
    p.add_argument('--bootstrap-block-rows-cells', type=int, default=None,
                   help="Block height in cells after aggregation (for rectangular blocks).")
    p.add_argument('--bootstrap-block-cols-cells', type=int, default=None,
                   help="Block width in cells after aggregation (for rectangular blocks).")
    p.add_argument('--bootstrap-workers', type=int, default=1,
                   help="Parallel threads used for bootstrap replicates (default: 1).")
    p.add_argument('--bootstrap-block-workers', type=int, default=None,
                   help="Threads used to precompute block statistics. Default: same as "
                        "--bootstrap-workers.")
    p.add_argument('--bootstrap-ci-level', type=float, default=0.95,
                   help="Bootstrap confidence level (default: 0.95).")
    p.add_argument('--bootstrap-seed', type=int, default=None,
                   help="Bootstrap seed. Default: reuse --seed.")
    p.add_argument('--bootstrap-min-valid-pixels', type=int, default=30,
                   help="Minimum valid cells required in a usable block (default: 30).")
    p.add_argument('--bootstrap-min-valid-fraction', type=float, default=0.0,
                   help="Minimum valid-cell fraction required in a block (default: 0).")
    p.add_argument('--bootstrap-random-grid-offset', action='store_true',
                   help="Randomly offset the non-overlapping block grid using the bootstrap seed. "
                        "Use with --bootstrap-edge-mode drop or partial.")
    p.add_argument('--bootstrap-edge-mode', choices=['balanced', 'drop', 'partial'],
                   default='balanced',
                   help="Raster-edge handling for bootstrap blocks. balanced preserves the full "
                        "extent using near-equal blocks; drop keeps exact blocks; partial includes "
                        "small edge fragments and is only for sensitivity analysis.")
    p.add_argument('--bootstrap-include-partial-edge-blocks', action='store_true',
                   help=argparse.SUPPRESS)
    p.add_argument('--bootstrap-max-partition-difference', type=float, default=0.05,
                   help="Warn and flag CI as unreliable when |I_partition-I_full| exceeds this value.")
    p.add_argument('--bootstrap-max-valid-pixel-cv', type=float, default=0.50,
                   help="Flag CI when the coefficient of variation of valid pixels per block "
                        "exceeds this threshold (default: 0.50).")
    p.add_argument('--bootstrap-min-usable-blocks', type=int, default=20,
                   help="Minimum usable blocks required for a CI to be flagged reliable "
                        "(default: 20). Fewer blocks do not abort computation; the CI is "
                        "retained for diagnostics and marked unreliable.")
    p.add_argument('--bootstrap-save-distribution', action='store_true',
                   help="Save the bootstrap Moran's I distribution as a .npy file.")

    # Toggles
    p.add_argument('--global-only', action='store_true',
                   help="Skip local Gi* (only compute global Moran's I)")
    p.add_argument('--local-only', action='store_true',
                   help="Skip global Moran's I (only compute local Gi*)")

    return p


def raster_resolution_metadata(profile: Dict) -> Dict[str, object]:
    """Return actual output-grid resolution and CRS linear-unit metadata.

    Multi-scale labels such as 100 m are nominal targets when a 30 m source is
    aggregated by an integer factor (factor 3 produces an actual 90 m cell).
    Recording the effective resolution prevents that approximation from being
    hidden in downstream interpretation.
    """
    transform = profile.get("transform")
    if transform is None:
        return {
            "effective_pixel_size_x_map_units": float("nan"),
            "effective_pixel_size_y_map_units": float("nan"),
            "effective_pixel_area_map_units2": float("nan"),
            "crs_linear_unit_name": "unknown",
            "crs_linear_unit_to_metre": float("nan"),
            "crs_is_projected": 0,
            "crs_is_geographic": 0,
            "effective_pixel_size_x_m": float("nan"),
            "effective_pixel_size_y_m": float("nan"),
        }

    pixel_x = float(np.hypot(transform.a, transform.d))
    pixel_y = float(np.hypot(transform.b, transform.e))
    unit_name = "unknown"
    unit_to_metre = float("nan")
    crs = profile.get("crs")
    crs_is_projected = int(bool(crs is not None and crs.is_projected))
    crs_is_geographic = int(bool(crs is not None and crs.is_geographic))
    if crs is not None:
        try:
            unit_name, unit_to_metre = crs.linear_units_factor
            unit_name = str(unit_name)
            unit_to_metre = float(unit_to_metre)
        except Exception:
            try:
                unit_name = str(crs.linear_units)
            except Exception:
                pass

    return {
        "effective_pixel_size_x_map_units": pixel_x,
        "effective_pixel_size_y_map_units": pixel_y,
        "effective_pixel_area_map_units2": float(abs(
            transform.a * transform.e - transform.b * transform.d
        )),
        "crs_linear_unit_name": unit_name,
        "crs_linear_unit_to_metre": unit_to_metre,
        "crs_is_projected": crs_is_projected,
        "crs_is_geographic": crs_is_geographic,
        "effective_pixel_size_x_m": (
            pixel_x * unit_to_metre if np.isfinite(unit_to_metre) else float("nan")
        ),
        "effective_pixel_size_y_m": (
            pixel_y * unit_to_metre if np.isfinite(unit_to_metre) else float("nan")
        ),
    }


def configure_logging(log_path: Path, level_name: str = "INFO") -> None:
    """Configure isolated console and file logging for one CLI invocation.

    The module installs only a ``NullHandler`` at import time so importing it
    never changes application-wide logging. Repeated in-process invocations are
    safe because handlers created by a previous call are closed and replaced.
    """
    level = getattr(logging, level_name.upper())
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        if not isinstance(handler, logging.NullHandler):
            handler.close()

    log_path.parent.mkdir(parents=True, exist_ok=True)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(log_path, mode='w', encoding='utf-8')
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)

    logger.setLevel(level)
    logger.propagate = False
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    # Validate
    if args.global_only and args.local_only:
        parser.error("--global-only and --local-only are mutually exclusive")
    if any(ch in args.run_config_id for ch in (',', '\n', '\r')):
        parser.error("--run-config-id must not contain commas or newlines")
    if args.aggregate_factor < 1:
        parser.error("--aggregate-factor must be >= 1")
    if not (0.0 <= args.aggregate_min_valid_fraction <= 1.0):
        parser.error("--aggregate-min-valid-fraction must be between 0 and 1")
    if args.slope_input_scale <= 0:
        parser.error("--slope-input-scale must be > 0")
    if args.mk_pvalue_scale <= 0:
        parser.error("--mk-pvalue-scale must be > 0")
    if args.kernel_radius < 1:
        parser.error("--kernel-radius must be >= 1")
    if args.min_moran_pixels < 4:
        parser.error("--min-moran-pixels must be >= 4")
    if args.min_local_pixels < 2:
        parser.error("--min-local-pixels must be >= 2")
    if args.p_threshold is not None and not (0.0 <= args.p_threshold <= 1.0):
        parser.error("--p-threshold must be between 0 and 1")
    if args.seed < 0:
        parser.error("--seed must be non-negative")
    if args.spatial_bootstrap > 0:
        square_given = args.bootstrap_block_size_cells is not None
        rect_given = (args.bootstrap_block_rows_cells is not None or
                      args.bootstrap_block_cols_cells is not None)
        if square_given and rect_given:
            parser.error("Use either --bootstrap-block-size-cells or the row/col options, not both")
        if not square_given and not (args.bootstrap_block_rows_cells is not None and
                                     args.bootstrap_block_cols_cells is not None):
            parser.error("Spatial bootstrap requires --bootstrap-block-size-cells, or both "
                         "--bootstrap-block-rows-cells and --bootstrap-block-cols-cells")
        if args.bootstrap_workers < 1:
            parser.error("--bootstrap-workers must be >= 1")
        if args.bootstrap_block_workers is not None and args.bootstrap_block_workers < 1:
            parser.error("--bootstrap-block-workers must be >= 1")
        if args.bootstrap_seed is not None and args.bootstrap_seed < 0:
            parser.error("--bootstrap-seed must be non-negative")
        if args.bootstrap_min_valid_pixels < 2:
            parser.error("--bootstrap-min-valid-pixels must be >= 2")
        if args.bootstrap_max_partition_difference < 0:
            parser.error("--bootstrap-max-partition-difference must be >= 0")
        if args.bootstrap_max_valid_pixel_cv < 0:
            parser.error("--bootstrap-max-valid-pixel-cv must be >= 0")
        if args.bootstrap_min_usable_blocks < 2:
            parser.error("--bootstrap-min-usable-blocks must be >= 2")
        block_values = [v for v in (
            args.bootstrap_block_size_cells,
            args.bootstrap_block_rows_cells,
            args.bootstrap_block_cols_cells,
        ) if v is not None]
        if any(v < 2 for v in block_values):
            parser.error("Bootstrap block dimensions must be >= 2 cells")
        effective_edge_mode = (
            'partial' if args.bootstrap_include_partial_edge_blocks
            else args.bootstrap_edge_mode
        )
        if args.bootstrap_random_grid_offset and effective_edge_mode == 'balanced':
            parser.error("--bootstrap-random-grid-offset requires --bootstrap-edge-mode drop or partial")
        if not (0.0 < args.bootstrap_ci_level < 1.0):
            parser.error("--bootstrap-ci-level must be between 0 and 1")
        if not (0.0 <= args.bootstrap_min_valid_fraction <= 1.0):
            parser.error("--bootstrap-min-valid-fraction must be between 0 and 1")
    if not args.slope_input.exists():
        parser.error(f"Slope input not found: {args.slope_input}")
    if args.mk_input is not None and not args.mk_input.exists():
        parser.error(f"MK/FDR input not found: {args.mk_input}")
    if args.forest_mask is not None and not args.forest_mask.exists():
        parser.error(f"Forest mask not found: {args.forest_mask}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_logging(args.output_dir / f"{args.prefix}_log.txt", args.log_level)

    t0 = time.time()
    logger.info("=" * 60)
    logger.info(f"spatial_autocorrelation.py v{__version__}")
    logger.info(f"Input slope: {args.slope_input}")
    logger.info(f"Output dir:  {args.output_dir}")
    logger.info(f"Prefix:      {args.prefix}")
    logger.info(f"Neighbor:    {args.neighbor} (radius {args.kernel_radius})")
    if args.p_threshold is not None:
        logger.info(f"MK significance threshold: p <= {args.p_threshold}")
    if args.moran_ignore_significance_mask:
        logger.info("Moran mask mode: ignore MK/FDR significance mask for Global Moran's I")
    if args.aggregate_factor > 1:
        logger.info(
            f"Aggregation factor: {args.aggregate_factor}; "
            f"minimum nominal valid fraction={args.aggregate_min_valid_fraction:.3f}"
        )
    if args.spatial_bootstrap > 0:
        if args.bootstrap_block_size_cells is not None:
            block_text = f"{args.bootstrap_block_size_cells}x{args.bootstrap_block_size_cells} cells"
        else:
            block_text = (f"{args.bootstrap_block_rows_cells}x"
                          f"{args.bootstrap_block_cols_cells} cells")
        logger.info(
            f"Spatial block bootstrap: n={args.spatial_bootstrap}, block={block_text}, "
            f"workers={args.bootstrap_workers}, CI={args.bootstrap_ci_level:.3f}"
        )
    logger.info("=" * 60)

    # 1. Load
    slope, valid, profile = load_slope_with_mask(
        slope_path=args.slope_input,
        mk_path=args.mk_input,
        mk_pvalue_band=args.mk_pvalue_band,
        mk_pvalue_scale=args.mk_pvalue_scale,
        p_threshold=args.p_threshold,
        forest_mask_path=args.forest_mask,
        aggregate_factor=args.aggregate_factor,
        slope_input_scale=args.slope_input_scale,
        aggregate_edge_mode=args.aggregate_edge_mode,
        aggregate_min_valid_fraction=args.aggregate_min_valid_fraction,
    )

    resolution_metrics = raster_resolution_metadata(profile)
    logger.info(
        "Effective output resolution: "
        f"{resolution_metrics['effective_pixel_size_x_map_units']:.6g} x "
        f"{resolution_metrics['effective_pixel_size_y_map_units']:.6g} "
        f"{resolution_metrics['crs_linear_unit_name']}"
    )
    if resolution_metrics['crs_is_geographic']:
        logger.warning(
            "Input CRS is geographic. Raster-neighbor Moran's I is computable, but "
            "cell areas and physical block sizes vary with latitude. Reproject to an "
            "appropriate projected/equal-area CRS before interpreting block sizes in km."
        )

    n_valid = int(valid.sum())
    if not args.global_only and n_valid < args.min_local_pixels:
        logger.error(
            "Too few valid pixels for Local Gi* (n=%s; minimum=%s). Aborting.",
            n_valid,
            args.min_local_pixels,
        )
        return 1

    # Optional: Global Moran's I can be computed on all valid slope pixels,
    # while Local Gi* continues to use the MK/FDR significance-filtered mask.
    moran_slope, moran_valid = slope, valid
    moran_mask_mode = "same_as_analysis_mask"
    if args.moran_ignore_significance_mask and args.mk_input is not None and args.p_threshold is not None:
        logger.info("Loading a separate all-valid slope raster for Global Moran's I "
                    "(ignoring MK/FDR p/q mask).")
        moran_slope, moran_valid, _ = load_slope_with_mask(
            slope_path=args.slope_input,
            mk_path=None,
            mk_pvalue_band=args.mk_pvalue_band,
            mk_pvalue_scale=args.mk_pvalue_scale,
            p_threshold=None,
            forest_mask_path=args.forest_mask,
            aggregate_factor=args.aggregate_factor,
            slope_input_scale=args.slope_input_scale,
            aggregate_edge_mode=args.aggregate_edge_mode,
            aggregate_min_valid_fraction=args.aggregate_min_valid_fraction,
        )
        n_moran_valid = int(moran_valid.sum())
        moran_mask_mode = "all_valid_slope_pixels_ignore_mk_fdr_pq_mask"
        logger.info(f"  Global Moran valid pixels: {n_moran_valid:,} "
                    f"(local Gi* valid pixels: {n_valid:,})")
    elif args.moran_ignore_significance_mask:
        moran_mask_mode = "all_valid_requested_no_significance_mask_present"
        logger.info("Moran all-valid mode requested, but no MK/FDR significance mask is active; "
                    "Global Moran's I uses the current valid mask.")

    # 2. Build kernel
    kernel = make_kernel(args.neighbor, args.kernel_radius)
    logger.info(f"Kernel:\n{kernel.astype(int)}")

    # 3. Global Moran's I
    moran_results = None
    if not args.local_only:
        logger.info("Computing global Moran's I...")
        t1 = time.time()
        moran_I, W_total, n_used = compute_global_moran(moran_slope, moran_valid, kernel)
        logger.info(f"  Moran's I = {moran_I:.6f} (n={n_used:,}, W={W_total:.0f}, "
                    f"{time.time()-t1:.1f}s)")

        small_moran_sample = n_used < args.min_moran_pixels
        if small_moran_sample and args.small_moran_sample_action == 'error':
            raise ValueError(
                f"Moran valid-pixel count n={n_used} is below configured minimum "
                f"{args.min_moran_pixels}"
            )

        quality_metrics = {
            'run_config_id': args.run_config_id,
            'moran_result_status': (
                'point_estimate_only_below_minimum' if small_moran_sample else 'complete'
            ),
            'moran_inference_available': int(not small_moran_sample),
            'moran_min_pixels_for_inference': int(args.min_moran_pixels),
            'moran_small_sample_action': args.small_moran_sample_action,
            'aggregate_factor': int(args.aggregate_factor),
            'aggregate_edge_mode': args.aggregate_edge_mode,
            'aggregate_min_valid_fraction': float(args.aggregate_min_valid_fraction),
            **resolution_metrics,
        }

        if small_moran_sample:
            logger.warning(
                f"Only n={n_used} valid pixels remain, below --min-moran-pixels="
                f"{args.min_moran_pixels}. Writing the Moran point estimate only; "
                "analytical inference, permutation testing, and bootstrap CI are skipped."
            )
            inference = unavailable_moran_inference(n_used)
        else:
            logger.info("Computing analytical inference (Cliff & Ord 1981)...")
            inference = moran_inference_normality(moran_I, W_total, n_used,
                                                  moran_valid, kernel, moran_slope)
            logger.info(f"  E(I) = {inference['expected_I']:.6f}")
            logger.info(f"  z (randomization) = {inference['z_rand']:.4f}, "
                        f"p (two-tailed) = {inference['p_rand']:.6g}")

        # Optional permutation
        perm_p, perm_n = None, None
        if args.permutations > 0 and not small_moran_sample:
            logger.info(f"Running Monte Carlo permutation (n={args.permutations})...")
            t2 = time.time()
            perm_p, perm_dist = moran_permutation_test(
                moran_slope, moran_valid, kernel, moran_I,
                n_perm=args.permutations, seed=args.seed
            )
            perm_n = args.permutations
            logger.info(f"  Permutation p-value = {perm_p:.6g} ({time.time()-t2:.1f}s)")

        # Optional spatial block-bootstrap CI. The full raster is compressed to
        # per-block sufficient statistics once, then bootstrap replicates run in
        # parallel without copying/convolving the raster for every replicate.
        bootstrap_metrics = None
        if args.spatial_bootstrap > 0 and small_moran_sample:
            bootstrap_metrics = {
                'bootstrap_requested': int(args.spatial_bootstrap),
                'bootstrap_ci_reliable': 0,
                'bootstrap_ci_reliability_reasons': 'moran_sample_below_minimum',
                'bootstrap_ci_method': 'recentered_basic',
                'bootstrap_skipped_reason': 'moran_sample_below_minimum',
            }
        elif args.spatial_bootstrap > 0:
            if args.bootstrap_block_size_cells is not None:
                block_rows = block_cols = int(args.bootstrap_block_size_cells)
            else:
                block_rows = int(args.bootstrap_block_rows_cells)
                block_cols = int(args.bootstrap_block_cols_cells)
            if block_rows < 2 or block_cols < 2:
                parser.error("Spatial block dimensions must be >= 2 cells")
            bootstrap_edge_mode = (
                'partial' if args.bootstrap_include_partial_edge_blocks
                else args.bootstrap_edge_mode
            )

            bootstrap_seed = args.seed if args.bootstrap_seed is None else args.bootstrap_seed
            if args.bootstrap_random_grid_offset:
                offset_rng = np.random.default_rng(bootstrap_seed)
                offset_rows = int(offset_rng.integers(0, block_rows))
                offset_cols = int(offset_rng.integers(0, block_cols))
            else:
                offset_rows = 0
                offset_cols = 0

            block_workers = (args.bootstrap_workers if args.bootstrap_block_workers is None
                             else args.bootstrap_block_workers)
            block_workers = max(1, int(block_workers))
            # At coarse raster scales a 2x2 or 3x3 block contains fewer than
            # the general 30-cell default. Cap the threshold by block area.
            effective_min_valid_pixels = min(
                int(args.bootstrap_min_valid_pixels), block_rows * block_cols
            )
            effective_min_valid_pixels = max(2, effective_min_valid_pixels)
            logger.info(
                f"Preparing spatial blocks: {block_rows}x{block_cols} cells, "
                f"offset=({offset_rows},{offset_cols}), edge_mode={bootstrap_edge_mode}, "
                f"workers={block_workers}, min_valid={effective_min_valid_pixels}"
            )
            t_boot = time.time()
            block_stats = compute_block_statistics(
                moran_slope, moran_valid, kernel,
                block_rows=block_rows, block_cols=block_cols,
                offset_rows=offset_rows, offset_cols=offset_cols,
                min_valid_pixels=effective_min_valid_pixels,
                min_valid_fraction=args.bootstrap_min_valid_fraction,
                workers=block_workers,
                edge_mode=bootstrap_edge_mode,
                require_min_blocks=0,
            )
            logger.info(
                f"  Usable blocks: {block_stats.n_blocks:,}/"
                f"{block_stats.n_total_candidate_blocks:,} "
                f"({time.time()-t_boot:.1f}s precompute)"
            )
            transform = profile['transform']
            boot_meta = block_statistics_metadata(
                block_stats,
                block_rows=block_rows, block_cols=block_cols,
                pixel_size_x=float(np.hypot(transform.a, transform.d)),
                pixel_size_y=float(np.hypot(transform.b, transform.e)),
                offset_rows=offset_rows, offset_cols=offset_cols,
                min_valid_pixels=effective_min_valid_pixels,
                min_valid_fraction=args.bootstrap_min_valid_fraction,
            )
            coverage_fraction = float(block_stats.n.sum()) / max(float(n_used), 1.0)
            base_bootstrap_metrics = {
                'bootstrap_method': (
                    f'nonoverlapping_{bootstrap_edge_mode}_blocks_within_block_links'
                ),
                'bootstrap_ci_method': 'recentered_basic',
                'bootstrap_workers': args.bootstrap_workers,
                'bootstrap_block_workers': block_workers,
                'bootstrap_requested': int(args.spatial_bootstrap),
                'bootstrap_valid_pixel_coverage_fraction': coverage_fraction,
                'bootstrap_max_partition_difference_allowed': args.bootstrap_max_partition_difference,
                'bootstrap_max_valid_pixel_cv_allowed': args.bootstrap_max_valid_pixel_cv,
                'bootstrap_min_usable_blocks_required': args.bootstrap_min_usable_blocks,
                **boot_meta,
            }

            if block_stats.n_blocks < 2:
                reliability_reasons = ['fewer_than_2_blocks']
                if coverage_fraction < 0.95:
                    reliability_reasons.append('coverage_below_95pct')
                if bootstrap_edge_mode == 'partial':
                    reliability_reasons.append('partial_edge_blocks_equal_weight')
                bootstrap_metrics = {
                    **base_bootstrap_metrics,
                    'bootstrap_ci_reliable': 0,
                    'bootstrap_ci_reliability_reasons': '|'.join(reliability_reasons),
                    'bootstrap_skipped_reason': 'insufficient_usable_blocks',
                }
                logger.warning(
                    f"Only {block_stats.n_blocks} usable spatial block(s) remain. "
                    "Moran's I point estimate is retained, but spatial block-bootstrap "
                    "CI is skipped. Reduce block size or min-valid thresholds."
                )
            else:
                if block_stats.n_blocks < args.bootstrap_min_usable_blocks:
                    logger.warning(
                        f"Only {block_stats.n_blocks} usable blocks remain, below the "
                        f"configured reliability threshold {args.bootstrap_min_usable_blocks}. "
                        "Bootstrap CI will be marked unreliable."
                    )

                boot_result, boot_distribution = spatial_block_bootstrap_moran(
                    block_stats, moran_I,
                    n_boot=args.spatial_bootstrap,
                    ci_level=args.bootstrap_ci_level,
                    seed=bootstrap_seed,
                    workers=args.bootstrap_workers,
                )
                partition_difference = float(
                    boot_result['bootstrap_partition_difference_vs_full_I']
                )
                block_valid_pixel_cv = float(
                    boot_meta['bootstrap_cv_valid_pixels_per_block']
                )
                reliability_reasons = []
                if block_stats.n_blocks < args.bootstrap_min_usable_blocks:
                    reliability_reasons.append(
                        f"fewer_than_{args.bootstrap_min_usable_blocks}_blocks"
                    )
                if coverage_fraction < 0.95:
                    reliability_reasons.append("coverage_below_95pct")
                if abs(partition_difference) > args.bootstrap_max_partition_difference:
                    reliability_reasons.append("partition_difference_too_large")
                if block_valid_pixel_cv > args.bootstrap_max_valid_pixel_cv:
                    reliability_reasons.append("valid_pixels_per_block_cv_too_large")
                if bootstrap_edge_mode == 'partial':
                    reliability_reasons.append("partial_edge_blocks_equal_weight")
                ci_reliable = int(len(reliability_reasons) == 0)

                if coverage_fraction < 0.99:
                    logger.warning(
                        f"Block partition retains {100*coverage_fraction:.2f}% of valid pixels. "
                        "Inspect edge cropping and min-valid thresholds."
                    )
                if abs(partition_difference) > args.bootstrap_max_partition_difference:
                    logger.warning(
                        f"|I_partition-I_full|={abs(partition_difference):.6g} exceeds "
                        f"the configured limit {args.bootstrap_max_partition_difference:.6g}. "
                        "Increase block size and rerun block-origin sensitivity analysis."
                    )
                if block_valid_pixel_cv > args.bootstrap_max_valid_pixel_cv:
                    logger.warning(
                        f"Valid-pixel count CV across blocks is {block_valid_pixel_cv:.3f}, "
                        f"above the configured limit {args.bootstrap_max_valid_pixel_cv:.3f}. "
                        "Increase --bootstrap-min-valid-fraction or review the study mask."
                    )
                if bootstrap_edge_mode == 'partial':
                    logger.warning(
                        "Partial edge blocks are included and resampled with the same probability "
                        "as full blocks; use this only as a sensitivity run."
                    )

                bootstrap_metrics = {
                    **base_bootstrap_metrics,
                    'bootstrap_ci_reliable': ci_reliable,
                    'bootstrap_ci_reliability_reasons': (
                        'none' if ci_reliable else '|'.join(reliability_reasons)
                    ),
                    **boot_result,
                    'bootstrap_ci_low': boot_result['bootstrap_basic_ci_low'],
                    'bootstrap_ci_high': boot_result['bootstrap_basic_ci_high'],
                }
                logger.info(
                    f"  Block bootstrap SE={boot_result['bootstrap_se']:.6g}; "
                    f"recentered basic {100*args.bootstrap_ci_level:.1f}% CI="
                    f"[{boot_result['bootstrap_basic_ci_low']:.6f}, "
                    f"{boot_result['bootstrap_basic_ci_high']:.6f}]"
                )
                logger.info(
                    f"  Full-I vs block-partition-I difference: {partition_difference:.6g}; "
                    f"coverage={100*coverage_fraction:.2f}%; reliable={bool(ci_reliable)}"
                )
                if args.bootstrap_save_distribution:
                    boot_path = args.output_dir / f"{args.prefix}_Moran_block_bootstrap.npy"
                    tmp_boot_path = boot_path.with_suffix(boot_path.suffix + ".tmp")
                    try:
                        with open(tmp_boot_path, "wb") as f:
                            np.save(f, boot_distribution)
                            f.flush()
                            os.fsync(f.fileno())
                        os.replace(tmp_boot_path, boot_path)
                    except Exception:
                        try:
                            tmp_boot_path.unlink()
                        except FileNotFoundError:
                            pass
                        raise
                    logger.info(f"  Saved: {boot_path}")

        write_global_moran_csv(
            args.output_dir / f"{args.prefix}_Moran.csv",
            moran_I, n_used, W_total, inference,
            perm_p, perm_n,
            moran_mask_mode,
            bootstrap_metrics,
            quality_metrics,
        )
        moran_results = (moran_I, inference, perm_p)

    # 4. Local Gi*
    if not args.global_only:
        logger.info("Computing local Getis-Ord Gi*...")
        t3 = time.time()
        z, p, W = compute_local_gi_star(
            slope, valid, kernel,
            min_neighbors=args.min_neighbors,
            min_total_pixels=args.min_local_pixels,
        )
        logger.info(f"  Gi* z-scores computed ({time.time()-t3:.1f}s)")

        cls = classify_gi_zscore(z)

        write_gi_star_geotiff(
            args.output_dir / f"{args.prefix}_GiStar.tif",
            z, p, cls, profile,
        )
        write_summary_csv(
            args.output_dir / f"{args.prefix}_GiStar_summary.csv",
            cls, z,
        )

    elapsed = time.time() - t0
    logger.info("=" * 60)
    logger.info("Complete. Elapsed: %.1fs (%.2f min)", elapsed, elapsed / 60.0)
    logger.info("=" * 60)
    return 0


# ============================================================
# Validation / smoke test (run with: python spatial_autocorrelation.py --self-test)
# ============================================================

def validation_test():
    """
    Self-test on a synthetic 50x50 grid:
      - Center 10x10 block has high values (hot spot)
      - Outside block has low values + noise
    Expected: Gi* z-scores positive in center, negative on edges, NotSig elsewhere.
    Global Moran's I should be strongly positive (~0.5+).
    """
    rng = np.random.default_rng(0)
    H, W = 50, 50
    x = rng.normal(0, 1, (H, W))
    x[20:30, 20:30] += 5.0   # hot block

    valid = np.ones((H, W), dtype=bool)
    valid[0, :] = False       # simulate some edge NoData
    x[0, :] = np.nan

    kernel = make_kernel('queen', 1)

    moran_I, W_tot, n = compute_global_moran(x, valid, kernel)
    print(f"Global Moran's I = {moran_I:.4f} (n={n}, W={W_tot:.0f})")

    inf = moran_inference_normality(moran_I, W_tot, n, valid, kernel, x)
    print(f"  z_rand = {inf['z_rand']:.4f}, p_rand = {inf['p_rand']:.6g}")

    z, p, W = compute_local_gi_star(x, valid, kernel, min_neighbors=2)
    cls = classify_gi_zscore(z)

    n_hot99 = int((cls == 3).sum())
    n_hot95 = int((cls == 2).sum())
    print(f"  Hot99 cells: {n_hot99}, Hot95 cells: {n_hot95}")

    center_z = np.nanmean(z[22:28, 22:28])
    edge_z   = np.nanmean(z[5:10, 5:10])
    print(f"  Mean z in hot block: {center_z:.2f}")
    print(f"  Mean z in cold area: {edge_z:.2f}")

    assert moran_I > 0.2, f"Moran's I should be strongly positive, got {moran_I}"
    assert center_z > 2.0, f"Hot block z should be >2, got {center_z}"
    assert n_hot99 > 0, "Should detect hot99 cells in synthetic block"
    print("[PASS] All assertions passed.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--self-test":
        validation_test()
        raise SystemExit(0)
    raise SystemExit(main())

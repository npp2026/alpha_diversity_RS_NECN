#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Memory-efficient spatial block bootstrap for Global Moran's I.

The raster is partitioned into non-overlapping rectangular blocks.  Each block
is reduced once to sufficient statistics for Moran's I, after which bootstrap
replicates resample only those compact statistics.  Cross-block spatial links
are deliberately excluded so independently resampled blocks do not create
artificial adjacency.

Important statistical limitation
--------------------------------
The bootstrap distribution is generated for the within-block (partition)
estimator, whereas the reported point estimate normally uses all raster links.
Confidence intervals are therefore *recentered* from the partition estimator to
the full-raster Moran's I.  The partition/full difference is always reported and
must be small relative to the scientific effect size; otherwise a larger block
or a block-origin sensitivity analysis is required.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import partial
from typing import Dict, List, Tuple

import numpy as np


@dataclass(frozen=True)
class BlockSpec:
    row0: int
    row1: int
    col0: int
    col1: int


@dataclass
class BlockStatistics:
    n: np.ndarray
    sum_x: np.ndarray
    sum_x2: np.ndarray
    w: np.ndarray
    edge_xx: np.ndarray
    edge_xsum: np.ndarray
    row0: np.ndarray
    col0: np.ndarray
    height: np.ndarray
    width: np.ndarray
    n_total_candidate_blocks: int
    n_dropped_blocks: int
    edge_mode: str

    @property
    def n_blocks(self) -> int:
        return int(self.n.size)


# Thread workers receive shared NumPy array references through a partial
# function.  No module-level mutable state is used, so the module remains
# re-entrant when multiple analyses run concurrently in one Python process.

def kernel_offsets(kernel: np.ndarray) -> Tuple[Tuple[int, int], ...]:
    """Return directed non-zero offsets for a binary symmetric kernel."""
    k = np.asarray(kernel, dtype=np.float64)
    if k.ndim != 2 or k.shape[0] % 2 != 1 or k.shape[1] % 2 != 1:
        raise ValueError("Spatial weights kernel must be a 2-D odd-sized array")
    if not np.all(np.isfinite(k)):
        raise ValueError("Spatial weights kernel contains non-finite values")
    if not np.all((k == 0.0) | (k == 1.0)):
        raise ValueError("Block bootstrap currently supports binary weights only")
    if not np.array_equal(k, np.flip(k, axis=(0, 1))):
        raise ValueError("Block bootstrap requires a symmetric spatial kernel")

    cy, cx = k.shape[0] // 2, k.shape[1] // 2
    offsets: List[Tuple[int, int]] = []
    for iy in range(k.shape[0]):
        for ix in range(k.shape[1]):
            if iy == cy and ix == cx:
                continue
            if k[iy, ix] != 0:
                offsets.append((iy - cy, ix - cx))
    if not offsets:
        raise ValueError("Spatial weights kernel has no non-center neighbors")
    return tuple(offsets)


def _balanced_intervals(length: int, nominal_size: int) -> List[Tuple[int, int]]:
    """Cover an axis with near-equal intervals no smaller than nominal on average."""
    n_groups = max(1, length // nominal_size)
    base, remainder = divmod(length, n_groups)
    sizes = [base + 1] * remainder + [base] * (n_groups - remainder)
    intervals: List[Tuple[int, int]] = []
    start = 0
    for size in sizes:
        intervals.append((start, start + size))
        start += size
    assert start == length
    return intervals


def generate_block_specs(
    shape: Tuple[int, int],
    block_rows: int,
    block_cols: int,
    offset_rows: int = 0,
    offset_cols: int = 0,
    *,
    edge_mode: str = "balanced",
) -> List[BlockSpec]:
    """Generate a non-overlapping block partition.

    edge_mode:
      - ``balanced``: preserve the complete extent by distributing remainder
        rows/columns across blocks; dimensions differ by at most one cell.
      - ``drop``: retain only exact full-size blocks and discard edge strips.
      - ``partial``: include potentially tiny edge fragments. This is provided
        only for sensitivity analysis because equal block resampling can
        overweight those fragments.
    """
    h, w = map(int, shape)
    block_rows = int(block_rows)
    block_cols = int(block_cols)
    if h < 1 or w < 1:
        raise ValueError("Raster shape must be positive")
    if block_rows < 1 or block_cols < 1:
        raise ValueError("block_rows and block_cols must be >= 1")
    if edge_mode not in {"balanced", "drop", "partial"}:
        raise ValueError("edge_mode must be balanced, drop, or partial")

    offset_rows = int(offset_rows) % block_rows
    offset_cols = int(offset_cols) % block_cols
    specs: List[BlockSpec] = []

    if edge_mode == "balanced":
        if offset_rows != 0 or offset_cols != 0:
            raise ValueError(
                "balanced edge mode requires zero grid offset; use drop mode "
                "for block-origin sensitivity runs"
            )
        row_intervals = _balanced_intervals(h, block_rows)
        col_intervals = _balanced_intervals(w, block_cols)
        for r0, r1 in row_intervals:
            for c0, c1 in col_intervals:
                specs.append(BlockSpec(r0, r1, c0, c1))
    elif edge_mode == "partial":
        for r0_raw in range(-offset_rows, h, block_rows):
            r0 = max(0, r0_raw)
            r1 = min(h, r0_raw + block_rows)
            if r1 <= r0:
                continue
            for c0_raw in range(-offset_cols, w, block_cols):
                c0 = max(0, c0_raw)
                c1 = min(w, c0_raw + block_cols)
                if c1 <= c0:
                    continue
                specs.append(BlockSpec(r0, r1, c0, c1))
    else:  # drop
        for r0 in range(offset_rows, h - block_rows + 1, block_rows):
            for c0 in range(offset_cols, w - block_cols + 1, block_cols):
                specs.append(BlockSpec(r0, r0 + block_rows, c0, c0 + block_cols))

    if not specs:
        raise ValueError(
            "No candidate blocks fit the raster. Reduce block dimensions or "
            "change the bootstrap edge mode."
        )
    return specs


def _shifted_pair_slices(
    h: int, w: int, dy: int, dx: int
) -> Tuple[Tuple[slice, slice], Tuple[slice, slice]]:
    """Source/target slices for one directed neighbor offset inside a block."""
    if dy >= 0:
        src_r = slice(0, h - dy)
        dst_r = slice(dy, h)
    else:
        src_r = slice(-dy, h)
        dst_r = slice(0, h + dy)

    if dx >= 0:
        src_c = slice(0, w - dx)
        dst_c = slice(dx, w)
    else:
        src_c = slice(-dx, w)
        dst_c = slice(0, w + dx)
    return (src_r, src_c), (dst_r, dst_c)


def _one_block_statistics(
    spec: BlockSpec,
    *,
    source_x: np.ndarray,
    source_valid: np.ndarray,
    offsets: Tuple[Tuple[int, int], ...],
    min_valid_pixels: int,
    min_valid_fraction: float,
) -> Tuple[float, ...] | None:
    x = source_x[spec.row0:spec.row1, spec.col0:spec.col1]
    valid = source_valid[spec.row0:spec.row1, spec.col0:spec.col1]
    n = int(valid.sum())
    area = int(valid.size)
    if n < min_valid_pixels:
        return None
    if area and n / area < min_valid_fraction:
        return None

    vals = x[valid]
    sum_x = float(vals.sum(dtype=np.float64))
    sum_x2 = float(np.square(vals, dtype=np.float64).sum(dtype=np.float64))

    h, w = x.shape
    total_w = 0.0
    edge_xx = 0.0
    edge_xsum = 0.0
    for dy, dx in offsets:
        if abs(dy) >= h or abs(dx) >= w:
            continue
        src_sl, dst_sl = _shifted_pair_slices(h, w, dy, dx)
        pair_valid = valid[src_sl] & valid[dst_sl]
        n_pair = int(pair_valid.sum())
        if n_pair == 0:
            continue
        src_vals = x[src_sl][pair_valid]
        dst_vals = x[dst_sl][pair_valid]
        total_w += n_pair
        edge_xx += float(np.multiply(src_vals, dst_vals).sum(dtype=np.float64))
        edge_xsum += float((src_vals + dst_vals).sum(dtype=np.float64))

    if total_w <= 0:
        return None

    return (
        float(n), sum_x, sum_x2, total_w, edge_xx, edge_xsum,
        float(spec.row0), float(spec.col0),
        float(spec.row1 - spec.row0), float(spec.col1 - spec.col0),
    )


def compute_block_statistics(
    x: np.ndarray,
    valid: np.ndarray,
    kernel: np.ndarray,
    block_rows: int,
    block_cols: int | None = None,
    *,
    offset_rows: int = 0,
    offset_cols: int = 0,
    min_valid_pixels: int = 30,
    min_valid_fraction: float = 0.0,
    workers: int = 1,
    edge_mode: str = "balanced",
    require_min_blocks: int = 2,
) -> BlockStatistics:
    """Compress raster blocks to sufficient statistics for Moran bootstrap."""
    x = np.asarray(x)
    valid = np.asarray(valid, dtype=bool)
    if x.ndim != 2 or valid.ndim != 2:
        raise ValueError("x and valid must be 2-D arrays")
    if x.shape != valid.shape:
        raise ValueError("x and valid must have identical shapes")
    if np.any(valid & ~np.isfinite(x)):
        raise ValueError("valid=True occurs at non-finite x values")

    block_rows = int(block_rows)
    block_cols = int(block_rows if block_cols is None else block_cols)
    workers = max(1, int(workers))
    min_valid_pixels = int(min_valid_pixels)
    min_valid_fraction = float(min_valid_fraction)
    require_min_blocks = int(require_min_blocks)
    if require_min_blocks < 0:
        raise ValueError("require_min_blocks must be >= 0")
    if min_valid_pixels < 2:
        raise ValueError("min_valid_pixels must be >= 2")
    if not (0.0 <= min_valid_fraction <= 1.0):
        raise ValueError("min_valid_fraction must be between 0 and 1")

    specs = generate_block_specs(
        x.shape,
        block_rows,
        block_cols,
        offset_rows,
        offset_cols,
        edge_mode=edge_mode,
    )
    offsets = kernel_offsets(kernel)

    block_worker = partial(
        _one_block_statistics,
        source_x=x,
        source_valid=valid,
        offsets=offsets,
        min_valid_pixels=min_valid_pixels,
        min_valid_fraction=min_valid_fraction,
    )
    if workers == 1:
        raw = [block_worker(spec) for spec in specs]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            raw = list(pool.map(block_worker, specs, chunksize=16))

    kept = [row for row in raw if row is not None]
    if len(kept) < require_min_blocks:
        raise ValueError(
            f"Only {len(kept)} usable spatial blocks remain; at least "
            f"{require_min_blocks} are required"
        )

    a = (
        np.asarray(kept, dtype=np.float64)
        if kept else np.empty((0, 10), dtype=np.float64)
    )
    return BlockStatistics(
        n=a[:, 0],
        sum_x=a[:, 1],
        sum_x2=a[:, 2],
        w=a[:, 3],
        edge_xx=a[:, 4],
        edge_xsum=a[:, 5],
        row0=a[:, 6].astype(np.int64),
        col0=a[:, 7].astype(np.int64),
        height=a[:, 8].astype(np.int64),
        width=a[:, 9].astype(np.int64),
        n_total_candidate_blocks=len(specs),
        n_dropped_blocks=len(specs) - len(kept),
        edge_mode=edge_mode,
    )


def moran_from_sufficient_statistics(
    n: float,
    sum_x: float,
    sum_x2: float,
    w: float,
    edge_xx: float,
    edge_xsum: float,
) -> float:
    """Compute Moran's I from node and directed-edge sufficient statistics."""
    if n < 2 or w <= 0:
        return float("nan")
    mean = sum_x / n
    denominator = sum_x2 - 2.0 * mean * sum_x + n * mean * mean
    numerator = edge_xx - mean * edge_xsum + mean * mean * w
    if denominator <= 0 or not np.isfinite(denominator):
        return float("nan")
    return float((n / w) * (numerator / denominator))


def block_partition_moran(stats: BlockStatistics) -> float:
    """Observed Moran's I using retained pixels and within-block links only."""
    return moran_from_sufficient_statistics(
        float(stats.n.sum()),
        float(stats.sum_x.sum()),
        float(stats.sum_x2.sum()),
        float(stats.w.sum()),
        float(stats.edge_xx.sum()),
        float(stats.edge_xsum.sum()),
    )


def _bootstrap_batch(
    seeds: np.ndarray,
    *,
    n_values: np.ndarray,
    sum_x_values: np.ndarray,
    sum_x2_values: np.ndarray,
    w_values: np.ndarray,
    edge_xx_values: np.ndarray,
    edge_xsum_values: np.ndarray,
) -> np.ndarray:
    """Generate a deterministic batch of bootstrap replicates.

    Arrays are shared by reference between threads.  Keeping them as explicit
    arguments makes the function safe for concurrent independent analyses.
    """
    n_blocks = int(n_values.size)
    seeds = np.asarray(seeds, dtype=np.uint64)
    out = np.empty(seeds.size, dtype=np.float64)
    for i, replicate_seed in enumerate(seeds):
        rng = np.random.default_rng(int(replicate_seed))
        sampled = rng.integers(0, n_blocks, size=n_blocks, endpoint=False)
        out[i] = moran_from_sufficient_statistics(
            float(n_values[sampled].sum(dtype=np.float64)),
            float(sum_x_values[sampled].sum(dtype=np.float64)),
            float(sum_x2_values[sampled].sum(dtype=np.float64)),
            float(w_values[sampled].sum(dtype=np.float64)),
            float(edge_xx_values[sampled].sum(dtype=np.float64)),
            float(edge_xsum_values[sampled].sum(dtype=np.float64)),
        )
    return out


def spatial_block_bootstrap_moran(
    stats: BlockStatistics,
    observed_moran: float,
    *,
    n_boot: int = 999,
    ci_level: float = 0.95,
    seed: int = 42,
    workers: int = 1,
) -> Tuple[Dict[str, float], np.ndarray]:
    """Generate block-bootstrap replicates and recentered confidence intervals.

    Let ``I_part`` be the observed within-block estimator and ``I_full`` the
    full-raster estimate.  Bootstrap errors are ``I* - I_part``.  Inverting
    those errors around ``I_full`` gives the recentered basic interval:

        [I_full + I_part - q_high, I_full + I_part - q_low]

    The former implementation used ``2*I_full - q`` and therefore shifted the
    interval by ``I_full - I_part`` whenever block boundaries mattered.
    """
    n_boot = int(n_boot)
    workers = max(1, min(int(workers), max(n_boot, 1)))
    seed = int(seed)
    if n_boot < 2:
        raise ValueError("n_boot must be >= 2")
    if seed < 0:
        raise ValueError("seed must be non-negative")
    if not np.isfinite(observed_moran):
        raise ValueError("observed_moran must be finite")
    if not (0.0 < ci_level < 1.0):
        raise ValueError("ci_level must be between 0 and 1")
    if stats.n_blocks < 2:
        raise ValueError("At least two usable spatial blocks are required")

    counts = np.full(workers, n_boot // workers, dtype=int)
    counts[: n_boot % workers] += 1
    replicate_seeds = np.random.SeedSequence(seed).generate_state(
        n_boot, dtype=np.uint64
    )
    split_at = np.cumsum(counts)[:-1]
    tasks = [part for part in np.split(replicate_seeds, split_at) if part.size]

    bootstrap_worker = partial(
        _bootstrap_batch,
        n_values=stats.n,
        sum_x_values=stats.sum_x,
        sum_x2_values=stats.sum_x2,
        w_values=stats.w,
        edge_xx_values=stats.edge_xx,
        edge_xsum_values=stats.edge_xsum,
    )
    if workers == 1:
        pieces = [bootstrap_worker(tasks[0])]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pieces = list(pool.map(bootstrap_worker, tasks))

    boot = np.concatenate(pieces)
    boot = boot[np.isfinite(boot)]
    min_success = max(2, int(np.ceil(0.8 * n_boot)))
    if boot.size < min_success:
        raise RuntimeError(
            f"Only {boot.size}/{n_boot} finite bootstrap replicates were produced"
        )

    alpha = 1.0 - ci_level
    q_low, q_high = np.quantile(boot, [alpha / 2.0, 1.0 - alpha / 2.0])
    boot_mean = float(boot.mean())
    boot_se = float(boot.std(ddof=1))
    block_i = float(block_partition_moran(stats))
    if not np.isfinite(block_i):
        raise RuntimeError("Block-partition Moran's I is not finite")

    # Raw quantiles target the partition estimator. Recenter them to I_full for
    # an interpretable interval while preserving the bootstrap spread.
    recentered_pct_low = observed_moran + float(q_low) - block_i
    recentered_pct_high = observed_moran + float(q_high) - block_i
    recentered_basic_low = observed_moran + block_i - float(q_high)
    recentered_basic_high = observed_moran + block_i - float(q_low)

    result: Dict[str, float] = {
        "bootstrap_n_requested": float(n_boot),
        "bootstrap_n_success": float(boot.size),
        "bootstrap_ci_level": float(ci_level),
        "bootstrap_mean": boot_mean,
        "bootstrap_se": boot_se,
        "bootstrap_block_partition_I": block_i,
        "bootstrap_partition_difference_vs_full_I": block_i - observed_moran,
        "bootstrap_bias_vs_partition_I": boot_mean - block_i,
        # Retained for backward compatibility; it conflates partition and
        # resampling bias, so do not use it as the primary bias diagnostic.
        "bootstrap_bias_vs_full_I": boot_mean - observed_moran,
        "bootstrap_partition_percentile_ci_low": float(q_low),
        "bootstrap_partition_percentile_ci_high": float(q_high),
        "bootstrap_percentile_ci_low": float(recentered_pct_low),
        "bootstrap_percentile_ci_high": float(recentered_pct_high),
        "bootstrap_basic_ci_low": float(recentered_basic_low),
        "bootstrap_basic_ci_high": float(recentered_basic_high),
    }
    return result, boot


def block_statistics_metadata(
    stats: BlockStatistics,
    *,
    block_rows: int,
    block_cols: int,
    pixel_size_x: float,
    pixel_size_y: float,
    offset_rows: int,
    offset_cols: int,
    min_valid_pixels: int,
    min_valid_fraction: float,
) -> Dict[str, float | str]:
    """Return flat block diagnostics suitable for the metric/value CSV."""
    return {
        "bootstrap_block_rows_cells": float(block_rows),
        "bootstrap_block_cols_cells": float(block_cols),
        "bootstrap_block_height_map_units": float(block_rows * pixel_size_y),
        "bootstrap_block_width_map_units": float(block_cols * pixel_size_x),
        "bootstrap_grid_offset_rows_cells": float(offset_rows),
        "bootstrap_grid_offset_cols_cells": float(offset_cols),
        "bootstrap_block_edge_mode": stats.edge_mode,
        "bootstrap_partial_edge_blocks_included": int(stats.edge_mode == "partial"),
        "bootstrap_min_valid_pixels_per_block": float(min_valid_pixels),
        "bootstrap_min_valid_fraction_per_block": float(min_valid_fraction),
        "bootstrap_candidate_blocks": float(stats.n_total_candidate_blocks),
        "bootstrap_usable_blocks": float(stats.n_blocks),
        "bootstrap_dropped_blocks": float(stats.n_dropped_blocks),
        "bootstrap_total_valid_pixels_across_blocks": float(stats.n.sum()),
        "bootstrap_min_block_rows_observed": float(stats.height.min()) if stats.height.size else float("nan"),
        "bootstrap_max_block_rows_observed": float(stats.height.max()) if stats.height.size else float("nan"),
        "bootstrap_min_block_cols_observed": float(stats.width.min()) if stats.width.size else float("nan"),
        "bootstrap_max_block_cols_observed": float(stats.width.max()) if stats.width.size else float("nan"),
        "bootstrap_median_valid_pixels_per_block": float(np.median(stats.n)) if stats.n.size else float("nan"),
        "bootstrap_cv_valid_pixels_per_block": float(
            np.std(stats.n, ddof=1) / np.mean(stats.n)
            if stats.n.size > 1 and np.mean(stats.n) > 0 else 0.0
        ),
        "bootstrap_min_valid_pixels_observed": float(stats.n.min()) if stats.n.size else float("nan"),
        "bootstrap_max_valid_pixels_observed": float(stats.n.max()) if stats.n.size else float("nan"),
    }

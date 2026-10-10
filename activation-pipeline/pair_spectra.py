"""
Detector spectra per decay of each (volume, isotope) pair, one file per event class.

Written by build_spectra.py (from the step-4 hit lists) and read by steps 6-8
(accumulate_spectra.py, activation_history.py, average_spectrum.py).

File: <dir>/spectra_<mode>.npz, with
    volumes, isotopes : str arrays (P,)       pair names, as in active_isotopes.pkl
    ndecays           : int64 (P,)            simulated decays of each pair
    edges             : float64 (C+1,)        channel edges [keV]
    counts            : sparse CSR (P, C)     selected events per channel (integers, stored as
                                              data/indices/indptr/shape)
    meta              : JSON string           mode, selection parameters, source directory, ...

The spectrum of a pair is counts / ndecays / channel width, in counts/keV/decay.

Library usage:

    from pair_spectra import save_pair_spectra, load_pair_spectra
    lib = load_pair_spectra("result_postact/spectra_compton.npz")
    lib.centres                                   # channel centres [keV]
    lib.spectrum("PV-Absorber_006", "Na22")       # counts/keV/decay, or None if absent
    lib.rate_per_decay("PV-Absorber_006", "Na22") # selected events per decay (all channels)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp


def spectra_path(directory: str | Path, mode: str) -> Path:
    return Path(directory) / f"spectra_{mode}.npz"


def save_pair_spectra(
    path:     str | Path,
    volumes:  list[str] | np.ndarray,
    isotopes: list[str] | np.ndarray,
    ndecays:  np.ndarray,
    edges:    np.ndarray,
    counts:   sp.spmatrix | np.ndarray,
    meta:     dict,
) -> None:
    """Writes one spectra file (see the module docstring)."""
    counts = sp.csr_matrix(counts)
    n_pairs = len(volumes)
    if not (len(isotopes) == n_pairs == len(ndecays) == counts.shape[0]):
        raise ValueError("volumes, isotopes, ndecays and counts rows must have the same length")
    if counts.shape[1] != len(edges) - 1:
        raise ValueError("counts must have len(edges) - 1 columns")
    np.savez_compressed(
        path,
        volumes=np.asarray(volumes, dtype=str),
        isotopes=np.asarray(isotopes, dtype=str),
        ndecays=np.asarray(ndecays, dtype=np.int64),
        edges=np.asarray(edges, dtype=np.float64),
        data=counts.data, indices=counts.indices, indptr=counts.indptr,
        shape=np.asarray(counts.shape),
        meta=json.dumps(meta),
    )


class PairSpectra:
    """Spectra per decay of all pairs of one event class."""

    def __init__(self, path: str | Path):
        with np.load(path, allow_pickle=False) as f:
            self.volumes = f["volumes"]
            self.isotopes = f["isotopes"]
            self.ndecays = f["ndecays"]
            self.edges = f["edges"]
            self.counts = sp.csr_matrix((f["data"], f["indices"], f["indptr"]),
                                        shape=tuple(f["shape"]))
            self.meta = json.loads(str(f["meta"]))
        self.path = Path(path)
        self.widths = np.diff(self.edges)
        self.centres = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.index = {(v, i): k for k, (v, i) in enumerate(zip(self.volumes, self.isotopes))}

    @property
    def mode(self) -> str:
        return self.meta.get("mode", "")

    def __contains__(self, pair: tuple[str, str]) -> bool:
        return pair in self.index

    def spectrum(self, volume: str, isotope: str) -> np.ndarray | None:
        """Counts/keV/decay of a pair, or None if the pair is not in the file."""
        k = self.index.get((volume, isotope))
        if k is None:
            return None
        return self.counts[k].toarray().ravel() / self.ndecays[k] / self.widths

    def rate_per_decay(self, volume: str, isotope: str) -> float | None:
        """Selected events per decay of a pair (all channels), or None if absent."""
        k = self.index.get((volume, isotope))
        if k is None:
            return None
        return float(self.counts[k].sum()) / self.ndecays[k]


def load_pair_spectra(path: str | Path) -> PairSpectra:
    return PairSpectra(path)

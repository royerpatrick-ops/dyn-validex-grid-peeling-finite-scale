#!/usr/bin/env python3
"""Render the 64 prespecified translation-level scores from the saved JSON.

This script performs no draw, peeling, oracle, PDE, or geometric campaign.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "evidence" / "v0.6" / "DECISION_RESULTAT_v0_6.json"
OUTPUT = ROOT / "figures" / "primary_scores.png"


def main() -> None:
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    primary = payload["primary_result"]
    rows = primary["row_scores_with_exact_identity"]
    x = [int(row["row_index"]) + 1 for row in rows]
    y = [float(row["W"]) * 1e7 for row in rows]
    mean = float(primary["mean_W"]) * 1e7
    lower = float(primary["lower_confidence_bound_nominal"]) * 1e7

    fig, ax = plt.subplots(figsize=(8.0, 4.25), constrained_layout=True)
    ax.scatter(x, y, s=28, color="#245B78", edgecolor="white", linewidth=0.5,
               label="Score par translation")
    ax.axhline(mean, color="#245B78", linewidth=1.8,
               label=f"Moyenne = {mean:.4f}")
    ax.axhline(lower, color="#9B6B16", linewidth=1.8, linestyle="--",
               label=f"Borne nominale = {lower:.4f}")
    ax.set_xlabel("Translation appariée (indice gelé)")
    ax.set_ylabel(r"$W_j$  (unités de $10^{-7}$)")
    ax.set_title("64 scores primaires à R = 512")
    ax.set_xlim(0, 65)
    ax.set_ylim(0.9, 1.4)
    ax.grid(axis="y", color="#D6DEE3", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=300, metadata={"Author": "Patrick Royer"})
    plt.close(fig)


if __name__ == "__main__":
    main()

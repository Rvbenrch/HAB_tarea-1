#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rubenscript.py — Enriquecimiento funcional con g:Profiler (ORA)
Requisitos: gprofiler-official, pandas, matplotlib, numpy

Resumen del método:
- Tipo: Over-Representation Analysis (hipergeométrico) vía g:Profiler.
- Múltiples contrastes: FDR por defecto de g:Profiler (Benjamini–Hochberg).
- Anotaciones IEA: incluidas por defecto (opción conmutables).
- Filtrado: se aplica en local con umbral de FDR (--fdr).

Entrada:
- Archivo de texto con un gen por línea (símbolo oficial).

Salida:
- TSV con columnas principales (source, term_id, term_name, p_value, adjusted_p_value, ...).
- PNG opcional con los TOP términos por -log10(FDR).
- JSON con metadatos de la ejecución.

Uso:
  python rubenscript.py \
    --input data/genes_input.txt \
    --organism hsapiens \
    --sources GO:BP GO:MF GO:CC REAC \
    --fdr 0.05 \
    --top 20 \
    --outdir results

Notas:
- Valida organismos y fuentes (muestra advertencias si algo no es soportado).
- Requiere conexión a Internet para consultar g:Profiler.
"""

from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

try:
    from gprofiler import GProfiler
except Exception as e:
    print(
        "ERROR: No se pudo importar 'gprofiler-official'. "
        "Instala dependencias con: pip install -r requirements.txt",
        file=sys.stderr,
    )
    raise

# --- Configuración/constantes ---
KNOWN_SOURCES = {
    "GO:BP", "GO:MF", "GO:CC",
    "REAC", "WP", "HP", "MIRNA", "TF", "CORUM",
    "KEGG"
}
KNOWN_ORGANISMS = {
    "hsapiens", "mmusculus", "rnorvegicus", "drerio",
    "athaliana", "scerevisiae", "dmelanogaster"
}

def read_gene_list(path: Path) -> list[str]:
    """Lee un archivo de texto con un gen por línea, ignora comentarios '#' y filas vacías."""
    if not path.exists():
        raise FileNotFoundError(f"No existe el archivo de entrada: {path}")
    genes: list[str] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g.replace(" ", "").upper())
    if not genes:
        raise ValueError("La lista de genes está vacía.")
    return genes

def ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)

def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")

def validate_sources(sources: list[str]) -> tuple[list[str], list[str]]:
    """Separa fuentes válidas de desconocidas para informar al usuario."""
    valid = [s for s in sources if s in KNOWN_SOURCES]
    unknown = [s for s in sources if s not in KNOWN_SOURCES]
    return (valid if valid else sources, unknown)

def clip_neglog10(pvals: pd.Series) -> pd.Series:
    """Devuelve -log10(p) con protección ante p=0."""
    safe = np.clip(pvals.astype(float), 1e-300, 1.0)
    return -np.log10(safe)

def plot_top_terms(df: pd.DataFrame, out_png: Path, top: int = 20) -> None:
    """
    Barplot horizontal de TOP términos por -log10(FDR).
    Eje Y muestra: name (truncada) + [id] automáticamente,
    con fallback si term_name/term_id no existen.
    """
    if df.empty:
        print("No hay términos significativos para graficar.")
        return

    need_col = "adjusted_p_value" if "adjusted_p_value" in df.columns else "p_value"
    if need_col not in df.columns:
        print("No encontré columnas de p-valor; omito la figura.")
        return

    plot_df = df.sort_values(need_col, ascending=True).head(top).copy()
    plot_df["neglog10_fdr"] = clip_neglog10(plot_df[need_col])

    # DETECCIÓN de columnas de nombre/id
    # Buscamos los nombres válidos de columna (usar 'name' y 'id' si existen)
    name_col = None
    id_col = None
    for candidate in ['term_name', 'name', 'description']:
        if candidate in plot_df.columns:
            name_col = candidate
            break
    for candidate in ['term_id', 'id']:
        if candidate in plot_df.columns:
            id_col = candidate
            break

    if not name_col:
        print("Advertencia: No se encontró columna de nombres de términos. Usando 'source'.")
        name_col = 'source'
    if not id_col:
        print("Advertencia: No se encontró columna de IDs de términos. Usando 'source'.")
        id_col = 'source'

    def shorten(txt: str, n: int = 80) -> str:
        return (txt[: n - 1] + "…") if isinstance(txt, str) and len(txt) > n else str(txt)

    ylabels = [
        f"{shorten(str(row[name_col])) if pd.notnull(row[name_col]) and str(row[name_col]).strip() else '(sin nombre)'}"
        f" [{str(row[id_col]) if pd.notnull(row[id_col]) and str(row[id_col]).strip() else '(sin ID)'}]"
        for _, row in plot_df.iterrows()
    ]

    plt.figure(figsize=(12, max(4, 0.5 * len(plot_df))))
    plt.barh(range(len(plot_df)), plot_df["neglog10_fdr"])
    plt.yticks(range(len(plot_df)), ylabels, fontsize=10)
    plt.xlabel("-log10(FDR)" if need_col == "adjusted_p_value" else "-log10(p)")
    plt.ylabel("Término")
    plt.title("Top términos enriquecidos")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()



def run_enrichment(
    genes: list[str],
    organism: str,
    sources: list[str],
    user_threshold_fdr: float,
    include_iea: bool = True,
) -> pd.DataFrame:
    """Lanza g:Profiler (profile/ORA) y devuelve DataFrame filtrado por FDR."""
    gp = GProfiler(return_dataframe=True)
    res = gp.profile(
        organism=organism,
        query=genes,
        sources=sources if sources else None,
        no_iea=not include_iea,
        user_threshold=1.0,
    )
    if res is None or len(res) == 0:
        return pd.DataFrame()
    # Mantén estas columnas si existen
    # Luego, en run_enrichment sustituye:
    keep_cols = [c for c in [
        "source", "term_id", "term_name", "name", "id",
        "p_value", "adjusted_p_value",
        "intersection_size", "effective_domain_size",
    "   query_size", "intersections",
        "precision", "recall"
    ] if c in res.columns]
    res = res[keep_cols].copy()
    
    if "adjusted_p_value" in res.columns:
        res = res[res["adjusted_p_value"] <= user_threshold_fdr].copy()
    sort_cols = [c for c in ["adjusted_p_value", "p_value"] if c in res.columns]
    if sort_cols:
        res = res.sort_values(sort_cols, ascending=True)
    if "intersections" in res.columns:
        res["genes_hits"] = res["intersections"].apply(
            lambda x: ",".join(x) if isinstance(x, (list, tuple)) else ("" if pd.isna(x) else str(x))
        )
    return res

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Enriquecimiento funcional (g:Profiler) — CLI reproducible."
    )
    p.add_argument("--input", type=str, default="data/genes_input.txt",
                   help="Ruta a lista de genes (un gen por línea).")
    p.add_argument("--organism", type=str, default="hsapiens",
                   help=f"Organismo (p.ej. {', '.join(sorted(KNOWN_ORGANISMS))}).")
    p.add_argument("--sources", nargs="*", default=["GO:BP", "GO:MF", "GO:CC", "REAC"],
                   help="Bases (p.ej. GO:BP GO:MF GO:CC REAC WP HP KEGG...).")
    p.add_argument("--fdr", type=float, default=0.05,
                   help="Umbral de FDR (ajustado).")
    p.add_argument("--top", type=int, default=20,
                   help="Nº de términos a mostrar en la figura.")
    p.add_argument("--outdir", type=str, default="results",
                   help="Directorio de salida.")
    p.add_argument("--no-plot", action="store_true",
                   help="No generar figura (solo TSV).")
    p.add_argument("--no-iea", action="store_true",
                   help="Excluir anotaciones IEA (por defecto se incluyen).")
    return p.parse_args(argv)

def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    if args.organism not in KNOWN_ORGANISMS:
        print(f"AVISO: Organismo '{args.organism}' no verificado en lista común; "
              f"g:Profiler puede admitirlo igualmente.", file=sys.stderr)
    sources, unknown = validate_sources(args.sources)
    if unknown:
        print(f"AVISO: Se ignoran/validarán en g:Profiler fuentes no reconocidas: {', '.join(unknown)}",
              file=sys.stderr)

    try:
        genes = read_gene_list(input_path)
    except Exception as e:
        print(f"Error leyendo genes: {e}", file=sys.stderr)
        return 2

    print(f"Genes leídos: {len(genes)} — organismo: {args.organism}")
    print(f"Fuentes solicitadas: {', '.join(args.sources)} — FDR ≤ {args.fdr}")

    try:
        df = run_enrichment(
            genes=genes,
            organism=args.organism,
            sources=sources,
            user_threshold_fdr=args.fdr,
            include_iea=not args.no_iea,
        )
    except Exception as e:
        print(
            "Fallo al ejecutar g:Profiler.\n"
            "¿Tienes conexión a Internet? ¿Instalaste gprofiler-official?\n"
            f"Detalle: {e}",
            file=sys.stderr,
        )
        return 3

    out_base = f"enrichment_{timestamp()}"
    tsv_path = outdir / f"{out_base}.tsv"
    png_path = outdir / f"{out_base}.png"
    meta_path = outdir / f"{out_base}.json"

    if df.empty:
        print("No se encontraron términos con FDR ≤ umbral.")
        pd.DataFrame(columns=[
            "source", "term_id", "term_name", "p_value", "adjusted_p_value", "genes_hits"
        ]).to_csv(tsv_path, sep="\t", index=False)
        print(f"Guardado TSV vacío: {tsv_path}")
    else:
        df.to_csv(tsv_path, sep="\t", index=False)
        print(f"Resultados guardados en: {tsv_path}")

    meta = {
        "input": str(input_path),
        "n_genes": len(genes),
        "organism": args.organism,
        "sources_requested": args.sources,
        "sources_used": sources,
        "fdr_threshold": args.fdr,
        "include_iea": not args.no_iea,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "tsv": str(tsv_path),
        "png": (str(png_path) if not args.no_plot else None),
        "script": "rubenscript.py",
        "version_info": {
            "python": sys.version.split()[0]
        }
    }
    with meta_path.open("w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)
    print(f"Metadatos guardados en: {meta_path}")

    if not args.no_plot and not df.empty:
        try:
            plot_top_terms(df, png_path, top=args.top)
            if png_path.exists():
                print(f"Figura guardada en: {png_path}")
        except Exception as e:
            print(f"No se pudo generar la figura: {e}", file=sys.stderr)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())

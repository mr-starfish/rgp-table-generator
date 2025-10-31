# plot_rgp_isobars_only.py
# Compara APENAS isóbaras (Original vs Gerado)
# Visual melhorado + ranges por dataset/tabela (config centralizada)

import re
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import LogFormatterMathtext
from matplotlib.ticker import MultipleLocator
from typing import List, Dict, Any, Optional

# ========= CONFIG =========
FILE_ORIG = "RGP_Exemplo_Ideal.rgp"          # <-- ajuste

# Exemplos de seleção de arquivo gerado:
FILE_GEN  = "CO2.rgp";       OUTDIR = "./plots_isobars_only/CO2"
# FILE_GEN  = "METH_09_CO2_91_800_1000_CORRIGIDO.rgp"; OUTDIR = "./plots_isobars_only/METHCO2"

N_ISOBARS   = 16                 # quantas isóbaras subamostrar
USE_LOG_FOR = {4, 5}             # y-log para Cv(4) e Cp(5)
MARK_EVERY  = 8                  # intervalo de marcadores na curva azul (None desliga)

PROPERTY_NAMES = [
    "Entalpia [ J/kg ]",
    "Velocidade do Som [ m/s ]",
    "Volume Específico [ m³/kg ]",
    "Cv [ J/(kg·K) ]",
    "Cp [ J/(kg·K) ]",
    "dP/dv [ Pa·m³/kg ]",
    "Entropia [ J/(kg·K) ]",
    "Viscosidade Dinâmica [ kg/(m·s) ]",
    "Condutividade Térmica [ W/(m·K) ]"
]

# --- Perfil: auto / CO2 / METHCO2 ---
PROFILE = "auto"  # "auto" tenta deduzir pelo nome do FILE_GEN; force "CO2" ou "METHCO2" se preferir

# ========= RANGES CENTRALIZADOS =========
# Convenções:
# - Índice de tabela = 1..9
# - X_RANGES[tab] = {'min': float|None, 'max': float|None} ou None (sem ajuste)
# - Y_RANGES[tab] pode ser:
#     None                               -> sem ajuste
#     {'mode':'raw',    'min':..., 'max':...}          -> valores no unit real
#     {'mode':'scaled', 'min':..., 'max':..., 'factor':...}  -> valores no "número exibido" com multiplicador ×10^n
#
# Exemplos de 'scaled':
#   {"mode":"scaled","min":0,   "max":1.4, "factor":1e6}  -> aplica 0..1.4×10^6
#   {"mode":"scaled","min":0.5, "max":3.0, "factor":1e3}  -> aplica 0.5..3.0×10^3
#   {"mode":"scaled","min":0,   "max":2.5, "factor":1e-4} -> aplica 0..2.5×10^-4
#
RANGE_CONFIG: Dict[str, Dict[str, Dict[int, Optional[Dict[str, Any]]]]] = {
    "CO2": {
        "X_RANGES": {
            # "muda nada" -> todas sem ajuste; deixe vazio ou ponha None
            # 1: None, 2: None, ...
            1: {"min": 250, "max": 1000},
            2: {"min": 250, "max": 1000},
            3: {"min": 250, "max": 1000},
            4: {"min": 250, "max": 1000},
            5: {"min": 250, "max": 1000},
            6: {"min": 250, "max": 1000},
            7: {"min": 250, "max": 1000},
            8: {"min": 250, "max": 1000},
            9: {"min": 250, "max": 1000},
        },
        "Y_RANGES": {
            1: {"mode": "scaled", "min": 0.0, "max": 1.4,  "factor": 1e6},
            2: {"mode": "scaled", "min": 0.4, "max": 1,  "factor": 1e3},
            3: {"mode": "raw",    "min": 0.0025, "max": 0.018},
            4: {"mode": "raw",    "min": 850, "max": 1100},     # log
            5: {"mode": "raw",    "min": 1000,"max": 4000},     # log
            6: {"mode": "scaled", "min": -1.0,"max": 0.0,  "factor": 1e12},
            7: {"mode": "scaled", "min": 0.5, "max": 3.0,  "factor": 1e3},
            8: {"mode": "scaled", "min": 0.0, "max": 2.5,  "factor": 1e-4},
            9: {"mode": "raw",    "min": 0.05,"max": 0.18},
        },
    },
    "METHCO2": {
        "X_RANGES": {
            # TODAS AS TABELAS COM X MÁXIMO DE 700
            1: {"min": 250, "max": 700},
            2: {"min": 250, "max": 700},
            3: {"min": 250, "max": 700},
            4: {"min": 250, "max": 700},
            5: {"min": 250, "max": 700},
            6: {"min": 250, "max": 700},
            7: {"min": 250, "max": 700},
            8: {"min": 250, "max": 700},
            9: {"min": 250, "max": 700},
        },
        "Y_RANGES": {
            1: {"mode": "scaled", "min": 0.0, "max": 1.0,  "factor": 1e6},
            2: {"mode": "scaled", "min": 0.4, "max": 1,  "factor": 1e3},
            3: {"mode": "raw",    "min": 0.0025, "max": 0.018},
            4: {"mode": "raw",    "min": 850, "max": 1100},     # log
            5: {"mode": "raw",    "min": 1000,"max": 4000},     # log
            6: {"mode": "scaled", "min": -1.0,"max": 0.0,  "factor": 1e12},
            7: {"mode": "scaled", "min": 0.5, "max": 3.0,  "factor": 1e3},
            8: {"mode": "scaled", "min": 0.0, "max": 2.5,  "factor": 1e-4},
            9: {"mode": "raw",    "min": 0.05,"max": 0.18},
        },
    },
}

# ===== utilidades =====
def safe_filename(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s)

def parse_sci_tokens(line: str) -> List[float]:
    vals = []
    for t in line.strip().split():
        try:
            vals.append(float(t.replace('D','E')))
        except Exception:
            pass
    return vals

def parse_rgp_tables(filepath: str) -> List[Dict[str, Any]]:
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    tabs, i = [], 0
    while i < len(lines):
        m = re.match(r"\$TABLE_(\d+)", lines[i].strip())
        if not m:
            i += 1; continue
        num = int(m.group(1))

        i += 1
        ints = re.findall(r"[-+]?\d+", lines[i])
        if len(ints) < 2:
            raise RuntimeError(f"Sem t_steps/p_steps em {filepath}, linha {i}")
        t_steps, p_steps = int(ints[0]), int(ints[1])

        i += 1
        T = []
        for _ in range((t_steps + 4)//5):
            T += parse_sci_tokens(lines[i]); i += 1
        T = np.array(T[:t_steps], float)

        P = []
        for _ in range((p_steps + 4)//5):
            P += parse_sci_tokens(lines[i]); i += 1
        P = np.array(P[:p_steps], float)

        need = t_steps * p_steps
        flat = []
        while len(flat) < need:
            flat += parse_sci_tokens(lines[i]); i += 1
        Z = np.array(flat[:need], float).reshape((p_steps, t_steps)).T  # (t, p)

        Tsat = []
        while len(Tsat) < p_steps:
            Tsat += parse_sci_tokens(lines[i]); i += 1
        Phi  = []
        while len(Phi) < p_steps:
            Phi  += parse_sci_tokens(lines[i]); i += 1

        tabs.append(dict(num=num, T=T, P=P, Z=Z))

    tabs.sort(key=lambda d: d["num"])
    return tabs

def pick_indices_evenly(n: int, k: int) -> np.ndarray:
    if k >= n:
        return np.arange(n)
    return np.unique(np.linspace(0, n-1, num=k, dtype=int))

def _interp_column(T_src: np.ndarray, y_src: np.ndarray, T_dst: np.ndarray) -> np.ndarray:
    out = np.full_like(T_dst, np.nan, dtype=float)
    mask = (T_dst >= T_src.min()) & (T_dst <= T_src.max())
    out[mask] = np.interp(T_dst[mask], T_src, y_src)
    return out

def ensure_same_grid(tab1, tab2):
    T1, P1, Z1 = tab1["T"], tab1["P"], tab1["Z"]
    T2, P2, Z2 = tab2["T"], tab2["P"], tab2["Z"]

    if np.allclose(T1, T2) and np.allclose(P1, P2) and Z1.shape == Z2.shape:
        return T1, P1, Z1, Z2

    cols = min(Z1.shape[1], Z2.shape[1])
    Z2i = np.empty_like(Z1)
    for j in range(cols):
        Z2i[:, j] = _interp_column(T2, Z2[:, j], T1)
    if cols < Z1.shape[1]:
        Z2i[:, cols:] = np.nan
    return T1, P1, Z1, Z2i

def detect_profile(file_gen: str) -> str:
    if PROFILE in ("CO2", "METHCO2"):
        return PROFILE
    base = os.path.basename(file_gen).upper()
    if "METH" in base:
        return "METHCO2"
    if "CO2" in base:
        return "CO2"
    # fallback: CO2 (ou crie um perfil "DEFAULT")
    return "CO2"

def get_legend_labels(profile_key: str):
    # Padrão: mantém "Original" e "Gerado"
    if profile_key == "METHCO2":
        return ("CO2", "METHCO2")   # Original=CO2, Gerado=METHCO2
    return ("Original", "Gerado")


def apply_axis_ranges(ax, T: np.ndarray, table_idx: int, profile_key: str, use_log: bool):
    cfg = RANGE_CONFIG.get(profile_key, {})
    xcfg = (cfg.get("X_RANGES") or {}).get(table_idx)
    ycfg = (cfg.get("Y_RANGES") or {}).get(table_idx)

    # X range
    if xcfg:
        xmin = xcfg.get("min", None)
        xmax = xcfg.get("max", None)
        # valores None deixam automático de um dos lados
        cur_xmin, cur_xmax = ax.get_xlim()
        if xmin is None: xmin = cur_xmin
        if xmax is None: xmax = cur_xmax
        ax.set_xlim(xmin, xmax)

    # Y range
    if ycfg:
        mode = ycfg.get("mode", "raw")
        if mode == "raw":
            ymin, ymax = ycfg["min"], ycfg["max"]
        elif mode == "scaled":
            factor = ycfg["factor"]
            ymin = ycfg["min"] * factor
            ymax = ycfg["max"] * factor
        else:
            ymin = ymax = None

        if ymin is not None and ymax is not None:
            # Para log, garanta que limites sejam > 0 (exceto se você quer incluir <=0)
            if use_log and (ymin <= 0 or ymax <= 0):
                # evita erro de escala log
                ymin = max(ymin, np.finfo(float).tiny)
                ymax = max(ymax, ymin * 10.0)
            ax.set_ylim(ymin, ymax)

def plot_isobars_only(tab_orig, tab_gen, title_base, use_log, outdir, table_idx, profile_key):
    # ====== Ajustes visuais globais por figura ======
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 14,
        'axes.labelsize': 14,
        'legend.fontsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'axes.labelpad': 10,
    })
    plt.figure(figsize=(11, 7))

    T, P, Z1, Z2 = ensure_same_grid(tab_orig, tab_gen)

    # subset de isóbaras
    idxP = pick_indices_evenly(len(P), N_ISOBARS)

    # isóbaras (fundo)
    for j in idxP:
        plt.plot(T, Z1[:, j], color=(0.5, 0, 0.8, 0.12), linewidth=1.0, zorder=1)
        plt.plot(T, Z2[:, j], color=(0.5, 0, 0.8, 0.12), linewidth=1.0, zorder=1)
    proxy_isobaras, = plt.plot([], [], color=(0.5, 0, 0.8, 0.3), linewidth=3)

    # curvas médias
    y1 = np.nanmean(Z1[:, idxP], axis=1)
    y2 = np.nanmean(Z2[:, idxP], axis=1)

    # vermelho com contorno branco
    line1, = plt.plot(
        T, y1, color='red', linewidth=2.8, zorder=3,
        path_effects=[pe.Stroke(linewidth=4.6, foreground='white'), pe.Normal()],
        label="Original"
    )
    # azul com marcadores
    line2, = plt.plot(
        T, y2, color='blue', linewidth=1.6, zorder=4,
        marker=('o' if MARK_EVERY else None), markevery=(MARK_EVERY if MARK_EVERY else None),
        markersize=3.8, label="Gerado"
    )

    ax = plt.gca()
    if use_log:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    else:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3))
        offy = ax.yaxis.get_offset_text()
        offy.set_fontsize(14)
        offy.set_fontweight('bold')
        # offy.set_x(-0.03)

    # diferença máxima
    dmax = float(np.nanmax(np.abs(y2 - y1)))

    # Extrai número e propriedade
    if " - " in title_base:
        table_label, prop_label = title_base.split(" - ", 1)
    else:
        table_label, prop_label = title_base, title_base

    # Título enxuto
    plt.title(f"{table_label} — Isobáricas (subset={len(idxP)}/{len(P)})  |  Δmáx(médias)={dmax:.6g}")

    # Eixos
    plt.xlabel("Temperatura [K]")
    plt.ylabel(prop_label)

    # === aplica ranges por perfil/tabela ===
    apply_axis_ranges(ax, T, table_idx, profile_key, use_log)
    
    # --- Ticks do eixo X a cada 50 unidades ---
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.ticklabel_format(axis='x', style='plain')
    
    plt.grid(True, alpha=0.35, linewidth=0.8)
    lbl_orig, lbl_gen = get_legend_labels(profile_key)
    plt.legend(handles=[line1, line2, proxy_isobaras],
           labels=[lbl_orig, lbl_gen, "Isobáricas"],
           loc="best", frameon=True)


    plt.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    fname = os.path.join(outdir, safe_filename(f"{title_base}_isobars.png"))
    plt.savefig(fname, dpi=180)
    plt.close()
    return fname, dmax

# ===== Main =====
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tabs_o = parse_rgp_tables(FILE_ORIG)
    tabs_g = parse_rgp_tables(FILE_GEN)

    if len(tabs_o) < 9 or len(tabs_g) < 9:
        print(f"[ERRO] Esperava 9 tabelas em cada arquivo. Achei {len(tabs_o)} e {len(tabs_g)}.")
        return

    profile_key = detect_profile(FILE_GEN)
    print(f"Profile ativo: {profile_key}")
    print(f"Saída: {OUTDIR}\n")

    for k in range(9):
        name  = PROPERTY_NAMES[k] if k < len(PROPERTY_NAMES) else f"Table_{k+1}"
        title = f"Tabela {k+1} - {name}"
        use_log = (k+1) in USE_LOG_FOR
        fn, dmax = plot_isobars_only(tabs_o[k], tabs_g[k], title, use_log, OUTDIR, k+1, profile_key)
        print(f"- {os.path.basename(fn)}   Δmáx(médias)={dmax:.6g}")

if __name__ == "__main__":
    main()

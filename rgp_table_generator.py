"""
Gerador de tabelas .rgp para ANSYS CFX usando CoolProp.
Inclui UX melhorada: barra de progresso, avisos coloridos,
resumo final, validação de fluido/mistura e **retomada de onde parou**.
"""
import io
import os
import math
import logging
import time
from typing import List, Tuple, Optional, Dict

from colorama import Fore, Style, init as colorama_init
from CoolProp.CoolProp import PropsSI
from tqdm import tqdm

colorama_init(autoreset=True)


# =========================
# Utilitários gerais
# =========================
class WarningCounter(logging.Handler):
    """Contador de warnings para monitorar valores NaN durante geração."""
    def __init__(self):
        super().__init__()
        self.warning_count = 0
        self.warnings_por_tabela = {}

    def emit(self, record):
        if record.levelno == logging.WARNING:
            self.warning_count += 1
            msg = record.getMessage()
            if "[Tabela " in msg:
                tabela_num = msg.split("[Tabela ")[1].split("]")[0]
                self.warnings_por_tabela[tabela_num] = self.warnings_por_tabela.get(tabela_num, 0) + 1
            elif "[SAT " in msg:
                tabela_num = f"SAT_{msg.split('[SAT ')[1].split(']')[0]}"
                self.warnings_por_tabela[tabela_num] = self.warnings_por_tabela.get(tabela_num, 0) + 1
        #print(self.format(record))


class _NullWriter:
    def write(self, *args, **kwargs):
        return 0


def format_value(x: float, decimals: int = 7) -> str:
    """Formata valor float para notação científica compatível com CFX."""
    if (
        x is None or not isinstance(x, (int, float)) or x != x
        or x == float('inf') or x == float('-inf')
    ):
        return "0.0000000E+000"
    s = f"{x:.{decimals}E}"
    mantissa, exp = s.split('E')
    sign = exp[0]
    digits = exp[1:].rjust(3, '0')
    return f"{mantissa}E{sign}{digits}"


def get_mixture_properties(mixture):
    """
    Retorna propriedades críticas e constantes do fluido/mistura.
    """
    t_crit = PropsSI("Tcrit", mixture)
    p_crit = PropsSI("pcrit", mixture)
    t_triple = PropsSI("Ttriple", mixture)
    p_triple = PropsSI("ptriple", mixture)
    molar_mass = PropsSI("M", mixture)
    universal_gas_constant = PropsSI("GAS_CONSTANT", mixture)
    r_especifico = universal_gas_constant / molar_mass
    return t_crit, p_crit, t_triple, p_triple, r_especifico


def write_header(file_obj, table_name, t_steps, p_steps, vals):
    """
    Escreve o HEADER e DATA do arquivo .rgp.
    """
    for block in ['HEADER', 'DATA']:
        file_obj.write(f'$$$${block} \n')
        file_obj.write(f'$$${table_name} \n\t1 \n$$PARAM \n\t26 \n')
        file_obj.write('DESCRIPTION \n')
        file_obj.write(f'{table_name} [HelmholtzEOSBackend, CoolProp 6.8.0] \n')
        file_obj.write('NAME \n')
        file_obj.write(f'{table_name} \n')
        file_obj.write('INDEX \n')
        file_obj.write(f'{table_name} \n')
        file_obj.write('MODEL \n\t3 \nUNITS \n\t1 \n')
        file_obj.write('PMIN_SUPERHEAT \n\t' + format_value(vals["pmin_superheat"]) + ' \n')
        file_obj.write('PMAX_SUPERHEAT \n\t' + format_value(vals["pmax_superheat"]) + ' \n')
        file_obj.write('TMIN_SUPERHEAT \n\t' + format_value(vals["tmin_superheat"]) + ' \n')
        file_obj.write('TMAX_SUPERHEAT \n\t' + format_value(vals["tmax_superheat"]) + ' \n')
        file_obj.write('TMIN_SATURATION \n\t' + format_value(vals["tmin_saturation"]) + ' \n')
        file_obj.write('TMAX_SATURATION \n\t' + format_value(vals["tmax_saturation"]) + ' \n')
        file_obj.write('SUPERCOOLING \n\t' + format_value(vals["supercooling"]) + ' \n')
        file_obj.write('P_CRTICAL \n\t' + format_value(vals["p_critical"]) + ' \n')
        file_obj.write('P_TRIPLE \n\t' + format_value(vals["p_triple"]) + ' \n')
        file_obj.write('T_CRTICAL \n\t' + format_value(vals["t_critical"]) + ' \n')
        file_obj.write('T_TRIPLE \n\t' + format_value(vals["t_triple"]) + '\n')
        file_obj.write('GAS_CONSTANT \n\t' + format_value(vals["gas_constant"]) + ' \n')
        for n in range(1, 10):
            file_obj.write(f'TABLE_{n} \n\t{t_steps}\t{p_steps} \n')
    file_obj.write('$$SUPER_TABLE \n\t9 \n')


def write_axis_headers(file_obj, t_vals, p_vals):
    """
    Escreve os cabeçalhos de T e P em blocos de 5 colunas.
    """
    for i in range(0, len(t_vals), 5):
        file_obj.write(
            '\t' + '\t'.join(format_value(t_vals[j]) for j in range(i, min(i + 5, len(t_vals)))) + ' \n'
        )
    for i in range(0, len(p_vals), 5):
        file_obj.write(
            '\t' + '\t'.join(format_value(p_vals[j]) for j in range(i, min(i + 5, len(p_vals)))) + ' \n'
        )


# =========================
# PROTEÇÃO DE FASE (SAFE)
# =========================
def make_get_prop_safe(
    mixture, tcrit, pcrit, t_step,
    eps_min=1e-4, eps_max=2.0, alpha=0.02, beta=0.2,
    dt_guard=0.0, crit_p_frac=0.90, crit_t_frac=0.96, eps_crit=0.5
):
    """
    Wrapper que evita estados ambíguos na vizinhança da saturação com ε ADAPTATIVO.
    """
    def wrap(symbol):
        def _f(t, p):
            t_use = t
            try:
                if (p < pcrit) and (t < tcrit):
                    try:
                        tsat = PropsSI('T', 'P', p, 'Q', 1, mixture)
                        if t <= tsat + dt_guard:
                            delta = tsat - t  # >= 0
                            eps = max(eps_min, alpha * delta, beta * t_step)
                            if (p > crit_p_frac * pcrit) or (tsat > crit_t_frac * tcrit):
                                eps = max(eps, eps_crit)
                            eps = min(eps, eps_max)
                            t_use = tsat + eps
                    except Exception:
                        t_use = t

                if symbol == '1/D':
                    d = PropsSI('D', 'T', t_use, 'P', p, mixture)
                    return 1.0 / d
                elif symbol == 'D':
                    return PropsSI('D', 'T', t_use, 'P', p, mixture)
                elif symbol == 'd(P)/d(D)|T':
                    return PropsSI('d(P)/d(D)|T', 'T', t_use, 'P', p, mixture)
                else:
                    return PropsSI(symbol, 'T', t_use, 'P', p, mixture)

            except Exception:
                try:
                    tsat = PropsSI('T', 'P', p, 'Q', 1, mixture)
                    t_use = tsat + max(eps_min, eps_crit)
                    if symbol == '1/D':
                        d = PropsSI('D', 'T', t_use, 'P', p, mixture)
                        return 1.0 / d
                    elif symbol == 'D':
                        return PropsSI('D', 'T', t_use, 'P', p, mixture)
                    elif symbol == 'd(P)/d(D)|T':
                        return PropsSI('d(P)/d(D)|T', 'T', t_use, 'P', p, mixture)
                    else:
                        return PropsSI(symbol, 'T', t_use, 'P', p, mixture)
                except Exception:
                    return float('nan')
        return _f
    return wrap


def generate_property_map(mixture, tcrit, pcrit, t_step):
    """Mapeia número da tabela para a função de propriedade (com proteção de fase)."""
    safe = make_get_prop_safe(mixture, tcrit, pcrit, t_step)
    return {
        1: ('H',      safe('H')),
        2: ('A',      safe('A')),
        3: ('1/D',    safe('1/D')),
        4: ('Cv',     safe('O')),                                  # Cv
        5: ('Cp',     safe('C')),                                  # Cp
        6: ('dP/dv',  (lambda t, p: -safe('D')(t, p)**2 * safe('d(P)/d(D)|T')(t, p))),
        7: ('S',      safe('S')),
        8: ('μ',      safe('V')),                                  # viscosidade dinâmica
        9: ('λ',      safe('L')),                                  # condutividade térmica
    }


# ======================================================
# LIMPEZA / CORREÇÃO DE “BURACOS” NAS TABELAS
# ======================================================

def _is_bad(v: float, zero_tol: float) -> bool:
    return (not math.isfinite(v)) or (abs(v) <= zero_tol)


def _detect_constant_runs(row: List[float], min_run: int, eps_equal: float) -> List[Tuple[int, int]]:
    runs = []
    n = len(row)
    i = 0
    while i < n:
        j = i + 1
        while j < n and abs(row[j] - row[i]) <= eps_equal:
            j += 1
        run_len = j - i
        if run_len >= min_run:
            runs.append((i, j - 1))
        i = j
    return runs


def _mask_bad_and_constant(row: List[float], zero_tol: float, min_const_run: int, eps_equal: float) -> List[bool]:
    n = len(row)
    mask = [_is_bad(v, zero_tol) for v in row]
    for a, b in _detect_constant_runs(row, min_const_run, eps_equal):
        for k in range(a, b + 1):
            mask[k] = True
    return mask


def _first_valid(mask: List[bool]) -> Optional[int]:
    for i, bad in enumerate(mask):
        if not bad:
            return i
    return None


def _last_valid(mask: List[bool]) -> Optional[int]:
    for i in range(len(mask) - 1, -1, -1):
        if not mask[i]:
            return i
    return None


def _safe_slope(y0: float, y1: float, x0: int, x1: int) -> float:
    if x1 == x0:
        return 0.0
    return (y1 - y0) / (x1 - x0)


def _mono_clamp(y: float, y_left: float, y_right: float, increasing: Optional[bool]) -> float:
    if increasing is None:
        return y
    if increasing:
        lo = min(y_left, y_right)
        hi = max(y_left, y_right)
        return min(max(y, lo), hi)
    else:
        hi = max(y_left, y_right)
        lo = min(y_left, y_right)
        return min(max(y, lo), hi)


def _infer_monotonicity(row: List[float], mask: List[bool]) -> Optional[bool]:
    last = None
    up = down = 0
    for i, bad in enumerate(mask):
        if bad:
            continue
        if last is None:
            last = row[i]
            continue
        if row[i] > last:
            up += 1
        elif row[i] < last:
            down += 1
        last = row[i]
    if up == 0 and down == 0:
        return None
    return True if up >= down else False


def _fix_row(row: List[float],
             zero_tol: float,
             min_const_run: int,
             eps_equal: float) -> List[float]:
    n = len(row)
    out = row[:]
    mask = _mask_bad_and_constant(out, zero_tol, min_const_run, eps_equal)
    if all(mask):
        return out
    mono = _infer_monotonicity(out, mask)
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j < n and mask[j]:
            j += 1
        left = _last_valid(mask[:i+1])
        right = _first_valid(mask[j:])
        if right is not None:
            right = j + right
        if left is None and right is None:
            i = j
            continue
        elif left is None:
            r0 = right
            r1 = None
            k = right + 1
            while k < n:
                if not mask[k]:
                    r1 = k
                    break
                k += 1
            if r1 is None:
                for k in range(i, j):
                    out[k] = out[r0]
                    out[k] = _mono_clamp(out[k], out[r0], out[r0], mono)
            else:
                slope = _safe_slope(out[r0], out[r1], r0, r1)
                for k in range(i, j):
                    out[k] = out[r0] - slope * (r0 - k)
                    out[k] = _mono_clamp(out[k], out[r0], out[r1], mono)
        elif right is None:
            l0 = left
            l1 = None
            k = left - 1
            while k >= 0:
                if not mask[k]:
                    l1 = k
                    break
                k -= 1
            if l1 is None:
                for k in range(i, j):
                    out[k] = out[l0]
                    out[k] = _mono_clamp(out[k], out[l0], out[l0], mono)
            else:
                slope = _safe_slope(out[l1], out[l0], l1, l0)
                for k in range(i, j):
                    out[k] = out[l0] + slope * (k - l0)
                    out[k] = _mono_clamp(out[k], out[l1], out[l0], mono)
        else:
            y0, x0 = out[left], left
            y1, x1 = out[right], right
            slope = _safe_slope(y0, y1, x0, x1)
            for k in range(i, j):
                out[k] = y0 + slope * (k - x0)
                out[k] = _mono_clamp(out[k], y0, y1, mono)
        for k in range(i, j):
            mask[k] = False
        i = j
    return out


def _reshape_stream_to_rows(values: List[float], row_len: int) -> List[List[float]]:
    return [values[i:i+row_len] for i in range(0, len(values), row_len)]


def _flatten_rows(rows: List[List[float]]) -> List[float]:
    return [v for r in rows for v in r]


def fix_grid(values_flat: List[float], n_t: int,
             zero_tol: float = 1e-9, min_const_run: int = 2, eps_equal: float = 1e-12) -> List[float]:
    """
    Corrige um grid achatado de tamanho n_p * n_t (ordem: p por fora, t por dentro),
    linha lógica = sequência de n_t ao longo de T para cada P.
    """
    rows = _reshape_stream_to_rows(values_flat, n_t)
    fixed_rows = [_fix_row(r, zero_tol, min_const_run, eps_equal) for r in rows]
    return _flatten_rows(fixed_rows)


def fix_vector(vec: List[float],
               zero_tol: float = 1e-9, min_const_run: int = 2, eps_equal: float = 1e-12) -> List[float]:
    """Corrige um vetor 1D (ex.: Tsat(p) ou phi_sat(p))."""
    return _fix_row(vec, zero_tol, min_const_run, eps_equal)


# ======================================================
# ESCRITA DOS BLOCOS (ORIGINAL)
# ======================================================

def write_table_block(
    file_obj, table_num, t_steps, p_steps, t_vals, p_vals, get_prop,
    table_progress=None, value_progress=None, collector=None
):
    """
    Grava exatamente como no RGP original: primeiro cabeçalhos T e P,
    depois o GRID achatado (p varre por fora, T por dentro).
    Se 'collector' for dict, guarda os valores numéricos brutos em collector['tables'][table_num].
    """
    file_obj.write(f'$TABLE_{table_num} \n\t{t_steps}\t{p_steps} \n')
    write_axis_headers(file_obj, t_vals, p_vals)

    flat_vals_numbers = []
    p_iter = p_vals if table_progress is None else table_progress(
        p_vals, desc=f"Tabela {table_num}: Pressões"
    )
    for p in p_iter:
        t_iter = t_vals if value_progress is None else value_progress(
            t_vals, desc=f"Tabela {table_num} - P={p:.2f}"
        )
        for t in t_iter:
            try:
                v = get_prop(t, p)
            except Exception as exc:
                logging.warning(
                    Fore.YELLOW + f"[Tabela {table_num}] Erro em P={p:.2f}, T={t:.2f}: {exc}" +
                    Style.RESET_ALL
                )
                v = float('nan')
            flat_vals_numbers.append(v)

    if isinstance(collector, dict):
        collector.setdefault('tables', {})[table_num] = list(flat_vals_numbers)

    flat_vals = [format_value(v) for v in flat_vals_numbers]
    for i in range(0, len(flat_vals), 5):
        file_obj.write('\t' + '\t'.join(flat_vals[i:i+5]) + ' \n')


def write_saturation_block(file_obj, tbl_num, p_vals, mixture, t_triple, t_crit, collector=None):
    """
    Gera linha de Tsat e phi_sat para cada tabela.
    Se 'collector' for dict, guarda em:
      collector['sat_t'][tbl_num] = lista Tsat
      collector['sat_phi'][tbl_num] = lista phi_sat
    """
    t_media = t_crit
    tsat_list = []
    for p in p_vals:
        try:
            tsat = PropsSI('T', 'P', p, 'Q', 1, mixture)
        except Exception:
            try:
                tsat = t_media
            except Exception:
                tsat = 0.0
        tsat_list.append(tsat)

    for i in range(0, len(tsat_list), 5):
        file_obj.write(
            '\t' + '\t'.join(format_value(tsat_list[j]) for j in range(i, min(i + 5, len(tsat_list)))) + '\n'
        )

    property_keys = {
        1: ('H', 'H'), 2: ('A', 'A'), 3: ('1/D', 'D'), 4: ('Cv', 'O'), 5: ('Cp', 'C'),
        6: ('dP/dv', 'd(P)/d(D)|T'), 7: ('S', 'S'), 8: ('μ', 'V'), 9: ('λ', 'L')
    }
    phi_key, prop_symbol = property_keys[tbl_num]

    phi_sat_list = []
    for p in p_vals:
        try:
            if phi_key == '1/D':
                d_sat = PropsSI('D', 'P', p, 'Q', 1, mixture)
                phi_sat = 1 / d_sat
            elif phi_key == 'dP/dv':
                d_sat = PropsSI('D', 'P', p, 'Q', 1, mixture)
                dpdd_sat = PropsSI('d(P)/d(D)|T', 'P', p, 'Q', 1, mixture)
                phi_sat = -d_sat ** 2 * dpdd_sat
            else:
                phi_sat = PropsSI(prop_symbol, 'P', p, 'Q', 1, mixture)
        except Exception as exc:
            try:
                if phi_key == '1/D':
                    d_sat = PropsSI('D', 'T', t_media, 'P', p, mixture)
                    phi_sat = 1 / d_sat
                elif phi_key == 'dP/dv':
                    d_sat = PropsSI('D', 'T', t_media, 'P', p, mixture)
                    dpdd_sat = PropsSI('d(P)/d(D)|T', 'T', t_media, 'P', p, mixture)
                    phi_sat = -d_sat ** 2 * dpdd_sat
                else:
                    phi_sat = PropsSI(prop_symbol, 'T', t_media, 'P', p, mixture)
            except Exception:
                phi_sat = 0.0
            logging.warning(
                Fore.YELLOW + f"[SAT {tbl_num}] Erro em P={p:.2f} Pa: {exc}" + Style.RESET_ALL
            )
        phi_sat_list.append(phi_sat)

    for i in range(0, len(phi_sat_list), 5):
        file_obj.write(
            '\t' + '\t'.join(format_value(phi_sat_list[j]) for j in range(i, min(i + 5, len(phi_sat_list)))) + '\n'
        )

    if isinstance(collector, dict):
        collector.setdefault('sat_t', {})[tbl_num] = list(tsat_list)
        collector.setdefault('sat_phi', {})[tbl_num] = list(phi_sat_list)


# ======================================================
# ESCRITA DE ARQUIVO CORRIGIDO (REWRITE)
# ======================================================

def write_table_block_from_values(file_obj, table_num, t_steps, p_steps, t_vals, p_vals, flat_values):
    file_obj.write(f'$TABLE_{table_num} \n\t{t_steps}\t{p_steps} \n')
    write_axis_headers(file_obj, t_vals, p_vals)
    strings = [format_value(v) for v in flat_values]
    for i in range(0, len(strings), 5):
        file_obj.write('\t' + '\t'.join(strings[i:i+5]) + ' \n')


def write_saturation_block_from_values(file_obj, tsat_list, phi_list):
    for i in range(0, len(tsat_list), 5):
        file_obj.write('\t' + '\t'.join(format_value(tsat_list[j]) for j in range(i, min(i + 5, len(tsat_list)))) + '\n')
    for i in range(0, len(phi_list), 5):
        file_obj.write('\t' + '\t'.join(format_value(phi_list[j]) for j in range(i, min(i + 5, len(phi_list)))) + '\n')


# ======================================================
# SUPORTE A RETOMADA
# ======================================================

def _header_string(table_name, t_steps, p_steps, vals):
    """Gera o HEADER+DATA exatamente como o arquivo, para comparação byte a byte."""
    buf = io.StringIO()
    write_header(buf, table_name, t_steps, p_steps, vals)
    return buf.getvalue()


def _expected_lines_per_table(t_steps: int, p_steps: int) -> int:
    # 1 linha $TABLE_n, 1 linha com t_steps/p_steps
    t_hdr = math.ceil(t_steps / 5)
    p_hdr = math.ceil(p_steps / 5)
    grid = math.ceil((t_steps * p_steps) / 5)
    sat = 2 * math.ceil(p_steps / 5)
    return 2 + t_hdr + p_hdr + grid + sat


def _detect_resume(file_path: str, header_str: str, t_steps: int, p_steps: int):
    """
    Retorna (start_table:int, truncate_pos:int) se o arquivo tiver o mesmo cabeçalho
    e houver tabela incompleta. Caso contrário, retorna (1, None) para começar do zero,
    ou (10, None) se todas as 9 tabelas já estiverem completas.
    """
    if not os.path.exists(file_path):
        return 1, None
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    if not content:
        return 1, None
    if not content.startswith(header_str):
        return 1, None  # cabeçalho diferente → reescreve tudo

    # Após o header, varre as tabelas contando linhas
    header_len = len(header_str)
    rem = content[header_len:]
    lines = rem.splitlines()
    idx = 0
    expected = _expected_lines_per_table(t_steps, p_steps)
    for n in range(1, 10):
        if idx >= len(lines):
            # início desta tabela é o ponto de truncagem
            seek_str = f'$TABLE_{n} '
            table_pos = content.find(seek_str, header_len)
            if table_pos == -1:
                table_pos = len(content)  # só header
            return n, table_pos
        # precisa começar com a linha $TABLE_n
        if not lines[idx].startswith(f'$TABLE_{n} '):
            seek_str = f'$TABLE_{n} '
            table_pos = content.find(seek_str, header_len)
            if table_pos == -1:
                table_pos = len(content)
            return n, table_pos
        # tabela completa?
        if idx + expected > len(lines):
            seek_str = f'$TABLE_{n} '
            table_pos = content.find(seek_str, header_len)
            return n, table_pos
        idx += expected
    # Se chegou aqui, as 9 tabelas existem completas
    return 10, None

# === Parser do .rgp base para coletores ===
# === Parser do .rgp base para coletores (tolerante a marcadores) ===
def _parse_base_to_collectors(filename: str, t_steps: int, p_steps: int) -> Dict[str, Dict[int, List[float]]]:
    """
    Lê <filename>.rgp (base) e reconstrói:
      collectors = {"tables": {tbl: flat_vals}, "sat_t": {tbl: list}, "sat_phi": {tbl: list}}
    Sem usar CoolProp. Aceita arquivos COM ou SEM marcadores $T / $P / $SAT_T / $SAT_PHI.
    Ordem esperada por tabela: $TABLE_n, dims, T(tt), P(pp), GRID(tt*pp), Tsat(pp), Phi(pp).
    """
    collectors = {"tables": {}, "sat_t": {}, "sat_phi": {}}

    def _next_nonempty(it):
        for ln in it:
            s = ln.strip()
            if s:
                return s
        return None

    def _collect_floats_from_line(line: str) -> List[float]:
        toks = line.split()
        out = []
        for tok in toks:
            tok = tok.replace(',', '')
            try:
                out.append(float(tok))
            except Exception:
                # ignora lixo eventual
                pass
        return out

    def _collect_n_numbers(it, n: int, maybe_first_line: Optional[str] = None) -> List[float]:
        """Coleta exatamente n números de linhas subsequentes (sem exigir marcadores)."""
        nums: List[float] = []
        if maybe_first_line:
            nums.extend(_collect_floats_from_line(maybe_first_line))
        while len(nums) < n:
            line = _next_nonempty(it)
            if line is None:
                raise ValueError(f"Fim de arquivo ao tentar ler {n} floats.")
            # Se aparecer um marcador ($ALGUMA_COISA) no meio, isso é erro de formato
            if line.startswith('$'):
                raise ValueError(f"Esperado {n} floats, mas encontrei marcador: {line}")
            nums.extend(_collect_floats_from_line(line))
        # corta o excesso se houver
        return nums[:n]

    with open(filename, 'r', encoding='utf-8') as f:
        it = iter(f.readlines())

        while True:
            line = _next_nonempty(it)
            if line is None:
                break
            if not line.startswith('$TABLE_'):
                # pula coisas do cabeçalho ($$$$ etc.) até achar uma tabela
                continue

            # número da tabela
            try:
                tbl_num = int(line.split('_')[1])
            except Exception as exc:
                raise ValueError(f"Não consegui ler número de tabela em: {line}") from exc

            # dimensões (t_steps p_steps)
            dims_line = _next_nonempty(it)
            if dims_line is None:
                raise ValueError(f"TABELA_{tbl_num}: EOF ao ler dimensões.")
            dims_tokens = dims_line.split()
            if len(dims_tokens) < 2:
                raise ValueError(f"Dimensões inválidas para TABLE_{tbl_num}: '{dims_line}'")
            tt = int(dims_tokens[0]); pp = int(dims_tokens[1])

            # ----- eixo T (tt floats): pode haver um marcador '$T' OU vir direto números
            line_t = _next_nonempty(it)
            if line_t is None:
                raise ValueError(f"TABELA_{tbl_num}: EOF antes de T-axis.")
            if line_t.startswith('$T'):
                # ler T começando da próxima linha
                t_axis = _collect_n_numbers(it, tt, None)
            else:
                # a própria line_t já é a 1ª linha com números de T
                t_axis = _collect_n_numbers(it, tt, line_t)

            # ----- eixo P (pp floats): idem
            line_p = _next_nonempty(it)
            if line_p is None:
                raise ValueError(f"TABELA_{tbl_num}: EOF antes de P-axis.")
            if line_p.startswith('$P'):
                p_axis = _collect_n_numbers(it, pp, None)
            else:
                p_axis = _collect_n_numbers(it, pp, line_p)

            # ----- GRID (tt*pp floats): sem marcador no teu formato
            flat_vals = _collect_n_numbers(it, tt * pp, None)
            collectors["tables"][tbl_num] = flat_vals

            # ----- Tsat (pp floats): pode ou não ter $SAT_T
            line_sat_t = _next_nonempty(it)
            if line_sat_t is None:
                raise ValueError(f"TABELA_{tbl_num}: EOF antes de Tsat.")
            if line_sat_t.startswith('$SAT_T'):
                tsat_list = _collect_n_numbers(it, pp, None)
            else:
                tsat_list = _collect_n_numbers(it, pp, line_sat_t)
            collectors["sat_t"][tbl_num] = tsat_list

            # ----- φ_sat (pp floats): pode ou não ter $SAT_PHI
            line_sat_phi = _next_nonempty(it)
            if line_sat_phi is None:
                raise ValueError(f"TABELA_{tbl_num}: EOF antes de φ_sat.")
            if line_sat_phi.startswith('$SAT_PHI'):
                phi_list = _collect_n_numbers(it, pp, None)
            else:
                phi_list = _collect_n_numbers(it, pp, line_sat_phi)
            collectors["sat_phi"][tbl_num] = phi_list

    return collectors

# ======================================================
# MAIN
# ======================================================

def main():
    # Função para padronizar prompts
    def prompt(msg):
        return input(Fore.BLUE + msg + Style.RESET_ALL)

    user_ranges = None

    print(Fore.BLUE + "\nTipo de entrada:")
    print(Fore.BLUE + "1 - Fluido puro")
    print(Fore.BLUE + "2 - Mistura de fluidos\n" + Style.RESET_ALL)
    while True:
        tipo = prompt("Digite 1 para Fluido puro ou 2 para Mistura: ").strip()
        if tipo in ['1', '2']:
            break
        print(Fore.RED + "Opção inválida! Digite 1 ou 2." + Style.RESET_ALL)

    if tipo == '1':
        # FLUIDO PURO
        print(Fore.BLUE + "\nExemplo de fluido puro: " + Fore.WHITE + "CO2" + Style.RESET_ALL)
        while True:
            fluid = prompt("Digite o nome do fluido: ").strip().upper()
            try:
                t_crit, p_crit, t_triple, p_triple, r_especifico = get_mixture_properties(fluid)
                break
            except Exception as exc:
                print(Fore.RED + f"Nome inválido ({exc})! Tente novamente." + Style.RESET_ALL)
        mixture = fluid
    else:
        # MISTURA
        print(Fore.BLUE + "\nExemplo de mistura de fluidos: " +
              Fore.WHITE + "Methane[0.4]&CO2[0.6]\n" + Style.RESET_ALL)
        while True:
            fluid1 = prompt("Digite o nome do primeiro fluido da mistura: ").strip().upper()
            try:
                _ = get_mixture_properties(fluid1)
                break
            except Exception as exc:
                print(Fore.RED + f"Nome inválido ({exc})! Tente novamente." + Style.RESET_ALL)
        while True:
            fluid2 = prompt("Digite o nome do segundo fluido da mistura: ").strip().upper()
            try:
                _ = get_mixture_properties(fluid2)
                break
            except Exception as exc:
                print(Fore.RED + f"Nome inválido ({exc})! Tente novamente." + Style.RESET_ALL)
        # Frações molares
        while True:
            while True:
                try:
                    frac1 = float(prompt(f"Digite a fração molar do {fluid1} (ex: 0.4): ").strip())
                    frac2 = float(prompt(f"Digite a fração molar do {fluid2} (ex: 0.6): ").strip())
                except ValueError:
                    print(Fore.RED + "Fração inválida! Digite um número, como 0.4 ou 0.6." + Style.RESET_ALL)
                    continue
                soma = frac1 + frac2
                if abs(soma - 1.0) > 1e-4:
                    print(Fore.RED + "A soma das frações deve ser igual a 1.0! Tente novamente." + Style.RESET_ALL)
                    continue
                break
            mixture = f"{fluid1}[{frac1}]&{fluid2}[{frac2}]"
            try:
                t_crit, p_crit, t_triple, p_triple, r_especifico = get_mixture_properties(mixture)
                break
            except Exception:
                mixture_inv = f"{fluid2}[{frac2}]&{fluid1}[{frac1}]"
                try:
                    t_crit, p_crit, t_triple, p_triple, r_especifico = get_mixture_properties(mixture_inv)
                    print(Fore.YELLOW + f"Ordem corrigida para: {mixture_inv}" + Style.RESET_ALL)
                    mixture = mixture_inv
                    break
                except Exception as exc:
                    print(Fore.RED + f"Sintaxe da mistura inválida: {exc}" + Style.RESET_ALL)
                    print(Fore.MAGENTA + "Tente novamente com frações diferentes." + Style.RESET_ALL)

    # === Mostrar ranges padrão ===
    print(Fore.GREEN + "\nPropriedades críticas do fluido " + Fore.WHITE + f"{mixture}:" + Style.RESET_ALL)
    print(Fore.CYAN + f"  T_triple = {t_triple:.2f} K,  T_crit = {t_crit:.2f} K" + Style.RESET_ALL)
    print(Fore.CYAN + f"  P_triple = {p_triple/1e6:.3f} MPa,  P_crit = {p_crit/1e6:.3f} MPa" + Style.RESET_ALL)

    ranges_padrao = {
        "tmin_superheat": round(t_triple + 2.0, 7),
        "tmax_superheat": round(t_crit + 100.0, 7),
        "pmin_superheat": max(p_triple + 1e4, p_crit * 0.01),
        "pmax_superheat": round(p_crit * 1.2, 7),
        "tmin_saturation": round(t_triple + 1.0, 7),
        "tmax_saturation": round(t_crit - 1.0, 7),
    }
    print(Fore.YELLOW + "\nRanges automáticos que serão usados:" + Style.RESET_ALL)
    print(Fore.WHITE + f"  Temperatura: {ranges_padrao['tmin_superheat']:.1f} K até {ranges_padrao['tmax_superheat']:.1f} K" + Style.RESET_ALL)
    print(Fore.WHITE + f"  Pressão:     {ranges_padrao['pmin_superheat']/1e6:.3f} MPa até {ranges_padrao['pmax_superheat']/1e6:.3f} MPa\n" + Style.RESET_ALL)

    # Range manual?
    valid_sim = ["s", "sim", "y", "yes"]
    valid_nao = ["n", "nao", "não", "no"]
    while True:
        range_manual = prompt("Deseja inserir manualmente os limites de Temperatura e Pressão? (S/N): ").strip().lower()
        if range_manual in valid_sim:
            range_manual = True
            break
        elif range_manual in valid_nao:
            range_manual = False
            break
        else:
            print(Fore.RED + "Opção inválida! Digite S/N, Sim/Não, Y/Yes/No." + Style.RESET_ALL)

    if range_manual:
        print(Fore.CYAN + "\nDefinir o intervalo de Temperatura e Pressão:" + Style.RESET_ALL)
        while True:
            try:
                tmin = float(prompt("Temperatura MÍNIMA [K]: ").strip())
                tmax = float(prompt("Temperatura MÁXIMA [K]: ").strip())
                pmin = float(prompt("Pressão MÍNIMA [Pa]: ").strip())
                pmax = float(prompt("Pressão MÁXIMA [Pa]: ").strip())
                if tmin >= tmax or pmin >= pmax:
                    print(Fore.RED + "Os valores mínimos devem ser menores que os máximos!" + Style.RESET_ALL)
                    continue
                break
            except ValueError:
                print(Fore.RED + "Digite valores numéricos válidos!" + Style.RESET_ALL)
        user_ranges = dict(
            tmin_superheat=tmin,
            tmax_superheat=tmax,
            pmin_superheat=pmin,
            pmax_superheat=pmax,
            tmin_saturation=tmin,
            tmax_saturation=tmax,
        )
    else:
        user_ranges = None

    # Número de pontos
    while True:
        try:
            t_steps = int(prompt("Digite o número de pontos de Temp. (ex: 15): ").strip())
            p_steps = int(prompt("Digite o número de pontos de Pres. (ex: 10): ").strip())
            if t_steps % 5 != 0 or p_steps % 5 != 0:
                print(Fore.RED + "ERRO: Número de pontos deve ser múltiplo de 5 para formato RGP!" + Style.RESET_ALL)
                continue
            break
        except ValueError:
            print(Fore.RED + "Digite um número inteiro!" + Style.RESET_ALL)

    # Nome da tabela RGP (máx 8 chars)
    while True:
        table_name = prompt("Nome do fluido/mistura na tabela RGP(máx 8 caracteres, ex: CO2/METHCO2): ").strip().upper()
        if len(table_name) <= 8:
            break
        print(Fore.RED + f"Nome muito longo ({len(table_name)} chars)! Máximo 8 caracteres." + Style.RESET_ALL)

    # Nome do arquivo
    file_name = prompt("Digite o nome do arquivo de saída (apenas o nome, sem extensão): ").strip()
    if not file_name:
        file_name = f'RGPTable_{mixture.replace("&", "_")}_v4'
    if not file_name.lower().endswith('.rgp'):
        file_name += '.rgp'

    start_time = time.time()
    warning_counter = WarningCounter()
    logging.basicConfig(level=logging.WARNING, handlers=[warning_counter])

    if user_ranges:
        vals = {
            "pmin_superheat": user_ranges["pmin_superheat"],
            "pmax_superheat": user_ranges["pmax_superheat"],
            "tmin_saturation": user_ranges["tmin_saturation"],
            "tmax_saturation": round(t_crit, 7),
            "tmin_superheat": user_ranges["tmin_superheat"],
            "tmax_superheat": user_ranges["tmax_superheat"],
            "supercooling": 0.0,
            "p_critical": round(p_crit, 7),
            "p_triple": round(p_triple, 7),
            "t_critical": round(t_crit, 7),
            "t_triple": round(t_triple, 7),
            "gas_constant": round(r_especifico, 7),
        }
    else:
        vals = {
            "pmin_superheat": round(p_triple + 1e5, 7),
            "pmax_superheat": round(p_crit - 1e5, 7),
            "tmin_saturation": round(t_triple + 1.0, 7),
            "tmax_saturation": round(t_crit - 1.0, 7),
            "tmin_superheat": round(t_triple + 2.0, 7),
            "tmax_superheat": round(t_crit - 1.0 + 50.0, 7),
            "supercooling": 0.0,
            "p_critical": round(p_crit, 7),
            "p_triple": round(p_triple, 7),
            "t_critical": round(t_crit, 7),
            "t_triple": round(t_triple, 7),
            "gas_constant": round(r_especifico, 7),
        }

    # Nomes amigáveis e unidades
    nomes_bonitos = {
        "pmin_superheat": "Pres. mín. superaquec.",
        "pmax_superheat": "Pres. máx. superaquec.",
        "tmin_saturation": "Temp. mín. saturação",
        "tmax_saturation": "Temp. máx. saturação",
        "tmin_superheat": "Temp. mín. superaquec.",
        "tmax_superheat": "Temp. máx. superaquec.",
        "supercooling": "Supercooling",
        "p_critical": "Pres. crítica",
        "p_triple": "Pres. no ponto tríplice",
        "t_critical": "Temp. crítica",
        "t_triple": "Temp. no ponto tríplice",
        "gas_constant": "Constante dos gases (R)"
    }
    unidades = {
        "pmin_superheat": "MPa",
        "pmax_superheat": "MPa",
        "p_critical": "MPa",
        "p_triple": "MPa",
        "tmin_saturation": "K",
        "tmax_saturation": "K",
        "tmin_superheat": "K",
        "tmax_superheat": "K",
        "t_critical": "K",
        "t_triple": "K",
        "gas_constant": "J/(kg·K)",
        "supercooling": "K"
    }
    pressao_mpa = {"pmin_superheat", "pmax_superheat", "p_critical", "p_triple"}

    print(Fore.GREEN + "Propriedades do fluido " + Fore.WHITE + f"{mixture}:\n" + Style.RESET_ALL)
    for key, val in vals.items():
        unidade = unidades.get(key, "")
        nome_bonito = nomes_bonitos.get(key, key)
        if key in pressao_mpa:
            val = val / 1_000_000
            print(
                Fore.YELLOW + f"  {nome_bonito:<32} " +
                Fore.YELLOW + "=" +
                Fore.WHITE + f" {val:8.3f} " +
                Fore.MAGENTA + f"{unidade}" +
                Style.RESET_ALL
            )
        elif unidade == "K":
            print(
                Fore.YELLOW + f"  {nome_bonito:<32} " +
                Fore.YELLOW + "=" +
                Fore.WHITE + f" {val:8.1f} " +
                Fore.MAGENTA + f"{unidade}" +
                Style.RESET_ALL
            )
        elif unidade == "J/(kg·K)":
            print(
                Fore.YELLOW + f"  {nome_bonito:<32} " +
                Fore.YELLOW + "=" +
                Fore.WHITE + f" {val:8.1f} " +
                Fore.MAGENTA + f"{unidade}" +
                Style.RESET_ALL
            )
        else:
            print(
                Fore.YELLOW + f"  {nome_bonito:<32} " +
                Fore.YELLOW + "=" +
                Fore.WHITE + f" {val} " +
                Fore.MAGENTA + f"{unidade}" +
                Style.RESET_ALL
            )

    # Malha
    t_resol = (vals["tmax_superheat"] - vals["tmin_superheat"]) / (t_steps - 1)
    p_resol = (vals["pmax_superheat"] - vals["pmin_superheat"]) / (p_steps - 1)
    t_vals = [vals["tmin_superheat"] + i * t_resol for i in range(t_steps)]
    p_vals = [vals["pmin_superheat"] + j * p_resol for j in range(p_steps)]

    # Coletores p/ pós-processamento
    collectors = {"tables": {}, "sat_t": {}, "sat_phi": {}}

    print("\n")
    # --- Suporte a retomada ---
    header_str = _header_string(table_name, t_steps, p_steps, vals)
    start_tbl, truncate_pos = _detect_resume(file_name, header_str, t_steps, p_steps)

    # Prepara mapa de propriedades
    property_map = generate_property_map(mixture, vals["t_critical"], vals["p_critical"], t_resol)

    if start_tbl == 1:
        # Sem arquivo anterior compatível → escrever do zero
        with open(file_name, 'w', encoding='utf-8') as file_obj:
            file_obj.write(header_str)
            for tbl_num in tqdm(range(1, 10), desc=Fore.CYAN + "Tabelas principais" + Style.RESET_ALL):
                get_prop = property_map[tbl_num][1]
                write_table_block(
                    file_obj, tbl_num, t_steps, p_steps, t_vals, p_vals, get_prop,
                    collector=collectors
                )
                write_saturation_block(
                    file_obj, tbl_num, p_vals, mixture, vals["t_triple"], vals["t_critical"],
                    collector=collectors
                )
    elif start_tbl <= 9:
        # Retomar: truncar no início da primeira tabela incompleta e continuar
        with open(file_name, 'r+', encoding='utf-8') as file_obj:
            if truncate_pos is not None:
                file_obj.truncate(truncate_pos)
                file_obj.seek(truncate_pos)
            for tbl_num in tqdm(range(start_tbl, 10), desc=Fore.CYAN + f"Retomando a partir da Tabela {start_tbl}" + Style.RESET_ALL):
                get_prop = property_map[tbl_num][1]
                write_table_block(
                    file_obj, tbl_num, t_steps, p_steps, t_vals, p_vals, get_prop,
                    collector=collectors
                )
                write_saturation_block(
                    file_obj, tbl_num, p_vals, mixture, vals["t_triple"], vals["t_critical"],
                    collector=collectors
                )
    else:
        # Arquivo já completo (9 tabelas). Não reescreve o original.
        print(Fore.YELLOW + "Arquivo original já completo. Pulando reescrita do .rgp base." + Style.RESET_ALL)

    # Passo 2/2: reescrita (sem recálculo) — início claro + barra de progresso
    print(Fore.MAGENTA + "\n==== PASSO 2/2: REESCRITA (sem recálculo) ====" + Style.RESET_ALL)

    # Garante que os coletores estão preenchidos mesmo se o arquivo base já estava pronto
    if not collectors["tables"]:
        collectors = _parse_base_to_collectors(file_name, t_steps, p_steps)

    corrected_name = file_name[:-4] + "_CORRIGIDO.rgp" if file_name.lower().endswith(".rgp") else file_name + "_CORRIGIDO.rgp"

    # Corrigir grid principal com barra de progresso
    corrected_tables = {}
    for tbl_num in tqdm(range(1, 10), desc="Reescrita: grid (1/2)", unit="tbl"):
        flat = collectors["tables"][tbl_num]
        corrected_tables[tbl_num] = fix_grid(flat, n_t=t_steps, zero_tol=1e-9, min_const_run=2, eps_equal=1e-12)

    # Corrigir saturações com barra de progresso (Tsat e φ), SEM PropsSI
    corrected_sat_t = {}
    corrected_sat_phi = {}
    for tbl_num in tqdm(range(1, 10), desc="Reescrita: saturações (2/2)", unit="tbl"):
        ts = collectors["sat_t"][tbl_num]
        ph = collectors["sat_phi"][tbl_num]
        corrected_sat_t[tbl_num]  = fix_vector(ts,  zero_tol=1e-9, min_const_run=2, eps_equal=1e-12)
        corrected_sat_phi[tbl_num] = fix_vector(ph, zero_tol=1e-9, min_const_run=2, eps_equal=1e-12)

    # Reescrever arquivo completo com os valores corrigidos
    with open(corrected_name, 'w', encoding='utf-8') as file_obj:
        write_header(file_obj, table_name, t_steps, p_steps, vals)
        for tbl_num in range(1, 10):
            write_table_block_from_values(
                file_obj, tbl_num, t_steps, p_steps, t_vals, p_vals, corrected_tables[tbl_num]
            )
            write_saturation_block_from_values(
                file_obj, corrected_sat_t[tbl_num], corrected_sat_phi[tbl_num]
            )

    elapsed = time.time() - start_time
    print(
        Fore.GREEN + "\nArquivos gerados com sucesso:" + Style.RESET_ALL
    )
    print(Fore.WHITE + f"  • {file_name}" + Style.RESET_ALL)
    print(Fore.WHITE + f"  • {corrected_name}  (pós-processado: buracos corrigidos)" + Style.RESET_ALL)
    print(Fore.BLUE + "Tempo total: " +
          Fore.WHITE + f"{elapsed:.1f} segundos" + Style.RESET_ALL)

    # Relatório de warnings
    if isinstance(logging.getLogger().handlers[0], WarningCounter):
        warning_counter = logging.getLogger().handlers[0]
        if warning_counter.warning_count == 0:
            print(
                Fore.GREEN + "Todas as propriedades foram calculadas com sucesso, sem erros numéricos." +
                Style.RESET_ALL
            )
        else:
            print(
                Fore.YELLOW +
                f"Foram encontrados {warning_counter.warning_count} valores não calculados (NaN):" +
                Style.RESET_ALL
            )
            for tabela, count in warning_counter.warnings_por_tabela.items():
                print(Fore.YELLOW + f"  Tabela {tabela}: {count} warnings" + Style.RESET_ALL)
            print(Fore.YELLOW + "Verifique os ranges e condições da mistura." + Style.RESET_ALL)

    print(
        Fore.CYAN +
        "\nVocê pode importar o arquivo diretamente no ANSYS CFX. O sufixo _CORRIGIDO traz a versão com buracos preenchidos por tendência.\n" +
        Style.RESET_ALL
    )

if __name__ == "__main__":
    main()

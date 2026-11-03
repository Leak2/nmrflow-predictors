"""文献风格报道值格式化。

排版约定（半角）：
  - 耦合写 ``J= 7.0 Hz``（半角等号、等号后一个空格）
  - 区间用 ASCII 连字符 ``-``
  - 峰之间用 ``, `` 分隔，行尾以 ``.`` 结束

> 曾经有一档 ``cn-pdf`` 风格，刻意复刻「中文 PDF 复制稿」的排版噪声
> （全角 ＝、整数与小数点之间插空格、U+2011 不换行连字符）。它是为了一次性
> 对照交付做的，实际使用频率极低，却会让默认输出出现 ``3 .64``、``J＝`` 这种
> 东西，容易被误当成产品行为。0.2.6 起整档删除。
"""

from __future__ import annotations

from ..core.types import NmrBlock, Peak

JOIN = ", "


def _num(value: float, digits: int) -> str:
    return f"{value:.{digits}f}"


_J_SYMBOL = {"19F": "JC-F", "31P": "JC-P", "1H": "J", "13C": "JC-C"}

# 积分后缀：2H / 3F / 1P。13C 不做积分，故不在表内。
_NUCLEUS_SYMBOL = {"1H": "H", "19F": "F", "31P": "P"}


def _j_symbol(nucleus: str) -> str:
    """13C 谱里异核耦合的记号：JC-F / JC-P / J（同核）。"""
    return _J_SYMBOL.get(nucleus, "J")


def _integration_text(p: Peak, nucleus: str) -> str:
    """`2F` / `1P`：峰对应的等价核个数，按核给出元素符号。"""
    symbol = _NUCLEUS_SYMBOL.get(nucleus)
    if not symbol or not p.integration:
        return ""
    return f"{int(round(p.integration))}{symbol}"


def _range(p: Peak, digits: int) -> str:
    if abs(p.shift_hi - p.shift_lo) < 10 ** -digits / 2:
        return _num(p.shift, digits)
    return f"{_num(p.shift_lo, digits)}-{_num(p.shift_hi, digits)}"


def format_peak(p: Peak, nucleus: str, *, digits: int = 2) -> str:
    """生成单个峰的描述，例如 4.36 (q, J = 7.2 Hz, 2H, -CH2)。

    只有 1H 写 δ 区间（`7.80-7.86 (m, 4H)`），因为同一个多重峰的
    各条谱线本来就散在一段区间里；13C / 19F / 31P 的文献写法是
    **一个中心值 + 多重度 + J**，写区间反而读不出来。
    """
    span = _range(p, digits) if nucleus == "1H" else _num(p.shift, digits)

    if nucleus == "13C":
        # 13C 通常只写位移，只有异核耦合才带 J（JC-F / JC-P ...）
        if p.j_values:
            j_txt = ", ".join(
                f" {_j_symbol(nuc)}= {j:.1f} Hz"
                for j, nuc in zip(p.j_values, p.j_partners)
            )
            return f"{span} ({p.multiplicity}{j_txt}, {p.label})".replace("  ", " ")
        return f"{span}"

    if nucleus in ("19F", "31P"):
        # 与 1H 一致，给出积分核数（2F / 1P ...），文献写法也是如此
        n_txt = _integration_text(p, nucleus)
        if p.j_values:
            j_txt = ", ".join(f"J= {j:.1f} Hz" for j in p.j_values)
            body = f"{p.multiplicity}, {j_txt}"
        else:
            body = p.multiplicity
        if n_txt:
            body += f", {n_txt}"
        return f"{span} ({body})"

    # 1H
    n_h = int(round(p.integration))
    n_h_txt = f"{n_h}H"
    if p.multiplicity == "s" and not p.j_values:
        if p.label:
            return f"{span} ({p.multiplicity}, {n_h_txt}, {p.label})"
        return f"{span} ({p.multiplicity}, {n_h_txt})"
    if p.j_values:
        js = ", ".join(f"J= {j:.1f} Hz" for j in p.j_values)
        body = f"{p.multiplicity}, {js}, {n_h_txt}"
    else:
        body = f"{p.multiplicity}, {n_h_txt}"
    if p.label:
        body += f", {p.label}"
    return f"{span} ({body})"


def format_block(block: NmrBlock, *, digits: int = 2,
                 nucleus_digits: dict[str, int] | None = None) -> str:
    """生成一整个核的报道行。"""
    nd = nucleus_digits or {}
    d = nd.get(block.nucleus, digits)
    parts = [format_peak(p, block.nucleus, digits=d) for p in block.peaks]
    freq = f"{block.freq_mhz:.0f} MHz"
    head = f"{block.nucleus} NMR ({block.solvent}, {freq}) \u03b4: "
    return head + JOIN.join(parts) + "."

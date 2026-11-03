"""从 NMRShiftDB2 的 SD 文件构建 HOSE 位移索引。

数据来源
--------
``nmrshiftdb2withsignals.sd``（151 MB，SourceForge `nmrshiftdb2/data/`，下载无需认证）。
数据许可为 *nmrshiftdb2 Database License*（ODbL 衍生版），其 4.5 条要求
「依赖该库实现功能的软件必须使用 OSI 认可的开源许可」——本项目据此开源。
分发索引时须按 4.2/4.3 附许可证与声明。

做法
----
1. 只取**实验谱**：带 ``NMRMethod`` / ``NMRProgram`` / ``NMRBasisSet`` / ``GeomMethod``
   任一属性的记录是 CDK 计算出来的，若混入索引等于拿它自己的 HOSE 预测结果训练 HOSE。
2. 谱字段格式 ``位移;强度;原子序号|...``，原子序号是 molfile 的 0-based 索引。
   **¹H 的序号指向带氢的重原子**（不是氢本身），因此四个核统一按
   「该原子自身的 HOSE 码」建键，不需要另造氢节点。
3. 每个原子只存最深层（默认 4 层）的码。查询时把码截断到第 k 层再做**前缀范围扫描**
   （``code >= prefix AND code < prefix + U+FFFF``），配合 ``level >= k`` 限制，
   就能复用同一批行做逐层降级，不必为每层各存一份。
4. 每行存 (n, Σshift, Σshift²)，合并时直接相加即得合并后的均值与标准差。

用法
----
    python tools/build_hose_index.py --sd /path/to/nmrshiftdb2withsignals.sd
"""

from __future__ import annotations

import argparse
import gzip
import os
import sqlite3
import sys
import time
from pathlib import Path
from typing import Callable

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from nmrflow.predict.hose import hose_codes  # noqa: E402

# 带这些属性的记录是计算谱，不是实验谱
CALC_PROPS = ("NMRMethod", "NMRProgram", "NMRBasisSet", "GeomMethod", "NMRLocalis")

SUPPORTED = ("1H", "13C", "19F", "31P")


def _iter_shift_entries(body: str):
    """解析 ``位移;强度;原子序号|...``，产出 (位移 ppm, 原子索引)。"""
    for chunk in body.split("|"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(";")
        if len(parts) < 3:
            continue
        try:
            shift = float(parts[0])
            idx = int(float(parts[2]))
        except ValueError:
            continue
        yield shift, idx


def _atom_matches_nucleus(atom: Chem.Atom, nucleus: str) -> bool:
    sym = atom.GetSymbol()
    if nucleus == "1H":
        # ¹H 的序号指向带氢的重原子
        return atom.GetTotalNumHs() > 0
    return sym == nucleus[2:]


def build(
    sd_path: Path,
    out_path: Path,
    spheres: int,
    nuclei: tuple[str, ...],
    limit: int | None,
    measure: bool,
    keep: "Callable[[Chem.Mol], bool] | None" = None,
    verbose: bool = True,
) -> dict:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()

    con = sqlite3.connect(out_path)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA temp_store=FILE")
    con.execute("CREATE TABLE raw(code TEXT NOT NULL, nucleus TEXT NOT NULL, level INTEGER NOT NULL, shift REAL NOT NULL)")

    supp = Chem.SDMolSupplier(str(sd_path), removeHs=False, sanitize=False)
    stats = {
        "mol": 0,
        "skip_none": 0,
        "skip_calc": 0,
        "skip_filter": 0,
        "skip_bad": 0,
        "skip_empty": 0,
    }
    per_nucleus: dict[str, int] = {n: 0 for n in nuclei}
    batch: list[tuple[str, str, int, float]] = []
    t0 = time.perf_counter()

    for mol in supp:
        stats["mol"] += 1
        if limit and stats["mol"] > limit:
            break
        if mol is None:
            stats["skip_none"] += 1
            continue
        if any(mol.HasProp(p) for p in CALC_PROPS):
            stats["skip_calc"] += 1
            continue
        if keep is not None and not keep(mol):
            stats["skip_filter"] += 1
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            stats["skip_bad"] += 1
            continue

        codes = hose_codes(mol, spheres)
        if not codes:
            stats["skip_empty"] += 1
            continue
        n_atoms = mol.GetNumAtoms()

        for nucleus in nuclei:
            prefix = f"Spectrum {nucleus} "
            for prop in mol.GetPropNames():
                if not prop.startswith(prefix):
                    continue
                for shift, idx in _iter_shift_entries(mol.GetProp(prop)):
                    if not (0 <= idx < n_atoms):
                        continue
                    if not _atom_matches_nucleus(mol.GetAtomWithIdx(idx), nucleus):
                        continue
                    code = codes.get(idx)
                    if not code:
                        continue
                    batch.append((code, nucleus, code.count(";"), shift))
                    per_nucleus[nucleus] += 1

        if len(batch) >= 50000:
            con.executemany("INSERT INTO raw VALUES (?,?,?,?)", batch)
            batch.clear()
        if verbose and stats["mol"] % 5000 == 0:
            print(
                f"  {stats['mol']:>6} 分子  行数 {sum(per_nucleus.values()):>9,}  "
                f"{time.perf_counter() - t0:.0f}s",
                flush=True,
            )

    if batch:
        con.executemany("INSERT INTO raw VALUES (?,?,?,?)", batch)
    con.commit()

    print("\n聚合中 ...", flush=True)
    con.execute(
        """
        CREATE TABLE hose(
            nucleus  TEXT    NOT NULL,
            code     TEXT    NOT NULL,
            level    INTEGER NOT NULL,
            n        INTEGER NOT NULL,
            sum_shift REAL   NOT NULL,
            sum_sq    REAL   NOT NULL,
            PRIMARY KEY (nucleus, code)
        ) WITHOUT ROWID
        """
    )
    con.execute(
        """
        INSERT INTO hose
        SELECT nucleus, code, MAX(level), COUNT(*), SUM(shift), SUM(shift*shift)
        FROM raw GROUP BY nucleus, code
        """
    )
    con.execute("DROP TABLE raw")
    con.commit()
    con.execute("VACUUM")
    con.commit()

    rows = con.execute("SELECT nucleus, COUNT(*), SUM(n) FROM hose GROUP BY nucleus ORDER BY nucleus").fetchall()
    levels = con.execute("SELECT level, COUNT(*) FROM hose GROUP BY level ORDER BY level").fetchall()
    con.close()

    db_bytes = out_path.stat().st_size
    gz_bytes = None
    if measure:
        print("压缩测量中 ...", flush=True)
        tmp = out_path.with_suffix(out_path.suffix + ".gz")
        with open(out_path, "rb") as src, gzip.open(tmp, "wb", compresslevel=9) as dst:
            while True:
                chunk = src.read(1 << 22)
                if not chunk:
                    break
                dst.write(chunk)
        gz_bytes = tmp.stat().st_size
        tmp.unlink()

    return {
        "stats": stats,
        "per_nucleus": per_nucleus,
        "rows": rows,
        "levels": levels,
        "db_bytes": db_bytes,
        "gz_bytes": gz_bytes,
        "seconds": time.perf_counter() - t0,
        "path": out_path,
    }


def _mb(n: int | None) -> str:
    return "—" if n is None else f"{n / 1024 / 1024:.1f} MB"


def main() -> int:
    ap = argparse.ArgumentParser(description="构建 HOSE 位移索引")
    ap.add_argument("--sd", required=True, type=Path, help="nmrshiftdb2withsignals.sd 路径")
    ap.add_argument("--out", type=Path, default=Path("build/hose_index.sqlite"))
    ap.add_argument("--spheres", type=int, default=4)
    ap.add_argument("--nuclei", default=",".join(SUPPORTED))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-measure", action="store_true", help="跳过 gzip 压缩测量")
    args = ap.parse_args()

    nuclei = tuple(n.strip() for n in args.nuclei.split(",") if n.strip())
    bad = [n for n in nuclei if n not in SUPPORTED]
    if bad:
        print(f"不支持的核: {bad}", file=sys.stderr)
        return 2

    print(f"输入: {args.sd}")
    print(f"球层: {args.spheres}   核: {', '.join(nuclei)}")
    r = build(args.sd, args.out, args.spheres, nuclei, args.limit, not args.no_measure)

    s = r["stats"]
    print("\n" + "=" * 58)
    print(f"分子记录      {s['mol']:,}")
    print(f"  解析失败    {s['skip_none']:,}")
    print(f"  计算谱跳过  {s['skip_calc']:,}")
    print(f"  sanitize 失败 {s['skip_bad']:,}")
    print(f"索引行数      {sum(x[1] for x in r['rows']):,}")
    print("\n按核:")
    for nuc, n_codes, n_shifts in r["rows"]:
        avg = n_shifts / n_codes if n_codes else 0
        print(f"  {nuc:4s} 码 {n_codes:>8,}   位移条目 {n_shifts:>9,}   平均每码 {avg:5.2f} 个")
    print("\n按球层:")
    for lv, cnt in r["levels"]:
        print(f"  球层 {lv}: {cnt:,} 个码")
    print("\n体积:")
    print(f"  SQLite      {_mb(r['db_bytes'])}")
    print(f"  gzip -9     {_mb(r['gz_bytes'])}")
    print(f"\n输出: {r['path']}")
    print(f"耗时: {r['seconds']:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

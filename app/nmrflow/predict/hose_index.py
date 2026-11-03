"""HOSE 位移索引的查询层（Tier B）。

索引由 ``tools/build_hose_index.py`` 从 NMRShiftDB2 实验谱构建，
结构为 ``hose(nucleus, code, level, n, sum_shift, sum_sq)``，主键 ``(nucleus, code)``。

逐层降级
--------
每个原子只存了最深一层（默认 4 层）的码。因为 HOSE 码有**前缀性质**
（``truncate(code, k)`` 等于直接按 k 层生成），查询第 k 层时把码截断到第 k 层，
再做前缀范围扫描即可复用同一批行：

    code >= 'prefix' AND code < 'prefix\\uffff' AND level >= k

先试最深一层；命中邻居数不足 ``min_neighbors`` 就退一层，直到 ``min_level``。
这与 CDK 的 HOSE 预测器行为一致，也是"命中不足就有记录地降级"的来源：
球层与邻居数可以直接写进 ``Provenance.note``，让报告说明这条预测是
"由 14 条同环境数据库条目取均值"还是"仅 1 条"。

统计量
------
每行存 (n, Σx, Σx²)，前缀扫描把多行相加即得合并后的精确均值与方差，
不需要把原始位移一条条存下来。
"""

from __future__ import annotations

import math
import sqlite3
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

from nmrflow.predict.hose import hose_codes, truncate

DEFAULT_INDEX_PATH = Path(__file__).resolve().parent.parent / "data" / "hose_index.sqlite"

# 开发期索引留在 build/ 下（32 MB，不进源码目录），发布时再拷进 nmrflow/data/
_BUILD_INDEX_PATH = Path(__file__).resolve().parents[2] / "build" / "hose_index.sqlite"

# 前缀范围扫描的上界：U+FFFF 大于任何码字符
_RANGE_HI = "\uffff"


def resolve_index_path(explicit: str | Path | None = None) -> Path | None:
    """按「显式指定 -> 包内 -> 开发期 build/」的顺序找索引；都没有则返回 None。"""
    for candidate in (explicit, DEFAULT_INDEX_PATH, _BUILD_INDEX_PATH):
        if candidate is None:
            continue
        path = Path(candidate)
        if path.exists():
            return path
    return None


@dataclass(frozen=True)
class HoseHit:
    """一次查询的结果。``level`` 是实际命中的球层数，``n`` 是同码邻居数。"""

    level: int
    n: int
    mean: float
    std: float

    @property
    def note(self) -> str:
        return f"HOSE {self.level} 层球, {self.n} 条同环境记录, σ={self.std:.2f}"


class HoseIndex:
    """只读的 HOSE 索引。"""

    def __init__(self, path: str | Path = DEFAULT_INDEX_PATH):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(
                f"未找到 HOSE 索引 {self.path}；先运行 tools/build_hose_index.py"
            )
        self._con = sqlite3.connect(f"file:{self.path.as_posix()}?mode=ro", uri=True)
        # 键必须含全部查询参数：同一码在不同 min_neighbors 下结果不同
        self._cache: dict[tuple[str, str, int, int, int], HoseHit | None] = {}

    def close(self) -> None:
        self._con.close()

    def __enter__(self) -> HoseIndex:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def lookup(
        self,
        code: str,
        nucleus: str,
        *,
        max_level: int = 4,
        min_level: int = 1,
        min_neighbors: int = 1,
    ) -> HoseHit | None:
        """从 ``max_level`` 起逐层降级，返回第一个邻居数达标的命中。

        ``min_neighbors`` 设成 3 可要求"至少 3 条同环境记录"，代价是球层更浅。
        若降到底仍不达标，返回**支持数最多、其次球层最深**的那次命中——
        球层更深意味着环境更具体，在支持数相同时优先保留。
        """
        key = (nucleus, code, max_level, min_level, min_neighbors)
        cached = self._cache.get(key)
        if cached is not None or key in self._cache:
            return cached

        best: HoseHit | None = None
        for level in range(max_level, min_level - 1, -1):
            prefix = truncate(code, level)
            hit = self._query_prefix(prefix, level, nucleus)
            if hit is None:
                continue
            if best is None or hit.n > best.n:
                best = hit
            if hit.n >= min_neighbors:
                self._cache[key] = hit
                return hit

        self._cache[key] = best
        return best

    def _query_prefix(self, prefix: str, level: int, nucleus: str) -> HoseHit | None:
        row = self._con.execute(
            """
            SELECT SUM(n), SUM(sum_shift), SUM(sum_sq)
            FROM hose
            WHERE nucleus = ? AND level >= ? AND code >= ? AND code < ?
            """,
            (nucleus, level, prefix, prefix + _RANGE_HI),
        ).fetchone()
        if not row or not row[0]:
            return None
        n, total, total_sq = int(row[0]), float(row[1]), float(row[2])
        mean = total / n
        var = total_sq / n - mean * mean
        return HoseHit(level=level, n=n, mean=mean, std=math.sqrt(max(0.0, var)))

    def predict(
        self,
        mol: Chem.Mol,
        nucleus: str = "1H",
        *,
        max_spheres: int = 4,
        min_neighbors: int = 1,
    ) -> dict[int, HoseHit]:
        """对分子中所有可预测原子给出位移。返回 ``{原子索引: HoseHit}``。

        ¹H 的键是**带氢的重原子**（与 NMRShiftDB2 的归属约定一致），
        因此 CH3 上的三个氢共用一条预测值。
        """
        codes = hose_codes(mol, max_spheres)
        out: dict[int, HoseHit] = {}
        for idx, code in codes.items():
            atom = mol.GetAtomWithIdx(idx)
            if nucleus == "1H":
                if atom.GetTotalNumHs() == 0:
                    continue
            elif atom.GetSymbol() != nucleus[2:]:
                continue
            hit = self.lookup(
                code, nucleus, max_level=max_spheres, min_neighbors=min_neighbors
            )
            if hit is not None:
                out[idx] = hit
        return out

    def stats(self) -> list[tuple[str, int, int]]:
        """``[(核, 码数, 位移条目数)]``。"""
        return self._con.execute(
            "SELECT nucleus, COUNT(*), SUM(n) FROM hose GROUP BY nucleus ORDER BY nucleus"
        ).fetchall()

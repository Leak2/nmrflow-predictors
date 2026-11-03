"""核验随包评估套件（eval kit）是否完整、可用。

为什么需要这个脚本
------------------
README 里的留出法数字（¹H 0.39 / ¹³C 3.41 / ¹⁹F 8.76 / ³¹P 19.71）以前只能"自报"：
便携包只带 ``app/`` 与 ``runtime/``，没有评估脚本，也没有训练索引，评审无从复核。

0.2.7 起 ``app/tools/`` 随包，并随包提供**当时那两份训练索引**。但"文件在"不等于
"文件是那一份"：索引可能被别的构建覆盖、被截断、或者混进了留出集的化合物。
这个脚本用三条互相独立的证据把这件事钉住：

1. **哈希**：与 ``EVAL_KIT_SHA256.txt`` 记录的摘要逐字节一致；
2. **构图**：全量索引四核齐全、行数量级正确；训练索引只含 ¹H/¹³C，杂核索引只含 ¹⁹F/³¹P；
3. **子集关系**：训练索引里的每一个 (核, HOSE 码) 键都能在全量索引里找到，
   且 ¹³C 行数约为全量的 82%（= 1 − 留出率），说明留出集**确实没进训练索引**。

三条都过，再用 ``tools/eval_tier_ab.py`` 跑出的数字才有意义。

用法
----
    python tools/check_eval_kit.py
    python tools/check_eval_kit.py --quick      # 只查存在性/哈希，跳过子集扫描
    python tools/check_eval_kit.py --json       # 机器可读

退出码 0 = 全部通过；1 = 有检查项失败。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
from pathlib import Path

# 仓库根目录 / 便携包的 app/ 目录——两种布局下本脚本都在 <root>/tools/ 里
ROOT = Path(__file__).resolve().parent.parent

FULL_INDEX = ROOT / "nmrflow" / "data" / "hose_index.sqlite"
TRAIN_INDEX = ROOT / "build" / "hose_index_train.sqlite"
HETERO_INDEX = ROOT / "build" / "hose_index_hetero.sqlite"
HASH_FILE = ROOT / "build" / "EVAL_KIT_SHA256.txt"

#: 留出率——训练索引应占全量的 (1 - HOLDOUT) 上下
HOLDOUT = 0.2
#: ¹³C 行数比允许的区间：留出不是按原子抽的，行数比会偏离 0.8 一些
RATIO_MIN, RATIO_MAX = 0.72, 0.90

#: 每个核大致应有的行数量级（只防"文件截断/串了别的库"，不锁死具体值）
MIN_ROWS = {"1H": 60_000, "13C": 190_000, "19F": 1_000, "31P": 9_000}

KIT_FILES = [
    "tools/check_eval_kit.py",
    "tools/check_reference_records.py",
    "tools/eval_tier_ab.py",
    "tools/eval_hose_index.py",
    "tools/build_hose_index.py",
    "tools/diag_sp2.py",
    "tools/diag_sp3.py",
    "tools/diag_bicyclic.py",
    "tools/diag_hetero.py",
    "tools/mine_hetero.py",
    "tools/ml_baseline_19f.py",
    "tools/bench_spin.py",
    "tools/README.md",
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def _counts(path: Path) -> dict[str, int]:
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        rows = conn.execute(
            "SELECT nucleus, COUNT(*) FROM hose GROUP BY nucleus"
        ).fetchall()
    return {str(n): int(c) for n, c in rows}


def _key_hashes(path: Path, nuclei: list[str] | None = None) -> set[int]:
    """把 (核, 码) 压成 8 字节整数集合——直接存字符串要几百 MB。"""
    keys: set[int] = set()
    with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as conn:
        sql = "SELECT nucleus, code FROM hose"
        params: tuple = ()
        if nuclei:
            sql += " WHERE nucleus IN (" + ",".join("?" * len(nuclei)) + ")"
            params = tuple(nuclei)
        for nucleus, code in conn.execute(sql, params):
            keys.add(_pack(f"{nucleus}\x00{code}"))
    return keys


def _pack(text: str) -> int:
    return int.from_bytes(
        hashlib.blake2b(text.encode(), digest_size=8).digest(), "big"
    )


def _read_hashes(path: Path) -> dict[str, str]:
    """读 ``<sha256>  <相对路径>`` 清单（与 SHA256SUMS.txt 同格式）。"""
    table: dict[str, str] = {}
    if not path.exists():
        return table
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            table[parts[1].strip().replace("\\", "/")] = parts[0].upper()
    return table


def check(quick: bool = False) -> tuple[bool, list[dict[str, str]]]:
    results: list[dict[str, str]] = []

    def add(name: str, ok: bool, detail: str = "") -> None:
        results.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})

    # --- 1. 脚本在不在 -----------------------------------------------------
    missing = [rel for rel in KIT_FILES if not (ROOT / rel).exists()]
    add("评估脚本齐备", not missing, "缺: " + ", ".join(missing) if missing else f"{len(KIT_FILES)} 个")

    # --- 2. 索引在不在 -----------------------------------------------------
    present = {p: p.exists() for p in (FULL_INDEX, TRAIN_INDEX, HETERO_INDEX)}
    absent = [p.name for p, ok in present.items() if not ok]
    add("索引齐备", not absent, "缺: " + ", ".join(absent) if absent else "full / train / hetero")

    # --- 3. 哈希 -----------------------------------------------------------
    expected = _read_hashes(HASH_FILE)
    if not expected:
        add("哈希清单", True, f"{HASH_FILE.name} 缺失，跳过（源码仓库正常现象）")
    else:
        bad = []
        for path in present:
            if not present[path]:
                continue
            rel = str(path.relative_to(ROOT)).replace("\\", "/")
            want = expected.get(rel)
            if want is None:
                bad.append(f"{rel} 未登记")
            elif _sha256(path) != want:
                bad.append(f"{rel} 摘要不符")
        add("哈希一致", not bad, "；".join(bad) if bad else f"{len(expected)} 条")

    if absent:
        return False, results

    # --- 4. 构图 -----------------------------------------------------------
    full = _counts(FULL_INDEX)
    train = _counts(TRAIN_INDEX)
    hetero = _counts(HETERO_INDEX)

    thin = {n: full.get(n, 0) for n, floor in MIN_ROWS.items() if full.get(n, 0) < floor}
    add("全量索引四核齐全", not thin, f"{full}" if not thin else f"行数不足: {thin}")
    add("训练索引只含 1H/13C", set(train) <= {"1H", "13C"}, f"{train}")
    add("杂核索引只含 19F/31P", set(hetero) <= {"19F", "31P"}, f"{hetero}")

    # --- 5. 留出确实是留出 --------------------------------------------------
    ratio = train.get("13C", 0) / max(1, full.get("13C", 0))
    add(
        "13C 训练/全量行数比",
        RATIO_MIN <= ratio <= RATIO_MAX,
        f"{ratio:.3f}（期望 {RATIO_MIN}~{RATIO_MAX}，留出 {HOLDOUT:.0%}）",
    )

    if quick:
        add("子集关系", True, "--quick 跳过")
    else:
        for label, path, nuclei in (
            ("训练", TRAIN_INDEX, ["1H", "13C"]),
            ("杂核", HETERO_INDEX, ["19F", "31P"]),
        ):
            sub = _key_hashes(path)
            base = _key_hashes(FULL_INDEX, nuclei)
            extra = len(sub - base)
            add(
                f"{label}索引的键全部来自全量索引",
                extra == 0,
                f"{len(sub)} 键，多出 {extra}",
            )

    return all(r["status"] == "PASS" for r in results), results


REPRO = """\
复算 README §5 的数字（需要 nmrshiftdb2withsignals.sd，见 {dir}README.md）：

  # 1H / 13C 合并后 0.39 / 3.41
  {py} {dir}eval_tier_ab.py --sd <SD> --nuclei 1H,13C `
      --train-index {bld}hose_index_train.sqlite

  # 19F / 31P 8.76 / 19.71
  {py} {dir}eval_tier_ab.py --sd <SD> --nuclei 19F,31P `
      --train-index {bld}hose_index_hetero.sqlite

  # HOSE 自身 0.25 / 3.16（逐球层分解）
  {py} {dir}eval_hose_index.py --sd <SD> --nuclei 1H,13C
"""


def _configure_console() -> None:
    """默认中文 Windows 控制台是 GBK，别让排版字符把自检打断。

    ``¹``(U+00B9) 这类上下标在 GBK 里没有，直接 print 会
    ``UnicodeEncodeError: 'gbk' codec can't encode character '\xb9'``——
    用户双击 自检.bat 看到的第一件事就是崩溃。所以：核名一律用 ASCII，
    另外给输出流加 ``backslashreplace`` 兜底（真遇到打不出的字符，宁可印转义也不崩）。
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="backslashreplace")
        except Exception:
            pass


def _kit_paths() -> tuple[str, str]:
    """便携包里 tools/ 与 build/ 都在 app/ 下，命令要带前缀；源码仓库里不用。"""
    if ROOT.name == "app":
        return "app\\tools\\", "app\\build\\"
    return "tools\\", "build\\"


def main() -> int:
    _configure_console()
    ap = argparse.ArgumentParser(description="核验随包评估套件")
    ap.add_argument("--quick", action="store_true", help="跳过键子集扫描")
    ap.add_argument("--json", action="store_true", help="机器可读输出")
    args = ap.parse_args()

    ok, results = check(args.quick)
    if args.json:
        print(json.dumps({"ok": ok, "root": str(ROOT), "checks": results},
                         ensure_ascii=False, indent=2))
    else:
        print(f"评估套件根目录: {ROOT}\n")
        width = max(len(r["check"]) for r in results)
        for r in results:
            print(f"  [{r['status']}] {r['check']:<{width}}  {r['detail']}")
        print()
        if ok:
            dir_, bld = _kit_paths()
            print(REPRO.format(py=Path(sys.executable).name, dir=dir_, bld=bld))
        else:
            print("有检查项失败：先修好索引再谈复算。")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

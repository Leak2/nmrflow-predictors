"""便携包内的回归入口：``python tests/run_tests.py``

包内只带两个文件：``run_tests.py`` + ``test_regressions.py``（stdlib unittest，
因为便携 runtime 没有 pytest）。仓库里另有 pytest 风格的测试（依赖 fixture），
所以这里**只跑 test_regressions.py**，不做全目录 discover——否则在仓库里执行会
把 pytest 专用文件也收进来，报一堆与包无关的错误。
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
APP = ROOT / "app"
if APP.exists():                      # 便携包布局：nmrflow 在 app/ 下
    sys.path.insert(0, str(APP))
else:                                 # 源码仓库布局：nmrflow/ 就在根目录
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))         # 同目录的 test_regressions.py

if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromName("test_regressions")
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)

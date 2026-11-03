"""可选：浏览器页面全部关闭后，自动停止 Streamlit 服务。

为什么不用 beforeunload 直接杀进程
--------------------------------
浏览器不能执行本地 bat；而且刷新/跳转/崩溃都会误触发。这里用心跳：

- 每个 Streamlit 会话每 30 秒 touch 一次；
- 所有页面都关闭后，心跳停止；
- 守护线程发现「最后一个心跳」超过宽限期，才 ``os._exit(0)``。

宽限期默认 180 秒，可用环境变量 ``NMRFLOW_AUTOSTOP_GRACE`` 调整（测试用）。
多个标签页共享心跳：只要还有一个页面开着，就不会停。
"""

from __future__ import annotations

import os
import threading
import time

#: 最后一个页面关闭后等多久再退出。默认 3 分钟，给刷新/切标签/短暂断网留余量。
GRACE_SECONDS = float(os.environ.get("NMRFLOW_AUTOSTOP_GRACE", "180"))
#: 守护线程检查间隔。
CHECK_SECONDS = float(os.environ.get("NMRFLOW_AUTOSTOP_CHECK", "5"))

_enabled = False
_started = False
_fallback = time.time()
_sessions: dict[str, float] = {}
_lock = threading.Lock()


def is_enabled() -> bool:
    with _lock:
        return _enabled


def set_enabled(value: bool) -> None:
    """打开/关闭自动停止；第一次打开时启动守护线程。"""
    global _enabled, _started, _fallback
    with _lock:
        _enabled = bool(value)
        if _enabled:
            _fallback = time.time()
            if not _started:
                _started = True
                threading.Thread(
                    target=_monitor, name="nmrflow-autostop", daemon=True
                ).start()


def touch(session_id: str) -> None:
    """页面还在时调用；多个标签页各用自己的 session_id。"""
    with _lock:
        _sessions[str(session_id)] = time.time()


def latest_heartbeat() -> float:
    with _lock:
        return max([_fallback, *_sessions.values()])


def _monitor() -> None:
    while True:
        time.sleep(max(0.2, CHECK_SECONDS))
        with _lock:
            enabled = _enabled
            latest = max([_fallback, *_sessions.values()])
        if not enabled:
            continue
        if time.time() - latest > max(1.0, GRACE_SECONDS):
            try:
                print(
                    f"[autostop] no heartbeat for {GRACE_SECONDS:.0f}s; exiting.",
                    flush=True,
                )
            except Exception:
                pass
            os._exit(0)

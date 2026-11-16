"""Stop flag for a calculation running on a background thread.

The flag is local to one job. CCP4 programs are started in their own
process group so a stop can signal the shell and the program it launched.
"""

from __future__ import annotations

import contextvars
import os
import signal
import subprocess
import threading

_CURRENT: contextvars.ContextVar[Interrupt | None] = contextvars.ContextVar(
    "xfel_calculation_interrupt",
    default=None,
)


class CalculationStopped(Exception):
    """The user interrupted this calculation."""


class Interrupt:
    """A stop request for one running calculation."""

    def __init__(self) -> None:
        self.event = threading.Event()
        self._lock = threading.Lock()
        self._proc: subprocess.Popen | None = None

    def check(self) -> None:
        if self.event.is_set():
            raise CalculationStopped("The calculation was stopped.")

    def attach(self, proc: subprocess.Popen) -> None:
        with self._lock:
            self._proc = proc
            armed = self.event.is_set()
        if armed:
            signal_process_group(proc, signal.SIGTERM)

    def detach(self, proc: subprocess.Popen) -> None:
        with self._lock:
            if self._proc is proc:
                self._proc = None

    def request(self) -> None:
        """Ask the calculation to stop, and signal its current program."""
        self.event.set()
        with self._lock:
            proc = self._proc
        if proc is not None:
            signal_process_group(proc, signal.SIGTERM)


def activate(interrupt: Interrupt) -> contextvars.Token:
    return _CURRENT.set(interrupt)


def deactivate(token: contextvars.Token) -> None:
    _CURRENT.reset(token)


def current() -> Interrupt | None:
    return _CURRENT.get()


def check() -> None:
    interrupt = _CURRENT.get()
    if interrupt is not None:
        interrupt.check()


def signal_process_group(proc: subprocess.Popen, sig: int) -> None:
    """Signal a child process group, never the process group of this server."""
    if proc.poll() is not None:
        return
    try:
        pgid = os.getpgid(proc.pid)
    except ProcessLookupError:
        return
    if pgid == os.getpgid(os.getpid()):
        return
    try:
        os.killpg(pgid, sig)
    except ProcessLookupError:
        return


def terminate_process(proc: subprocess.Popen) -> None:
    """End a child process group. SIGKILL follows if it ignores SIGTERM."""
    signal_process_group(proc, signal.SIGTERM)
    try:
        proc.wait(timeout=2)
        return
    except subprocess.TimeoutExpired:
        pass
    signal_process_group(proc, signal.SIGKILL)
    try:
        proc.wait(timeout=2)
    except subprocess.TimeoutExpired:
        pass

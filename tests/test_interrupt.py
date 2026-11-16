"""Stopping a calculation kills its process group and leaves the job stopped."""

import os
import shlex
import signal
import tempfile
import threading
import time
import unittest
from pathlib import Path

from xfel_web.ccp4 import Ccp4Session, CommandFailed
from xfel_web.interrupt import (
    CalculationStopped,
    Interrupt,
    activate,
    check,
    deactivate,
)


def _ensure_dead(pid: int) -> None:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return
    try:
        os.kill(pid, signal.SIGKILL)
    except ProcessLookupError:
        return


class InterruptTests(unittest.TestCase):
    def test_check_does_nothing_until_a_calculation_is_stopped(self):
        check()
        interrupt = Interrupt()
        token = activate(interrupt)
        self.addCleanup(deactivate, token)
        check()
        interrupt.request()
        with self.assertRaises(CalculationStopped):
            check()

    def test_stream_returns_command_output(self):
        session = Ccp4Session.__new__(Ccp4Session)
        session._sink = lambda _message: None
        lines = session._stream(["/bin/echo", "hello"], cwd=Path("/"), env=None)
        self.assertEqual(lines, ["hello"])

    def test_a_failed_command_is_still_a_failure(self):
        interrupt = Interrupt()
        token = activate(interrupt)
        self.addCleanup(deactivate, token)
        session = Ccp4Session.__new__(Ccp4Session)
        session._sink = lambda _message: None
        with self.assertRaises(CommandFailed):
            session._stream(["/bin/bash", "-c", "echo nope; exit 3"], cwd=Path("/"), env=None)

    def test_stop_before_a_command_does_not_run_it(self):
        interrupt = Interrupt()
        token = activate(interrupt)
        self.addCleanup(deactivate, token)
        interrupt.request()
        session = Ccp4Session.__new__(Ccp4Session)
        session._sink = lambda _message: None
        with self.assertRaises(CalculationStopped):
            session._stream(["/bin/sleep", "30"], cwd=Path("/"), env=None)

    def test_stop_kills_the_shell_and_its_child(self):
        interrupt = Interrupt()
        token = activate(interrupt)
        self.addCleanup(deactivate, token)
        session = Ccp4Session.__new__(Ccp4Session)
        session._sink = lambda _message: None
        with tempfile.TemporaryDirectory() as folder:
            child_path = Path(folder) / "child"
            script = f"sleep 30 & echo $! > {shlex.quote(str(child_path))}; wait"

            def stop_when_ready() -> None:
                deadline = time.time() + 5
                while time.time() < deadline:
                    if child_path.is_file() and child_path.read_text().strip().isdigit():
                        interrupt.request()
                        return
                    time.sleep(0.02)
                interrupt.request()

            threading.Thread(target=stop_when_ready, daemon=True).start()
            started = time.time()
            with self.assertRaises(CalculationStopped):
                session._stream(["/bin/bash", "-c", script], cwd=Path(folder), env=None)
            self.assertLess(time.time() - started, 10)
            child = int(child_path.read_text().strip())
            self.addCleanup(_ensure_dead, child)
            deadline = time.time() + 3
            while time.time() < deadline:
                try:
                    os.kill(child, 0)
                except ProcessLookupError:
                    return
                time.sleep(0.05)
            self.fail(f"child process {child} was still running")


if __name__ == "__main__":
    unittest.main()

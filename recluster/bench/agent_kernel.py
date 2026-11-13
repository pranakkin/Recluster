"""A persistent Python process for the agent case study (run as a script, not imported).

Reads one JSON object per line on stdin ({"code": ...}), runs the code in a single namespace
that persists between calls, and answers with one JSON line ({"output": ...}). Output written
at the file-descriptor level (by C extensions as well as print) is captured, so the protocol
uses a duplicate of the original stdout that cell output never reaches.

Importing ``recluster`` is blocked so the agent cannot use the tool that later audits it.
"""

import ast
import importlib.abc
import json
import os
import sys
import tempfile
import traceback


class _BlockRecluster(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None) -> None:  # noqa: ANN001
        if name == "recluster" or name.startswith("recluster."):
            raise ModuleNotFoundError(f"No module named {name!r}")
        return None


def _run(code: str, namespace: dict) -> None:
    tree = ast.parse(code, "<cell>", "exec")
    last = None
    if tree.body and isinstance(tree.body[-1], ast.Expr):
        last = ast.Expression(tree.body[-1].value)
        tree.body.pop()
    exec(compile(tree, "<cell>", "exec"), namespace)
    if last is not None:
        value = eval(compile(last, "<cell>", "eval"), namespace)
        if value is not None:
            print(repr(value))


def main() -> None:
    sys.meta_path.insert(0, _BlockRecluster())
    channel = os.fdopen(os.dup(1), "w", buffering=1)
    namespace: dict = {"__name__": "__main__"}
    for line in sys.stdin:
        code = json.loads(line)["code"]
        with tempfile.TemporaryFile(mode="w+b") as sink:
            saved = os.dup(1), os.dup(2)
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(sink.fileno(), 1)
            os.dup2(sink.fileno(), 2)
            try:
                _run(code, namespace)
            except BaseException:  # a cell's SystemExit or KeyboardInterrupt is its own error
                traceback.print_exc()
            finally:
                sys.stdout.flush()
                sys.stderr.flush()
                os.dup2(saved[0], 1)
                os.dup2(saved[1], 2)
                os.close(saved[0])
                os.close(saved[1])
            sink.seek(0)
            output = sink.read().decode("utf-8", errors="replace")
        channel.write(json.dumps({"output": output}) + "\n")


if __name__ == "__main__":
    main()

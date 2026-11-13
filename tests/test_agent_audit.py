"""The agent case study's harness: the Python session and the tool loop, with a scripted model."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from recluster.bench.agent_audit import Kernel, render_transcript, run_agent, truncate


def test_kernel_keeps_state_and_captures_output(tmp_path: Path):
    k = Kernel(tmp_path, timeout=60)
    try:
        assert k.run("x = 21") == ""
        assert k.run("print('twice', x * 2)\nx") == "twice 42\n21\n"
        assert "ZeroDivisionError" in k.run("1 / 0")
        assert k.run("x") == "21\n"  # an error does not lose the session
        assert "No module named 'recluster'" in k.run("import recluster")
        assert k.run("import os; os.environ.get('ANTHROPIC_API_KEY')") == ""
        assert k.run("import os; os.write(1, b'fd-level\\n'); None") == "fd-level\n"
    finally:
        k.close()


def test_kernel_restarts_after_a_timeout(tmp_path: Path):
    k = Kernel(tmp_path, timeout=1)
    try:
        k.run("y = 1")
        out = k.run("import time; time.sleep(5)")
        assert "longer than 1 s" in out and "restarted" in out
        assert "NameError" in k.run("y")
    finally:
        k.close()


def test_truncate_keeps_head_and_tail():
    text = "a" * 50 + "b" * 50
    out = truncate(text, limit=20)
    assert out.startswith("a" * 10) and out.endswith("b" * 10) and "80 characters" in out


class _Block(SimpleNamespace):
    def model_dump(self, mode: str = "json") -> dict:
        return dict(self.__dict__)


class _Stream:
    def __init__(self, response: SimpleNamespace) -> None:
        self.response = response

    def __enter__(self) -> _Stream:
        return self

    def __exit__(self, *exc: object) -> None:
        return None

    def get_final_message(self) -> SimpleNamespace:
        return self.response


class _ScriptedClient:
    """Returns one tool call, then a final answer; records what it was sent."""

    def __init__(self) -> None:
        self.sent: list[list] = []
        usage = _Block(input_tokens=10, output_tokens=5)
        self.responses = [
            SimpleNamespace(
                id="m1",
                model="scripted",
                stop_reason="tool_use",
                usage=usage,
                content=[
                    _Block(type="text", text="Loading."),
                    _Block(type="tool_use", id="t1", name="python", input={"code": "print(6*7)"}),
                ],
            ),
            SimpleNamespace(
                id="m2",
                model="scripted",
                stop_reason="end_turn",
                usage=usage,
                content=[_Block(type="text", text="ISG15 is up in monocytes.")],
            ),
        ]
        self.messages = SimpleNamespace(stream=self._stream)

    def _stream(self, **kwargs: object) -> _Stream:
        self.sent.append(list(kwargs["messages"]))  # type: ignore[arg-type]
        return _Stream(self.responses[len(self.sent) - 1])


def test_loop_runs_tools_and_saves_everything(tmp_path: Path):
    client = _ScriptedClient()
    record = run_agent(client, tmp_path, tmp_path / "out", model="scripted", max_turns=5)
    out = tmp_path / "out"
    assert record["stop"] == "end_turn" and record["n_code_cells"] == 1
    assert record["usage_total"]["output_tokens"] == 10
    assert (out / "cells" / "001.py").read_text() == "print(6*7)"
    assert (out / "cells" / "001.out.txt").read_text() == "42\n"
    assert (out / "final_answer.md").read_text() == "ISG15 is up in monocytes.\n"
    transcript = json.loads((out / "transcript.json").read_text())
    assert [m["role"] for m in transcript] == ["user", "assistant", "user", "assistant"]
    assert transcript[2]["content"][0]["content"] == "42\n"
    second_call = client.sent[1]
    assert second_call[2]["content"][0]["tool_use_id"] == "t1"
    md = render_transcript(transcript)
    assert "print(6*7)" in md and "ISG15 is up" in md


def test_turn_limit_is_recorded(tmp_path: Path):
    client = _ScriptedClient()
    record = run_agent(client, tmp_path, tmp_path / "out", model="scripted", max_turns=1)
    assert record["stop"] == "turn limit" and record["n_api_calls"] == 1


@pytest.mark.parametrize("bad", [{"code": 3}, {}])
def test_malformed_tool_input_is_an_error_result(tmp_path: Path, bad: dict):
    client = _ScriptedClient()
    client.responses[0].content[1].input = bad
    run_agent(client, tmp_path, tmp_path / "out", model="scripted", max_turns=5)
    result = client.sent[1][2]["content"][0]
    assert result["is_error"] is True

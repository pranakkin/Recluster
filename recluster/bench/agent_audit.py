"""Run an off-the-shelf LLM agent on Kang 2018 and keep everything it did.

The agent is Claude through the Anthropic Messages API with one tool, a persistent Python
session in a working directory that holds the data. There is no system prompt and no
instruction about replicates, donors or statistics; the only text the agent receives is
``PROMPT`` and the tool description. The transcript, every code cell, the final answer and
the files the agent wrote are saved unedited under ``bench/agent_audit/run``.

Usage: ``uv run python -m recluster.bench.agent_audit [--max-turns N]`` (needs
ANTHROPIC_API_KEY). Nothing here is used by the rest of Recluster.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import selectors
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

MODEL = "claude-opus-5-5"
API_URL = "https://api.anthropic.com"
MAX_TOKENS = 32000
CELL_TIMEOUT_S = 1200
OUTPUT_LIMIT = 20000
INPUT_NAME = "pbmc_ifnb.h5ad"

PROMPT = f"""\
I have a single-cell RNA-seq dataset of human PBMCs in {INPUT_NAME} (raw UMI counts in .X). \
Cells from 8 donors were either left untreated or stimulated with interferon-beta. In .obs, \
`label` is the condition (ctrl or stim) and `replicate` is the donor ID.

Please analyze it: cluster the cells and annotate the cell types, compare stimulated versus \
control, and tell me which cell types change in abundance and which genes change in which \
cell types. Finish with a short list of your main findings and the numbers that support them. \
Save the AnnData with your cell-type annotation as annotated.h5ad."""

TOOL = {
    "name": "python",
    "description": (
        "Run Python code in a persistent Python 3.11 session; variables persist between calls. "
        "The working directory contains the user's files. Printed output and the value of the "
        "last expression are returned (long output is truncated). Common scientific and "
        "single-cell analysis packages are installed. There is no internet access and packages "
        "cannot be installed. Figures are not displayed; save them to files if needed."
    ),
    "input_schema": {
        "type": "object",
        "properties": {"code": {"type": "string", "description": "Python code to run"}},
        "required": ["code"],
        "additionalProperties": False,
    },
    "eager_input_streaming": True,
}

DEFAULT_OUT = Path("bench/agent_audit/run")


def prepare_input(source: Path, target: Path) -> str:
    """Raw counts with only donor and condition: the authors' labels and embeddings removed."""
    import anndata as ad

    full = ad.read_h5ad(source)
    obs = full.obs[["label", "replicate"]].copy()
    for c in obs.columns:
        obs[c] = obs[c].astype(str).astype("category")
    slim = ad.AnnData(X=full.X.copy(), obs=obs, var=full.var[[]].copy())
    target.parent.mkdir(parents=True, exist_ok=True)
    slim.write_h5ad(target)
    return hashlib.sha256(target.read_bytes()).hexdigest()


def truncate(text: str, limit: int = OUTPUT_LIMIT) -> str:
    if len(text) <= limit:
        return text
    head, tail = text[: limit // 2], text[-limit // 2 :]
    return f"{head}\n... [{len(text) - limit:,} characters truncated] ...\n{tail}"


class Kernel:
    """A persistent Python process (see agent_kernel.py) with a per-cell timeout."""

    def __init__(self, cwd: Path, timeout: float = CELL_TIMEOUT_S) -> None:
        self.cwd, self.timeout = cwd, timeout
        self.proc: subprocess.Popen[str] | None = None
        self.start()

    def start(self) -> None:
        env = {k: v for k, v in os.environ.items() if not k.startswith("ANTHROPIC_")}
        env["MPLBACKEND"] = "Agg"
        script = Path(__file__).with_name("agent_kernel.py")
        self.proc = subprocess.Popen(
            [sys.executable, "-I", "-u", str(script)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            cwd=self.cwd,
            env=env,
            text=True,
        )

    def run(self, code: str) -> str:
        assert self.proc is not None and self.proc.stdin and self.proc.stdout
        self.proc.stdin.write(json.dumps({"code": code}) + "\n")
        self.proc.stdin.flush()
        sel = selectors.DefaultSelector()
        sel.register(self.proc.stdout, selectors.EVENT_READ)
        ready = sel.select(timeout=self.timeout)
        sel.close()
        line = self.proc.stdout.readline() if ready else ""
        if not line:
            reason = (
                f"the cell ran longer than {self.timeout:.0f} s"
                if not ready
                else "the Python process exited"
            )
            self.close()
            self.start()
            return f"Error: {reason}; the Python session was restarted and its state is lost."
        return json.loads(line)["output"]

    def close(self) -> None:
        if self.proc is not None:
            self.proc.kill()
            self.proc.wait()
            self.proc = None


def _dump(content: Any) -> Any:
    if isinstance(content, str):
        return content
    return [b if isinstance(b, dict) else b.model_dump(mode="json") for b in content]


def run_agent(
    client: Any,
    workdir: Path,
    out: Path,
    *,
    model: str = MODEL,
    max_turns: int = 60,
    kernel: Kernel | None = None,
) -> dict[str, Any]:
    """The tool-use loop. Writes transcript.json, cells/, final_answer.md and run.json to out."""
    kernel = kernel or Kernel(workdir)
    cells_dir = out / "cells"
    cells_dir.mkdir(parents=True, exist_ok=True)
    messages: list[dict[str, Any]] = [{"role": "user", "content": PROMPT}]
    turns: list[dict[str, Any]] = []
    n_cells = 0
    stop = "turn limit"
    started = datetime.now(UTC).isoformat(timespec="seconds")

    def save() -> None:
        (out / "transcript.json").write_text(
            json.dumps(
                [{"role": m["role"], "content": _dump(m["content"])} for m in messages], indent=1
            )
        )

    try:
        for _ in range(max_turns):
            t0 = time.perf_counter()
            with client.messages.stream(
                model=model,
                max_tokens=MAX_TOKENS,
                tools=[TOOL],
                messages=messages,
                cache_control={"type": "ephemeral"},
            ) as stream:
                response = stream.get_final_message()
            turns.append(
                {
                    "id": response.id,
                    "model": response.model,
                    "stop_reason": response.stop_reason,
                    "usage": response.usage.model_dump(mode="json"),
                    "seconds": round(time.perf_counter() - t0, 1),
                }
            )
            messages.append({"role": "assistant", "content": response.content})
            calls = [b for b in response.content if b.type == "tool_use"]
            if response.stop_reason != "tool_use" or not calls:
                stop = response.stop_reason
                save()
                break
            results = []
            for call in calls:
                code = call.input.get("code") if isinstance(call.input, dict) else None
                if call.name != "python" or not isinstance(code, str):
                    results.append(
                        {
                            "type": "tool_result",
                            "tool_use_id": call.id,
                            "content": "Error: expected the python tool with a string `code`.",
                            "is_error": True,
                        }
                    )
                    continue
                n_cells += 1
                (cells_dir / f"{n_cells:03d}.py").write_text(code)
                t1 = time.perf_counter()
                output = kernel.run(code)
                (cells_dir / f"{n_cells:03d}.out.txt").write_text(output)
                results.append(
                    {
                        "type": "tool_result",
                        "tool_use_id": call.id,
                        "content": truncate(output) or "(no output)",
                    }
                )
                turns[-1].setdefault("cell_seconds", []).append(round(time.perf_counter() - t1, 1))
            messages.append({"role": "user", "content": results})
            save()
    finally:
        kernel.close()
        save()

    final = messages[-1]
    final_text = (
        "\n".join(b.text for b in final["content"] if getattr(b, "type", "") == "text")
        if final["role"] == "assistant" and not isinstance(final["content"], str)
        else ""
    )
    (out / "final_answer.md").write_text(final_text + "\n")
    (out / "analysis.py").write_text(
        "\n\n".join(
            f"# --- cell {i:03d} ---\n{(cells_dir / f'{i:03d}.py').read_text()}"
            for i in range(1, n_cells + 1)
        )
        + "\n"
    )
    usage = {
        k: sum(t["usage"].get(k) or 0 for t in turns)
        for k in (
            "input_tokens",
            "output_tokens",
            "cache_creation_input_tokens",
            "cache_read_input_tokens",
        )
    }
    record = {
        "model_requested": model,
        "models_returned": sorted({t["model"] for t in turns}),
        "started": started,
        "finished": datetime.now(UTC).isoformat(timespec="seconds"),
        "stop": stop,
        "n_api_calls": len(turns),
        "n_code_cells": n_cells,
        "usage_total": usage,
        "turns": turns,
    }
    (out / "run.json").write_text(json.dumps(record, indent=1))
    return record


def render_transcript(transcript: list[dict[str, Any]]) -> str:
    """A readable Markdown view of transcript.json (the JSON is the record)."""
    lines = ["# Agent transcript", ""]
    for m in transcript:
        content = m["content"]
        if isinstance(content, str):
            lines += [f"## {m['role']}", "", content, ""]
            continue
        for b in content:
            kind = b.get("type")
            if kind == "text":
                lines += [f"## {m['role']}", "", b["text"], ""]
            elif kind == "tool_use":
                lines += ["### code", "", "```python", b["input"].get("code", ""), "```", ""]
            elif kind == "tool_result":
                body = b["content"] if isinstance(b["content"], str) else json.dumps(b["content"])
                lines += ["### output", "", "```", body, "```", ""]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--data", type=Path, default=Path("data/kang_2018.h5ad"))
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--max-turns", type=int, default=60)
    p.add_argument("--model", default=MODEL)
    args = p.parse_args(argv)

    import anthropic

    from recluster.bench.kang import download_kang
    from recluster.llm.common import get_client

    get_client()  # loads .env and fails clearly without a key
    # explicit, so an ANTHROPIC_BASE_URL set by another tool cannot reroute the run
    client = anthropic.Anthropic(base_url=API_URL)
    out: Path = args.out
    if (out / "transcript.json").exists():
        sys.exit(f"{out} already holds a run; move it away first (runs are never overwritten)")
    workdir = out / "workspace"
    digest = prepare_input(download_kang(args.data), workdir / INPUT_NAME)
    (out / "prompt.md").write_text(PROMPT + "\n")
    (out / "setup.json").write_text(
        json.dumps(
            {
                "input": str(workdir / INPUT_NAME),
                "input_sha256": digest,
                "source": str(args.data),
                "api_url": API_URL,
                "tool": TOOL,
                "system_prompt": None,
                "max_tokens": MAX_TOKENS,
                "max_turns": args.max_turns,
                "cell_timeout_s": CELL_TIMEOUT_S,
                "anthropic_sdk": anthropic.__version__,
                "python": platform.python_version(),
                "platform": platform.platform(),
            },
            indent=1,
        )
    )
    record = run_agent(client, workdir, out, model=args.model, max_turns=args.max_turns)
    (out / "transcript.md").write_text(
        render_transcript(json.loads((out / "transcript.json").read_text()))
    )
    print(json.dumps({k: v for k, v in record.items() if k != "turns"}, indent=1))


if __name__ == "__main__":
    main()

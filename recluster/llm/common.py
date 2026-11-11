"""Shared plumbing for the optional LLM layer.

Everything in ``recluster.llm`` is optional. The rest of the tool never imports it, and a
missing API key or SDK produces a clear message pointing at the hand-written YAML route.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

DEFAULT_MODEL = "claude-opus-5-5"
MODEL_ENV = "RECLUSTER_LLM_MODEL"
KEY_ENV = "ANTHROPIC_API_KEY"


class LLMUnavailable(RuntimeError):
    """Raised when the LLM layer cannot run; the message says how to proceed without it."""


def load_dotenv(path: str | Path = ".env") -> None:
    """Load ``KEY=value`` lines from a local .env into the environment (no override)."""
    p = Path(path)
    if not p.exists():
        return
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        key, value = key.strip(), value.strip().strip("'\"")
        if key.startswith("ANTHROPIC_") or key == MODEL_ENV:
            os.environ.setdefault(key, value)


def model_name(override: str | None = None) -> str:
    return override or os.environ.get(MODEL_ENV) or DEFAULT_MODEL


def get_client() -> Any:
    """An Anthropic client, or a clear LLMUnavailable error."""
    load_dotenv()
    try:
        import anthropic
    except ImportError as exc:  # pragma: no cover - depends on the extra being installed
        raise LLMUnavailable(
            "the anthropic SDK is not installed; install with `uv sync --extra llm` "
            "(or write the claim YAML by hand and use `recluster run`)"
        ) from exc
    if not os.environ.get(KEY_ENV) and not os.environ.get("ANTHROPIC_AUTH_TOKEN"):
        raise LLMUnavailable(
            f"{KEY_ENV} is not set (put it in the environment or in a .env file, which is "
            "gitignored). The LLM layer is optional: write the claim YAML by hand and run "
            "`recluster run --data X.h5ad --claim claim.yaml`."
        )
    return anthropic.Anthropic()


def text_of(response: Any) -> str:
    """Concatenate the text blocks of a Messages API response."""
    parts = [b.text for b in getattr(response, "content", []) if getattr(b, "type", "") == "text"]
    return "\n".join(parts).strip()


# a number is digits (with an optional decimal part) not glued to letters: "ISG15" and "log2"
# contain no numbers, "7.3." at the end of a sentence contains 7.3, "6.2e-61" is one number
_NUM = re.compile(r"(?<![A-Za-z0-9_])\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?![A-Za-z0-9_])")


def numbers_in(text: str) -> set[str]:
    """Every number token in ``text``, normalised (``1.0`` -> ``1``, ``0.50`` -> ``0.5``)."""
    out = set()
    for m in _NUM.findall(text):
        try:
            f = float(m)
        except ValueError:
            continue
        out.add(f"{f:g}")
    return out

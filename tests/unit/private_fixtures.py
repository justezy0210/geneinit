"""Optional local regression data; synthetic tests do not depend on these files."""

import json
from pathlib import Path

import pytest


DATA = Path(__file__).resolve().parents[1] / "data"


def load_private_json(name):
    try:
        text = (DATA / name).read_text()
    except FileNotFoundError:
        pytest.skip(f"Private regression fixture unavailable: {name}")
    return json.loads(text)


def private_pairs(name):
    try:
        text = (DATA / name).read_text()
    except FileNotFoundError:
        return [pytest.param(
            None,
            marks=pytest.mark.skip(reason=f"Private regression fixture unavailable: {name}"),
            id="private-data-unavailable",
        )]
    return [pytest.param(pair, id=pair["case"]) for pair in json.loads(text)["pairs"]]

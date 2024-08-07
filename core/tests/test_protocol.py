import re

import pytest
from matplotlib.testing.decorators import image_comparison
from matplotlib.testing.exceptions import ImageComparisonFailure

import core

protocol_types = ["original", "manual"]


@pytest.mark.parametrize(
    "protocol_type, remove_false_nodes, known_lines, known_verts",
    zip(protocol_types, [True, False], [25, 16], [14, 8], strict=True),
    ids=protocol_types,
)
@pytest.mark.xfail(
    reason="Needs updated Manual & Parenex -- See GH#131.",
    raises=AssertionError,
)
def test_generate_case(protocol_type, remove_false_nodes, known_lines, known_verts):
    _case = core.protocol.protocol_cases["04"]

    city = _case["city"]

    in_data = core.utils.read_original(city)
    if protocol_type != "original":
        in_data = core.utils.read_manual(city, in_data.crs)

    observed_lines, observed_verts = core.protocol.generate_case(
        in_data,
        _case["coordinates"],
        _case["buffer"],
        remove_false_nodes=remove_false_nodes,
    )

    assert known_lines == observed_lines.shape[0]
    assert in_data.crs == observed_lines.crs
    assert known_verts == observed_verts.shape[0]
    assert in_data.crs == observed_verts.crs


def test_process_case_invalid_type():
    __case = "01"
    _case = core.protocol.protocol_cases[__case]

    with pytest.raises(
        ValueError,
        match=re.escape(
            "Invalid protocol type passed in: protocol_type='invalid'. "
            "Must be in ['original', 'manual']."
        ),
    ):
        core.protocol.process_case(
            __case,
            "invalid",
            _case["city"],
            _case["coordinates"],
            _case["buffer"],
            _case["title"],
        )


@pytest.mark.xfail(
    reason="Needs updated Manual & Parenex -- See GH#131.",
    raises=ImageComparisonFailure,
)
@image_comparison(
    baseline_images=["test_process_case_original.png"],
    style="mpl20",
    tol=5.1 if pytest.ENV_TYPE == "dev" else 1,
)
def test_process_case_original():
    __case = "04"
    _case = core.protocol.protocol_cases[__case]

    core.protocol.process_case(
        __case,
        "original",
        _case["city"],
        _case["coordinates"],
        _case["buffer"],
        _case["title"],
    )


@image_comparison(
    baseline_images=["test_process_case_manual.png"],
    style="mpl20",
    tol=7.4 if pytest.ENV_TYPE == "dev" else 1,
)
def test_process_case_manual():
    __case = "04"
    _case = core.protocol.protocol_cases[__case]

    core.protocol.process_case(
        __case,
        "manual",
        _case["city"],
        _case["coordinates"],
        _case["buffer"],
        _case["title"],
    )

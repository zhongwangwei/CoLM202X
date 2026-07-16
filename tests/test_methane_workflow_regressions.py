from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_run_configuration_changes_trigger_ci() -> None:
    workflow = (ROOT / ".github" / "workflows" / "build_CoLM_gnu.yml").read_text(
        encoding="utf-8"
    )

    assert "- 'run/**'" not in workflow

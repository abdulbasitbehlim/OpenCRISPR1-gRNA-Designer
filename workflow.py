"""Design identity, state reset and reproducible export helpers shared by the UI."""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json

from opencrispr_designer import guide_row, guide_format_variants
from validation import validate_opencrispr_guide, validation_summary_row

APP_VERSION = "1.4.0"


def guide_key(record, guide) -> str:
    payload = [record.sequence_sha256, record.source_record_version, guide.segment_id,
               guide.start, guide.end, guide.strand, guide.spacer, guide.pam]
    return hashlib.sha256(json.dumps(payload).encode()).hexdigest()


def clear_design_state(state) -> None:
    for key in list(state):
        if key in {"oc_record", "oc_guides", "oc_settings", "oc_screens", "oc_panel_config"} or key.startswith("oc_local::"):
            del state[key]


def reset_panel_if_changed(state, raw, radius, intended, display_limit) -> None:
    key = hashlib.sha256(json.dumps([raw, radius, asdict(intended) if intended else None, display_limit]).encode()).hexdigest()
    if state.get("oc_panel_config") != key:
        state["oc_screens"] = {}
        state["oc_panel_config"] = key


def export_bundle(record, guides, screens, settings):
    rows = []
    details = []
    for guide in guides:
        screen = screens.get(guide_key(record, guide))
        report = validate_opencrispr_guide(
            guide, min_sequence_score=settings["min_score"],
            min_specificity_review=settings["min_spec_review"], panel_screen=screen)
        row = {**guide_row(guide), **validation_summary_row(report)}
        row.update({"MIT specificity": screen.specificity_score if screen else None,
                    "Specificity method": screen.method if screen else "Not screened",
                    "Off-target hits": screen.total_hits if screen else None,
                    "Reference SHA-256": screen.panel_sha256 if screen else None,
                    "Reference bp": screen.panel_bp if screen else None,
                    "Mismatch radius": screen.max_mismatches if screen else None,
                    "Intended locus excluded": screen.intended_site_excluded if screen else None,
                    "Displayed hits truncated": screen.truncated if screen else None})
        rows.append(row)
        details.append({**row, "guide_id": guide_key(record, guide),
                        "expression_formats": [asdict(x) for x in guide_format_variants(guide.spacer)],
                        "validation_checks": [asdict(c) for c in report.checks],
                        "panel_screen": screen.to_dict() if screen else None})
    return rows, {"app": "OpenCRISPR-1 gRNA Designer", "version": APP_VERSION,
                  "target": record.provenance_dict(), "settings": settings,
                  "coordinate_system": "1-based inclusive within each supplied segment; not chromosome coordinates",
                  "evidence_status": "Sequence compatibility and heuristic scores are not experimentally calibrated efficacy predictions.",
                  "guides": details}

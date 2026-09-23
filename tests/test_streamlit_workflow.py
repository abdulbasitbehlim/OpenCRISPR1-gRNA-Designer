# ============================================================================
# TEST STREAMLIT WORKFLOW
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: This file contains automated checks for the OpenCRISPR-1 gRNA
# Designer. The tests help make sure scientific and interface behaviour does
# not change accidentally when the program is edited.
#
# HOW TO READ A TEST:
# 1. Prepare sample input or a reusable fixture.
# 2. Call the function or workflow being tested.
# 3. Use assertions to compare the actual result with the expected result.
# 4. Existing test logic and expected scientific outcomes are unchanged.
#
# MAIN TOP-LEVEL TESTS / HELPERS:
# - function: label
# - function: manual_app
# - function: test_screen_result_persists_then_clears_when_reference_or_design_changes
# - function: test_verified_coordinate_updates_validation_and_setting_change_hides_exports
# - function: test_offline_example_loads_and_designs_without_network
# ============================================================================

"""Exercise real Streamlit reruns, rather than inspecting UI source strings."""
from pathlib import Path
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parents[1] / 'app.py')
LOCUS = 'ACGTACGTACGTACGTACGTAGG'


def label(elements, name):
    return next(x for x in elements if x.label == name)


def manual_app():
    at=AppTest.from_file(APP, default_timeout=20).run()
    at.radio[0].set_value('Manual sequence / FASTA').run()
    label(at.text_area,'DNA / RNA / FASTA').set_value('>target\n'+LOCUS)
    label(at.button,'Design OpenCRISPR guides').click().run()
    assert not at.exception
    return at


def test_screen_result_persists_then_clears_when_reference_or_design_changes():
    at=manual_app()
    label(at.text_area,'Reference FASTA').set_value('>other\n'+LOCUS)
    label(at.button,'Run local specificity').click().run()
    assert not at.exception
    results=list(at.session_state['oc_screens'].values())
    assert len(results) == 1 and results[0].total_hits == 1
    at.run()
    assert at.session_state['oc_screens']
    label(at.text_area,'Reference FASTA').set_value('>changed\n'+LOCUS).run()
    assert at.session_state['oc_screens'] == {}
    label(at.button,'Run local specificity').click().run()
    assert at.session_state['oc_screens']
    label(at.text_area,'DNA / RNA / FASTA').set_value('>new_target\n'+LOCUS)
    label(at.button,'Design OpenCRISPR guides').click().run()
    assert not at.exception
    assert at.session_state['oc_screens'] == {}
    # A failed design cannot expose the previous successful design as its result.
    label(at.text_area,'DNA / RNA / FASTA').set_value('>x\nACGT\n>x\nTTTT')
    label(at.button,'Design OpenCRISPR guides').click().run()
    assert not at.exception
    assert 'oc_guides' not in at.session_state
    assert any('Duplicate FASTA' in x.value for x in at.error)


def test_verified_coordinate_updates_validation_and_setting_change_hides_exports():
    at=manual_app()
    label(at.text_area,'Reference FASTA').set_value('>target\n'+LOCUS).run()
    label(at.checkbox,'Exclude a verified intended locus').check().run()
    label(at.text_input,'Intended reference FASTA ID').set_value('target')
    label(at.button,'Run local specificity').click().run()
    assert not at.exception
    result=list(at.session_state['oc_screens'].values())[0]
    assert result.intended_site_excluded and result.total_hits == 0
    label(at.slider,'Maximum mismatches').set_value(2).run()
    assert at.session_state['oc_screens'] == {}
    label(at.slider,'Maximum guides').set_value(25).run()
    assert any('Design settings changed' in x.value for x in at.warning)
    assert not at.get('download_button')


def test_offline_example_loads_and_designs_without_network():
    at=AppTest.from_file(APP, default_timeout=20).run()
    label(at.button,'Load offline example').click().run()
    assert at.radio[0].value == 'Manual sequence / FASTA'
    label(at.button,'Design OpenCRISPR guides').click().run()
    assert not at.exception
    assert {g.strand for g in at.session_state['oc_guides']} == {'+', '-'}

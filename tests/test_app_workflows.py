# ============================================================================
# TEST APP WORKFLOWS
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Contains automated tests that protect the Plant MultiGene gRNA Designer workflow from accidental behaviour changes.
#
# HOW TO READ THIS FILE:
# 1. Tests first prepare an input or fixture.
# 2. The relevant program function is called.
# 3. Assertions check that the result still matches the expected behaviour.
# 4. Test logic and expected scientific results are intentionally unchanged.
#
# MAIN TOP-LEVEL PARTS:
# - function: manual_app
# - function: test_saved_settings_and_panel_evidence_survive_correctly
# - function: test_zero_guides_still_allows_custom_check_and_fallback
# - function: test_empty_panel_errors_and_failed_submission_clears_previous_run
# ============================================================================

"""Exercise actual Streamlit reruns, not just source markers."""
from pathlib import Path
from streamlit.testing.v1 import AppTest
S='ACGTACGTACGTACGTACGA'
APP=str(Path(__file__).resolve().parents[1]/'app.py')



# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: manual_app
# ----------------------------------------------------------------------------
def manual_app(raw):
    at=AppTest.from_file(APP, default_timeout=20).run()
    at.radio[0].set_value('Manual multi-FASTA').run()
    next(w for w in at.text_area if w.label=='One FASTA record per gene').set_value(raw)
    next(w for w in at.button if w.label=='Design shared guides').click().run()
    assert not at.exception
    return at



# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: test_saved_settings_and_panel_evidence_survive_correctly
# ----------------------------------------------------------------------------
def test_saved_settings_and_panel_evidence_survive_correctly():
    at=manual_app('>A\n'+S+'AGG\n>B\n'+S+'TGG')
    assert at.session_state['plant_snapshot']['settings']['max_mismatches_per_gene']==2
    next(w for w in at.slider if w.label=='Max mismatches per target gene').set_value(0).run()
    assert at.session_state['plant_snapshot']['settings']['max_mismatches_per_gene']==2
    assert any('Settings have changed' in w.value for w in at.warning)
    next(w for w in at.text_area if w.label=='Reference panel FASTA').set_value('>panel\n'+S+'AGG')
    next(w for w in at.button if w.label=='Screen selected guide').click().run()
    assert not at.exception
    assert any('1 total hits' in m.value for m in at.markdown)
    next(w for w in at.text_area if w.label=='Reference panel FASTA').set_value('>panel\nAAAA').run()
    assert any('No saved screen matches' in c.value for c in at.caption)



# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: test_zero_guides_still_allows_custom_check_and_fallback
# ----------------------------------------------------------------------------
def test_zero_guides_still_allows_custom_check_and_fallback():
    at=manual_app('>A\n'+S+'AGA\n>B\n'+'A'*35)
    assert at.session_state['plant_guides']==[]
    custom=next(w for w in at.text_input if w.label=='Custom spacer (20 nt, no PAM)')
    custom.set_value(S)
    next(w for w in at.button if w.label=='Validate custom guide').click().run()
    assert not at.exception
    assert any(m.label=='Custom guide validation' and m.value=='FAIL' for m in at.metric)



# ----------------------------------------------------------------------------
# TEST / HELPER SECTION: test_empty_panel_errors_and_failed_submission_clears_previous_run
# ----------------------------------------------------------------------------
def test_empty_panel_errors_and_failed_submission_clears_previous_run():
    at=manual_app('>A\n'+S+'AGG\n>B\n'+S+'TGG')
    next(w for w in at.button if w.label=='Screen selected guide').click().run()
    assert any('nonempty' in e.value for e in at.error)
    next(w for w in at.text_area if w.label=='One FASTA record per gene').set_value('>A\nAAA\n>A\nCCC')
    next(w for w in at.button if w.label=='Design shared guides').click().run()
    assert any('Duplicate' in e.value for e in at.error)
    assert 'plant_guides' not in at.session_state

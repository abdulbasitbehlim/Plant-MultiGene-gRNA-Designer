# ============================================================================
# CONFTEST
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Defines shared pytest setup used by the Plant MultiGene test suite.
#
# HOW TO READ THIS FILE:
# 1. Tests first prepare an input or fixture.
# 2. The relevant program function is called.
# 3. Assertions check that the result still matches the expected behaviour.
# 4. Test logic and expected scientific results are intentionally unchanged.
# ============================================================================

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# ============================================================================
# CONFTEST
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
# ============================================================================

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

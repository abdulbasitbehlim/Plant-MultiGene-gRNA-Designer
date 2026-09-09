# External-tool comparison protocol

External tools compared: CRISPys, CRISPR MultiTargeter.

A valid head-to-head comparison must use the **same frozen input sequences/assembly, PAM model, spacer length, target segments and off-target scope**. Export each external tool's guide list to CSV with a `guide` column, then use `concordance.py` to calculate:
- guides found by both;
- guides only this project finds;
- guides only the comparator finds.

Differences must then be classified (PAM model, strand handling, exon-boundary policy, mismatch/consensus strategy, activity score threshold, off-target filtering, or database/version differences).

**Build status:** external web/CLI execution was not reproducibly available in the artifact environment, therefore v1.3.1 does not invent a CHOPCHOP/MultiTargeter/CRISPys concordance count. The protocol and machine-readable comparison utility are included so the exact run can be committed before publication.

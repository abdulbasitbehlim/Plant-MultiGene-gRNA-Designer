# Benchmark evidence resolution - v1.3.1

## Question
Was `CTCTACTTTCTCCCTCATCT` proven to be a single guide from the 5,635-guide Multi-Knock transportome library assigned to PUP7/PUP8/PUP21?

## Evidence checked

1. **Hu et al. 2023, Figure 5 legend:** CR8/21 is the PUP8/PUP21 double-mutant line, and its sgRNA was explicitly not designed to target PUP7 because PUP7 does not contain the corresponding PAM for that guide.
2. **Hu et al. 2023, CRISPR/Cas9 cloning Methods:** `CTCTACTTTCTCCCTCATCT` is explicitly described as a 20-nt protospacer picked to target PUP7, PUP8 and PUP21 at once. This appears in the follow-up validation-cloning section, not as a demonstrated transportome-library assignment.
3. **Nature article landing page:** Supplementary Data 1 is described as the all-family sgRNA library. The binary assignment table was not directly retrievable in this artifact environment.

## Resolution

The previous v1.3 wording was too strong. v1.3.1 treats `CTCTACTTTCTCCCTCATCT` as a **published three-gene validation-construct spacer**, while treating the **CR8/21 screen guide as a distinct PUP8/PUP21 library guide**.

We do **not** state that the validation spacer is a verified member of the 5,635-guide transportome library. We also do **not** state that the CR7/8/21 triple mutant necessarily required two guides: the Methods themselves describe one 20-nt validation protospacer intended to target all three genes at once. Direct Supplementary Data 1 inspection is required to settle library membership.

## Terminology policy

- `exact_spacer_recovered = true` only when the tool emits the identical 20-nt published spacer.
- Per-gene mismatch tolerance is a separate target-compatibility property.
- `library_guide_rediscovered` may be used only when both the exact spacer sequence and intended gene assignment have been verified from the published library data.

## Manuscript-safe wording

> In a locked regression using PUP7/PUP8/PUP21 target contexts, the designer regenerated the exact 20-nt spacer reported by Hu et al. for a follow-up three-gene validation construct. The generated spacer was rank 1 in the fixture, with two exact target matches and one one-mismatch target. Because the publisher's all-family library assignment table was not directly verified in this build, this result is not presented as rediscovery of a 5,635-guide transportome-library entry.

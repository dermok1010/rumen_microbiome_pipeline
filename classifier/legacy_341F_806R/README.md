# Superseded: 341F/806R classifier

These two scripts built the classifier this pipeline originally shipped
with -- trained on reference reads extracted with 341F/806R primers.

That was wrong: the actual PCR amplification of these samples used
515F/806R primers. Classifying against a 341F/806R-trained classifier
produced fewer confident genus-level calls than the region actually
sequenced supports (observed directly: the resulting genus table had far
fewer retained genera than the corrected run on the same ASVs).

Kept here for provenance only -- not part of the active pipeline. The
correct classifier build is `classifier/build/`, and
`nextflow.config`'s `params.classifier` default points at its output.

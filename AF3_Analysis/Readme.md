# Cohesin and condensin hinge analysis

These scripts calculate confidence, interface, contact, energetic and conformational-state metrics from AlphaFold 3 hinge models. They were written for four species (human, mouse, cerevisiae, and pombe) and assume chain A is Smc1/SMC2 and chain B is Smc3/SMC4.

Included scripts

hinge_pipeline.py — the main analysis pipeline. For every model it calculates AlphaFold confidence metrics, ipSAE, buried surface area (BSA), interface residues, heavy-atom contacts, hydrogen-bond-like contacts, salt bridges, PRODIGY ΔG/Kd, and interface spread/elongation for conformational-state classification.

ipsae_lib.py — implementation of the published Dunbrack ipSAE interface-confidence metric. It is imported by hinge_pipeline.py and must remain in the same directory.

cb_contacts.py — calculates the canonical inter-chain Cβ-Cβ contact map at ≤8 Å (Cα is used for glycine), together with the corresponding 5 Å count. This is the residue-level, comparatively size-unbiased contact metric used for contact-density comparisons.

contact_sensitivity.py — compares 5 Å heavy-atom contacts with residue-pair definitions at 5 and 8 Å. Include and run this script when the contact-cutoff sensitivity figure or supplementary analysis is reported; otherwise it is optional.

northsouth.py — analyses the retained contact patch in open models for north/south open-substate classification. Include and run it when the north/south open-state analysis is reported; otherwise it is optional.

hinge_report.py — the reporting and plotting script. It performs the state-matched comparison, calculates per-species medians and Kruskal-Wallis statistics, and generates the summary tables and figures.

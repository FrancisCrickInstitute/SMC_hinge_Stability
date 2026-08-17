Cohesin and condensin hinge analysis

These scripts calculate confidence, interface, contact, energetic and conformational-state metrics from AlphaFold 3 hinge models. They were written for four species (human, mouse, cerevisiae, and pombe) and assume chain A is Smc1/SMC2 and chain B is Smc3/SMC4.

Keep all .py files in the same directory because the main scripts import hinge_pipeline.py and ipsae_lib.py locally.

Expected AlphaFold 3 layout

The fast drivers expect files arranged as:

<root>/<species>/results/seed_<N>/<complex>/seed-<N>_sample-<k>/
    <complex>_seed-<N>_sample-<k>_model.cif
    <complex>_seed-<N>_sample-<k>_confidences.json
    <complex>_seed-<N>_sample-<k>_summary_confidences.json

Here, <complex> is <species>_cohesin_hinge_harmonized or <species>_condensin_hinge_harmonized. Seeds are numbered from 1 and diffusion samples from 0.

Main analysis

The supplied root must contain only one complex type; analyse cohesin and condensin from separate roots and write them to separate output files:

python hinge_pipeline.py --root /path/to/cohesin_results --out cohesin_metrics.csv --budget 600
python hinge_pipeline.py --root /path/to/condensin_results --out condensin_metrics.csv --budget 600

--budget is the maximum runtime in seconds. The scripts append to the CSV and are resumable, so repeat the command until all models have been processed.

For the fixed split-per-seed layout, the faster parallel drivers are:

python fast_rerun_par.py --root /path/to/cohesin_root --out cohesin_metrics.csv --budget 600 --workers 4
python fast_rerun_condensin.py --root /path/to/condensin_root --out condensin_metrics.csv --budget 600 --workers 4

fast_rerun.py is the slower serial cohesin alternative. fast_rerun.py.tmp.92414.86a85ee8fd85 is an obsolete temporary copy and should not be used.

The main CSV contains one row per model, including AF3 confidence scores, ipSAE, buried surface area, 5 Å heavy-atom contacts, interface residues, hydrogen-bond-like contacts, salt bridges, interface spread/elongation, and PRODIGY ΔG/Kd.

Contact analyses

Run the canonical residue-contact calculation separately:

python cb_contacts.py --root /path/to/cohesin_root --out cohesin_cb_contacts.csv --budget 600 --workers 8
python cb_contacts_condensin.py --root /path/to/condensin_root --out condensin_cb_contacts.csv --budget 600 --workers 8

These scripts count inter-chain Cβ-Cβ contacts at 5 and 8 Å, using Cα for glycine. For a direct comparison of the original 5 Å heavy-atom definition with residue-pair definitions, run:

python contact_sensitivity.py --root /path/to/cohesin_root --out contact_sensitivity.csv --budget 600 --workers 8

The contact scripts are also resumable.

Reports and optional conformational analysis

Generate summary tables and figures from a main metrics CSV with:

python hinge_report.py --master cohesin_metrics.csv --outdir cohesin_report

This writes master_metrics.csv, closed_reliable_medians.csv, SUMMARY.md, and three PNG figures. The current report script applies an ipSAE-based reliability selection before its closed-state comparison; check that this matches the final analysis protocol before using the figures in a manuscript. Although --exp is accepted, the supplied script does not currently use that file.

To explore the retained interface patch among models classified as open:

python northsouth.py --master cohesin_metrics.csv --root /path/to/cohesin_root --out northsouth.csv --budget 600

Important limits

fast_rerun_par.py, fast_rerun_condensin.py, cb_contacts.py, cb_contacts_condensin.py, and contact_sensitivity.py assume 50 seeds × 5 samples. For the 500-seed condensin analysis, use hinge_pipeline.py or change the hard-coded seed limit to 500.

Species names, complex names, chain IDs and directory structure are hard-coded in several scripts.

Existing output CSVs are appended to. Use a new filename (or deliberately remove the old CSV) when analysing a different dataset.

The open/closed classification used by hinge_report.py is based on interface spread: <13 Å is open and ≥13 Å is closed.

PRODIGY energies and AlphaFold confidence scores are computational estimates; AF3 sampling frequencies should not be interpreted as equilibrium conformational populations.

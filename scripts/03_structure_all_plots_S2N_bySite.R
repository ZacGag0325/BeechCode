The revised script became too strict and now rejects all STRUCTURE Q files.

Current console output:
  
  [03_structure] Preparing final STRUCTURE individual barplots...
[03_structure] Searching folder: /Users/zacharygagnon/Desktop/BeechCode/data/structure/Q_Files
[03_structure] Ignored helper files: 1
[03_structure] Saved run inventory: /Users/zacharygagnon/Desktop/BeechCode/outputs/tables/supplementary/structure_run_inventory.csv
[03_structure] No usable K>=2 Q runs after validation.
[03_structure] Done.

This means the script no longer makes any figures at all. Please fix this properly.

Main goal:
  Make the script robust enough to correctly detect and parse valid STRUCTURE Q files, while still validating them safely. Right now the validation/parser is over-rejecting all runs.

Please revise:
  scripts/03_structure_all_plots_S2N_bySite.R

What I need:
  
  1) Diagnose why every run is being rejected
For every candidate Q file, print:
  - filename
- inferred K
- number of rows
- number of columns
- which columns are being treated as Q columns
- whether those columns are numeric after parsing
- min and max of Q columns
- row-sum range
- exact reason the file is accepted or rejected

2) Relax validation enough to accept real STRUCTURE files
The script should:
  - robustly parse whitespace-delimited Q files
- safely coerce Q columns to numeric
- ignore non-Q helper columns if present
- accept files if they contain the correct number of individuals and K ancestry columns that are valid proportions
- allow tiny floating-point tolerance in row sums and bounds
- not reject a whole run because of harmless formatting quirks

3) Keep validation informative, not destructive
- do not silently reject everything
- if a file is rejected, state the exact reason
- if at least one valid run exists for a K, use it
- if multiple valid runs exist, keep current logic for selecting the preferred run

4) Restore output generation
The script must again produce:
  - one standalone figure per K
- one combined all-K figure
- run inventory CSV

5) Keep my requested figure behavior
- order sites from south to north, left to right
- group individuals by site
- keep the same individual order across all K
- save outputs to outputs/figures/
  
  6) Return the FULL corrected script from line 1 to end
Do not give partial patches only.
Also explain briefly:
  - why all runs were being rejected
- what you changed to make parsing/validation robust again

Important:
  Do not just remove validation. Keep validation, but make it realistic for my actual STRUCTURE Q files.
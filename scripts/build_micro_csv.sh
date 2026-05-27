#!/usr/bin/env bash

set -euo pipefail

# build_micro_csv.sh -- collect microscale Fortran run parameters into one CSV
#
# Scans a directory tree of completed *microscale* Fortran runs and extracts
# the parameters each run was invoked with, writing them to a single CSV for
# side-by-side comparison. Each run's parameters are recovered from its
# captured Fortran stdout log: micro_rates echoes its command-line arguments
# as lines like "command arg 5 = --nodes" / "command arg 6 = 5", and prints
# the seed it used as "seed= <value>". This is handy for auditing or
# recovering the inputs of legacy runs whose parameters live only in their
# log files.
#
# Usage:
#   ./build_micro_csv.sh SIM_DIRECTORY OUTPUT.csv
#
# Example:
#   ./build_micro_csv.sh ./simulations micro_runs.csv
#
# Input  (SIM_DIRECTORY): a directory with one subfolder per run, each holding
#   a Fortran microscale log *.txt. Subfolders named "macro*" are skipped
#   (those are macroscale runs, not microscale). Expects one log .txt per run
#   folder; the folder's name becomes that run's column header.
#
# Output (OUTPUT.csv): one row per parameter (fixed order; see PARAM_ORDER),
#   one column per run folder. Column 1 is the parameter name and the header
#   row lists the folder names. Missing values are left blank.
#
# Notes:
#   * Fortran CLI names are translated to their lysis/Python names via NAME_MAP
#     (e.g. radius -> fiber_radius, nodes -> nodes_in_micro_row,
#     simulations -> micro_simulations). See
#     docs/source/usage/fortran_microscale.rst for the full parameter reference.
#   * Canonical units (UNITS) are appended to values where applicable.
#   * micro_seed is read from the log's "seed=" line -- the seed actually used,
#     which matters when it was randomly drawn (--seed 0).
#   * runCode / outFileCode and any unmapped arguments are ignored.
#   * Requires bash 4+ (associative arrays), plus find, awk, and grep.

ROOT_DIR="$1"
OUTPUT_CSV="$2"

[[ ! -d "$ROOT_DIR" ]] && {
    echo "Error: directory not found: $ROOT_DIR"
    exit 1
}

# ----------------------------------
# Fortran name -> Python name
# ----------------------------------
declare -A NAME_MAP=(
    [radius]="fiber_radius"
    [KdtPAyesplg]="diss_const_tPA_wPLG"
    [KdtPAnoplg]="diss_const_tPA_woPLG"
    [KdPLGintact]="diss_const_PLG_intact"
    [KdPLGnicked]="diss_const_PLG_nicked"
    [ktPAon]="bind_rate_tPA"
    [kplgon]="bind_rate_PLG"
    [freeplg]="conc_free_PLG"
    [kdeg]="deg_rate_fibrin"
    [kplioff]="unbind_rate_PLi"
    [kapcat]="activation_rate_PLG"
    [kncat]="exposure_rate_binding_site"
    [nodes]="nodes_in_micro_row"
    [snap_proportion]="snap_proportion"
    [simulations]="micro_simulations"
    [seed]="micro_seed"
)

# ----------------------------------
# Define parameter order for output CSV
# ----------------------------------
PARAM_ORDER=(
    fiber_radius
    diss_const_tPA_wPLG
    diss_const_tPA_woPLG
    diss_const_PLG_intact
    diss_const_PLG_nicked
    bind_rate_tPA
    bind_rate_PLG
    conc_free_PLG
    deg_rate_fibrin
    unbind_rate_PLi
    activation_rate_PLG
    exposure_rate_binding_site
    nodes_in_micro_row
    snap_proportion
    micro_simulations
    micro_seed
)

# ----------------------------------
# Define parameter units
# ----------------------------------
declare -A UNITS=(
    [fiber_radius]="microns"
    [diss_const_tPA_wPLG]="micromolar"
    [diss_const_tPA_woPLG]="micromolar"
    [diss_const_PLG_intact]="micromolar"
    [diss_const_PLG_nicked]="micromolar"
    [bind_rate_tPA]="(micromolar*sec)^-1"
    [bind_rate_PLG]="(micromolar*sec)^-1"
    [conc_free_PLG]="micromolar"
    [deg_rate_fibrin]="sec^-1"
    [unbind_rate_PLi]="sec^-1"
    [activation_rate_PLG]="sec^-1"
    [exposure_rate_binding_site]="sec^-1"
    [nodes_in_micro_row]=""
    [snap_proportion]=""
    [micro_simulations]=""
    [micro_seed]=""
)

# Storage:
# DATA["folder|parameter"] = value
declare -A DATA

# Ordered list of experiment columns
FOLDERS=()

# ----------------------------------
# Process simulation files
# Skip any folder named macro*
# ----------------------------------
while IFS= read -r txtfile; do

    declare -A CURRENT=()

    # Use parent folder name as unique column header
    folder_name="$(basename "$(dirname "$txtfile")")"

    FOLDERS+=("$folder_name")

    # ----------------------------------
    # Extract command argument pairs
    #
    # Example:
    # command arg 5 = --nodes
    # command arg 6 = 5
    #
    # command arg 7 = --kdeg
    # command arg 8 = 50
    # ----------------------------------
    while IFS='|' read -r raw_key raw_value; do

        # Remove leading --
        key="${raw_key#--}"

        # Skip metadata args
        [[ "$key" == "runCode" ]] && continue
        [[ "$key" == "outFileCode" ]] && continue

        mapped="${NAME_MAP[$key]:-}"

        # Skip unmapped parameters
        [[ -z "$mapped" ]] && continue

        unit="${UNITS[$mapped]}"

        if [[ -n "$unit" ]]; then
            CURRENT["$mapped"]="$raw_value $unit"
        else
            CURRENT["$mapped"]="$raw_value"
        fi

    done < <(
        awk '
            /command arg/ {

                value=$NF

                if (value ~ /^--/) {
                    key=value

                    getline
                    val=$NF

                    print key "|" val
                }
            }
        ' "$txtfile"
    )

    # ----------------------------------
    # Extract seed
    #
    # Example:
    # seed= 1350410594
    # ----------------------------------
    seed=$(
        grep "seed=" "$txtfile" |
        awk -F= '
            {
                gsub(/ /,"",$2)
                print $2
            }
        ' |
        head -n1
    )

    if [[ -n "${seed:-}" ]]; then
        CURRENT["micro_seed"]="$seed"
    fi

    # ----------------------------------
    # Save experiment data
    # ----------------------------------
    for param in "${!CURRENT[@]}"; do
        DATA["$folder_name|$param"]="${CURRENT[$param]}"
    done

done < <(
    find "$ROOT_DIR" \
        -type d -name "macro*" -prune -o \
        -type f -name "*.txt" -print | sort
)

# ----------------------------------
# Write CSV
# ----------------------------------
{

    # Header row
    printf "parameter"

    for folder in "${FOLDERS[@]}"; do
        printf ",%s" "$folder"
    done

    printf "\n"

    # Parameter rows
    for param in "${PARAM_ORDER[@]}"; do

        printf "%s" "$param"

        for folder in "${FOLDERS[@]}"; do
            printf ",%s" "${DATA[$folder|$param]:-}"
        done

        printf "\n"

    done

} > "$OUTPUT_CSV"

echo "Created: $OUTPUT_CSV"
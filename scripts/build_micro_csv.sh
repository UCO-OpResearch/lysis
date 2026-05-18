#!/usr/bin/env bash

set -euo pipefail

# Usage:
# ./build_micro_csv.sh SIM_DIRECTORY OUTPUT.csv
#
# Example:
# ./build_micro_csv.sh ./simulations micro_runs.csv

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
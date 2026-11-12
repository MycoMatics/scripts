#!/usr/bin/env bash
# Run as "bash run_analysis.sh" in tmux, never via "source run_analysis.sh".
# Each step is checked without changing interactive shell options.

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    printf 'ERROR: Run this script with: bash run_analysis.sh\n' >&2
    return 2
fi

if [[ -z "${PROJECT_ROOT:-}" ]]; then
    printf 'ERROR: Set PROJECT_ROOT to the UMIvsRAW project root.\n' >&2
    exit 2
fi
export REVIEW_DIR="${REVIEW_DIR:-$PROJECT_ROOT/REVIEW}"
THREADS="${THREADS:-63}"
WORKERS="${WORKERS:-20}"

cd "$(dirname "$0")" || exit 2

run_step() {
    local label="$1"
    shift
    printf '\n===== %s =====\n' "$label"
    "$@"
    local status=$?
    printf 'Finished %s: exit status %s\n' "$label" "$status"
    return "$status"
}

run_step '01 SINTAX classifications' python3 scripts/01_run_sintax.py --project-root "$PROJECT_ROOT" --threads "$THREADS" || exit $?
run_step '02 Build SINTAX tables' python3 scripts/02_build_sintax_tables.py || exit $?
run_step '03 Pair D1 and D3' python3 scripts/03_pair_d1_d3.py || exit $?
run_step '04 Build CD2' python3 scripts/04_build_cd2.py || exit $?
run_step '05 Four-dataset community analysis' Rscript scripts/05_community_analysis.R || exit $?
run_step '06 Classification success and A1 cohorts' python3 scripts/06_classification_success.py || exit $?
run_step '07 Chimeric mismatch audit: all 97' python3 scripts/07_mismatch_analysis.py \
  --review-dir "$REVIEW_DIR" \
  --cohort-file "$REVIEW_DIR/SINTAX_equal_ITS_v1/analysis/appendixC_current/cohorts/all97.txt" || exit $?
run_step '07 Chimeric mismatch audit: A1 90' python3 scripts/07_mismatch_analysis.py \
  --review-dir "$REVIEW_DIR" \
  --cohort-file "$REVIEW_DIR/SINTAX_equal_ITS_v1/analysis/appendixC_current/cohorts/A1filtered90.txt" || exit $?
run_step '08 SH-matching rate extraction' python3 scripts/08_extract_sh_matching_rates.py || exit $?
run_step '09 Cross-algorithm contingency' Rscript scripts/09_compare_annotation_methods.R || exit $?
run_step '10 Figure C3' Rscript scripts/10_plot_figure_C3.R || exit $?

IDENTITY_FILE="$REVIEW_DIR/raw_vs_umi_identity_harmonized/raw_vs_umi_consensus_identity.tsv"
if [[ ! -s "$IDENTITY_FILE" ]]; then
    run_step '11 Raw-consensus identity' python3 scripts/11_raw_vs_umi_identity.py --root "$PROJECT_ROOT" --workers "$WORKERS" || exit $?
else
    printf '\nSkipping identity calculation: existing input %s\n' "$IDENTITY_FILE"
fi

run_step '12 Figure B3' Rscript scripts/12_plot_figure_B3.R || exit $?
run_step '13 Appendix C audit summary' python3 scripts/13_compile_appendixC_summary.py || exit $?
printf '\nALL STEPS COMPLETE\n'

#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TCGA_DIR="$REPO_ROOT/data/raw/TCGA_PanCan"
GENE_LIST="$REPO_ROOT/docs/gene_lists/panel_B_pancan_drivers_53.txt"

STUDIES=(brca luad lusc coadread stad hnsc blca prad ov lihc)

VERIFY=0
[[ "${1:-}" == "--verify" ]] && VERIFY=1

if [[ ! -f "$GENE_LIST" ]]; then
  echo "ERROR: gene list not found: ${GENE_LIST#"$REPO_ROOT/"}" >&2
  exit 1
fi

n_genes=$(grep -v '^\s*#' "$GENE_LIST" | grep -c '[A-Za-z0-9]' || true)
echo "==========================================================================="
echo "  Panel B driver-subset builder$([[ $VERIFY -eq 1 ]] && echo '  [--verify]')"
echo "  Driver panel: $n_genes genes from ${GENE_LIST#"$REPO_ROOT/"}"
echo "==========================================================================="

emit_header() {
  cat <<'EOF'
EOF
}

build_one() {
  local full="$1"
  emit_header
  awk -F'\t' -v OFS='\t' '
    FNR == NR {                                    # pass 1: the gene list
      sub(/#.*/, "", $0); gsub(/[ \t\r]/, "", $0)
      if (length($0)) D[$0] = 1
      next
    }
    !hdr {                                         # pass 2: locate the header
      if ($1 == "Hugo_Symbol") {
        hdr = 1
        for (i = 1; i <= NF; i++) idx[$i] = i
        hg  = idx["Hugo_Symbol"]
        vc  = idx["Variant_Classification"]
        tsb = idx["Tumor_Sample_Barcode"]
        if (!hg || !vc || !tsb) {
          print "ERROR: MAF is missing a required column" > "/dev/stderr"; exit 2
        }
        split("Frame_Shift_Del Frame_Shift_Ins Splice_Site Translation_Start_Site \
               Nonsense_Mutation Nonstop_Mutation In_Frame_Del In_Frame_Ins \
               Missense_Mutation", ns, /[ \t\n]+/)
        for (i in ns) if (length(ns[i])) NS[ns[i]] = 1   # maftools vc_nonSyn default
        print
      }
      next
    }
    {
      keep = ($hg in D)
      if (!($tsb in seen_any)) { seen_any[$tsb] = 1; keep = 1 }
      if (($vc in NS) && !($tsb in seen_ns)) { seen_ns[$tsb] = 1; keep = 1 }
      if (keep) print
    }
    END { if (!hdr) { print "ERROR: no MAF header row found" > "/dev/stderr"; exit 2 } }
  ' "$GENE_LIST" "$full"
}

missing=0
rebuilt=0
differs=0

for cancer in "${STUDIES[@]}"; do
  study="${cancer}_tcga_pan_can_atlas_2018"
  full="$TCGA_DIR/$study/data_mutations.txt"
  out="$TCGA_DIR/$study/data_mutations.txt.gz"

  if [[ ! -s "$full" ]]; then
    echo "  skip (full MAF absent): $study"
    missing=$((missing + 1))
    continue
  fi

  if [[ $VERIFY -eq 1 ]]; then
    if [[ ! -s "$out" ]]; then
      echo "  FAIL  $study — no committed subset to verify against"
      differs=$((differs + 1))
      continue
    fi
    tmp="$(mktemp)"
    build_one "$full" > "$tmp"
    if gzip -cd "$out" | diff -q - "$tmp" >/dev/null; then
      echo "  ok    $study — committed subset matches a fresh rebuild"
    else
      echo "  FAIL  $study — committed subset differs from a fresh rebuild"
      differs=$((differs + 1))
    fi
    rm -f "$tmp"
  else
    build_one "$full" | gzip -9 -n > "$out"
    echo "  wrote $(printf '%7.2f MB' "$(echo "$(wc -c < "$out") / 1048576" | bc -l)")  ${out#"$TCGA_DIR/"}"
  fi
  rebuilt=$((rebuilt + 1))
done

echo ""
if [[ $missing -eq ${#STUDIES[@]} ]]; then
  echo "Nothing to do: no full study MAF is present."
  echo "This script is provenance only. Panel B reads the committed subsets and"
  echo "needs no download. To re-derive them, fetch the full MAFs first:"
  echo "    bash scripts/download_input_data.sh --full-tcga"
  exit 0
fi

if [[ $VERIFY -eq 1 ]]; then
  if [[ $differs -gt 0 ]]; then
    echo "VERIFY FAILED: $differs of $rebuilt subset(s) do not match a rebuild."
    exit 1
  fi
  echo "VERIFY OK: all $rebuilt rebuilt subset(s) match the committed files."
else
  echo "Built $rebuilt subset(s)."
  [[ $missing -gt 0 ]] && echo "($missing study MAF(s) were absent and were left alone.)"
  echo "Next: Rscript scripts/00_check_input_data.R, then re-run Panel B and"
  echo "confirm output/tables/Fig1A_*.csv are unchanged."
fi

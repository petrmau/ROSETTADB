#!/usr/bin/env bash
# download.sh — fetch / update all three AMR source databases
#
# Usage:
#   ./download.sh              # update all sources
#   ./download.sh card         # update CARD only
#   ./download.sh resfinder    # update ResFinder only
#   ./download.sh amrfinder    # update AMRFinder only
#
# Re-entrant / idempotent: safe to run repeatedly.
# Each run writes a VERSION file recording the download timestamp.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCES="${SCRIPT_DIR}/sources"

# ── colours ──────────────────────────────────────────────────────────────────
BOLD="\033[1m"; GREEN="\033[32m"; YELLOW="\033[33m"; RESET="\033[0m"
info()  { echo -e "${BOLD}[INFO]${RESET} $*"; }
ok()    { echo -e "${GREEN}[OK]${RESET}   $*"; }
warn()  { echo -e "${YELLOW}[WARN]${RESET} $*"; }

# ── CARD ─────────────────────────────────────────────────────────────────────
download_card() {
    info "Downloading CARD broadstreet data …"
    local dest="${SOURCES}/CARD"
    mkdir -p "${dest}"

    local url="https://card.mcmaster.ca/latest/data"
    local archive="${dest}/card-data.tar.bz2"

    wget --quiet --show-progress \
         --no-use-server-timestamps \
         -O "${archive}" \
         "${url}"

    info "Extracting CARD archive …"
    # Extract into dest; -C changes directory, overwriting existing files.
    tar xjf "${archive}" -C "${dest}"

    # Record the ARO version (present inside the tarball as card.json "dataVersion")
    local ver=""
    if command -v python3 &>/dev/null && [[ -f "${dest}/card.json" ]]; then
        ver=$(python3 -c "
import json, sys
with open('${dest}/card.json') as f:
    d = json.load(f)
print(d.get('_comment', {}).get('dataVersion', '') or
      d.get('_version', ''))
" 2>/dev/null || true)
    fi
    {
        echo "downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "source_url: ${url}"
        [[ -n "${ver}" ]] && echo "data_version: ${ver}"
    } > "${dest}/VERSION"

    ok "CARD done  (${dest}/VERSION)"
}

# ── ResFinder ─────────────────────────────────────────────────────────────────
download_resfinder() {
    info "Updating ResFinder database …"
    local dest="${SOURCES}/resfinder_db"

    if [[ -d "${dest}/.git" ]]; then
        info "  git pull …"
        git -C "${dest}" pull --ff-only --quiet
    else
        info "  git clone …"
        # Remove any partial directory before cloning
        rm -rf "${dest}"
        git clone --quiet \
            https://bitbucket.org/genomicepidemiology/resfinder_db.git \
            "${dest}"
    fi

    local commit
    commit=$(git -C "${dest}" rev-parse --short HEAD)
    {
        echo "downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "git_commit: ${commit}"
        echo "repo: https://bitbucket.org/genomicepidemiology/resfinder_db"
        [[ -f "${dest}/VERSION" ]] && echo "db_version: $(cat "${dest}/VERSION")"
    } > "${dest}/ROSETTA_VERSION"

    ok "ResFinder done  (commit ${commit})"
}

# ── AMRFinder ─────────────────────────────────────────────────────────────────
download_amrfinder() {
    info "Downloading AMRFinder Plus database …"
    local dest="${SOURCES}/amr_finder_plus"
    mkdir -p "${dest}"

    local base="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"

    # Files required by ingest.py + harmonise.py
    local files=(
        "AMR_CDS.fa"
        "ReferenceGeneCatalog.txt"
        "AMRProt.fa"
        "fam.tsv"
        "md5sum.txt"
    )

    for fname in "${files[@]}"; do
        info "  downloading ${fname} …"
        wget --quiet --show-progress \
             --no-use-server-timestamps \
             -O "${dest}/${fname}" \
             "${base}/${fname}"
    done

    # Verify checksums (md5sum.txt covers all files in the release)
    info "Verifying checksums …"
    local ok_count=0
    local fail_count=0
    while IFS= read -r line; do
        expected_hash="${line%% *}"
        fname="${line##* }"              # last space-separated token
        fname="${fname##*/}"             # basename only
        local fpath="${dest}/${fname}"
        if [[ ! -f "${fpath}" ]]; then
            continue   # skip files we didn't download
        fi
        actual_hash=$(md5sum "${fpath}" | awk '{print $1}')
        if [[ "${actual_hash}" == "${expected_hash}" ]]; then
            ok_count=$((ok_count + 1))
        else
            warn "Checksum MISMATCH for ${fname}"
            fail_count=$((fail_count + 1))
        fi
    done < "${dest}/md5sum.txt"

    if [[ ${fail_count} -gt 0 ]]; then
        echo "ERROR: ${fail_count} checksum failure(s). Aborting." >&2
        exit 1
    fi
    info "  ${ok_count} checksum(s) verified."

    # Read database version from md5sum.txt parent directory listing via latestdir symlink
    # The FTP index does not expose the version directly in these files; use
    # the build date embedded in ReferenceGeneCatalog.txt header comment if present,
    # or fall back to the parent-directory name resolved via HTTP redirect.
    local db_ver=""
    db_ver=$(wget -qO- "${base}/../" 2>/dev/null \
        | grep -oP '(?<=href=")[0-9]{4}-[0-9]{2}-[0-9]{2}\.[0-9]+(?=/)' \
        | sort -V | tail -1 || true)
    {
        echo "downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "source_url: ${base}"
        [[ -n "${db_ver}" ]] && echo "db_version: ${db_ver}"
    } > "${dest}/VERSION"

    ok "AMRFinder done  (${dest}/VERSION)"
}

# ── Dispatch ──────────────────────────────────────────────────────────────────
TARGET="${1:-all}"

case "${TARGET}" in
    all)
        download_card
        download_resfinder
        download_amrfinder
        ;;
    card)       download_card ;;
    resfinder)  download_resfinder ;;
    amrfinder)  download_amrfinder ;;
    *)
        echo "Usage: $0 [all|card|resfinder|amrfinder]" >&2
        exit 1
        ;;
esac

echo ""
ok "All requested sources up to date."

#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
dockerfile="$script_dir/Dockerfile"

docker_image=${GENEINIT_PREDICTOR_IMAGE:-geneinit-predictor:ubuntu24.04-etp1.02-augustus3.5.0}
output_dir=${GENEINIT_OUTPUT_DIR:-$repo_root/dist}
sif_name=${GENEINIT_SIF_NAME:-geneinit_predictor_ubuntu24.04_etp1.02_augustus3.5.0.sif}
source_cache_dir=${GENEINIT_SOURCE_CACHE_DIR:-$output_dir/sources}
augustus_archive=${AUGUSTUS_ARCHIVE:-$source_cache_dir/Augustus-v3.5.0.tar.gz}
augustus_url=${AUGUSTUS_URL:-https://github.com/Gaius-Augustus/Augustus/archive/refs/tags/v3.5.0.tar.gz}
genemark_archive=${GENEMARK_ETP_ARCHIVE:-}
genemark_key=${GENEMARK_KEY_PATH:-}
standalone_archive=${GENEMARK_STANDALONE_ARCHIVE:-}

if [[ "$sif_name" == */* || "$sif_name" != *.sif ]]; then
    echo "GENEINIT_SIF_NAME must be a .sif basename" >&2
    exit 2
fi
if [[ -z "$genemark_archive" || -z "$genemark_key" || -z "$standalone_archive" ]]; then
    echo "GENEMARK_ETP_ARCHIVE, GENEMARK_KEY_PATH, and GENEMARK_STANDALONE_ARCHIVE are required" >&2
    exit 2
fi

sif="$output_dir/$sif_name"
sif_partial="$sif.partial"
version_report="${sif%.sif}.versions.txt"
source_stage=
smoke_dir=
smoke_home=

cleanup() {
    [[ -z "$source_stage" ]] || rm -rf -- "$source_stage"
    [[ -z "$smoke_dir" ]] || rm -rf -- "$smoke_dir"
    [[ -z "$smoke_home" ]] || rm -rf -- "$smoke_home"
    rm -f -- "$sif_partial"
}
trap cleanup EXIT

mkdir -p "$output_dir" "$source_cache_dir"
test -s "$dockerfile"
test -s "$genemark_archive"
test -s "$genemark_key"
test -s "$standalone_archive"
test ! -e "$sif"
test ! -e "$sif_partial"
command -v docker >/dev/null
command -v apptainer >/dev/null

if [[ ! -s "$augustus_archive" ]]; then
    command -v curl >/dev/null
    augustus_partial="$augustus_archive.partial"
    rm -f -- "$augustus_partial"
    curl -L --fail --retry 3 --output "$augustus_partial" "$augustus_url"
    tar -tzf "$augustus_partial" >/dev/null
    mv "$augustus_partial" "$augustus_archive"
fi
tar -tzf "$augustus_archive" >/dev/null
tar -xOzf "$augustus_archive" Augustus-3.5.0/include/types.hh \
    | grep -Fx '#define VERSION "3.5.0"'

source_stage=$(mktemp -d "$output_dir/predictor_source_stage.XXXXXX")
install -m 644 "$genemark_archive" "$source_stage/GeneMark-ETP.tar.gz"
install -m 644 "$augustus_archive" "$source_stage/Augustus-v3.5.0.tar.gz"

DOCKER_BUILDKIT=1 docker build \
    --build-context "sources=$source_stage" \
    --file "$dockerfile" \
    --tag "$docker_image" \
    "$script_dir"

apptainer build "$sif_partial" "docker-daemon://$docker_image"
mv "$sif_partial" "$sif"

apptainer exec --cleanenv "$sif" /bin/bash -lc \
    "cat /etc/os-release; /opt/Augustus/bin/augustus --version; grep -m1 'my .version' /opt/ETP/bin/gmetp.pl" \
    > "$version_report" 2>&1
grep -F 'VERSION_ID="24.04"' "$version_report"
grep -F '3.5.0' "$version_report"
grep -F '1.02' "$version_report"

smoke_dir=$(mktemp -d "$output_dir/gmcore_smoke.XXXXXX")
smoke_home=$(mktemp -d "$output_dir/gmcore_home.XXXXXX")
install -m 600 "$genemark_key" "$smoke_home/.gm_key"
tar -xzf "$standalone_archive" -C "$smoke_dir" \
    gmes_linux_64_4/GeneMark-E-tests/GeneMark.hmm/input/athaliana.mod \
    gmes_linux_64_4/GeneMark-E-tests/GeneMark.hmm/input/sequence.fasta

model="$smoke_dir/gmes_linux_64_4/GeneMark-E-tests/GeneMark.hmm/input/athaliana.mod"
sequence="$smoke_dir/gmes_linux_64_4/GeneMark-E-tests/GeneMark.hmm/input/sequence.fasta"
prediction="$smoke_dir/genemark.gff3"
apptainer exec --cleanenv --home "$smoke_home" --bind "$smoke_dir:$smoke_dir" "$sif" \
    /opt/ETP/bin/gmes/gmhmme3_original \
    -o "$prediction" -m "$model" -f gff3 "$sequence"
test -s "$prediction"

sha256sum "$sif" > "$sif.sha256"
printf 'predictor_sif\t%s\n' "$sif"
printf 'version_report\t%s\n' "$version_report"
printf 'licensed_core_smoke\tPASS\n'

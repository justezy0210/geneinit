#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
dockerfile="$script_dir/Dockerfile"

predictor_image=${GENEINIT_PREDICTOR_IMAGE:-geneinit-predictor:ubuntu24.04-etp1.02-augustus3.5.0}
build_revision=${GENEINIT_BUILD_REVISION:-}
default_image=geneinit-allinone:ubuntu24.04-etp1.02-augustus3.5.0
if [[ -n "$build_revision" ]]; then
    default_image="${default_image}-${build_revision}"
fi
docker_image=${GENEINIT_ALLINONE_IMAGE:-$default_image}
output_dir=${GENEINIT_OUTPUT_DIR:-$repo_root/dist}
sif_name=${GENEINIT_SIF_NAME:-geneinit_allinone_ubuntu24.04_etp1.02_augustus3.5.0.sif}
genemark_key=${GENEMARK_KEY_PATH:-}
standalone_archive=${GENEMARK_STANDALONE_ARCHIVE:-}

if [[ "$sif_name" == */* || "$sif_name" != *.sif ]]; then
    echo "GENEINIT_SIF_NAME must be a .sif basename" >&2
    exit 2
fi
if [[ -z "$genemark_key" || -z "$standalone_archive" ]]; then
    echo "GENEMARK_KEY_PATH and GENEMARK_STANDALONE_ARCHIVE are required" >&2
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

mkdir -p "$output_dir"
test -s "$dockerfile"
test -s "$genemark_key"
test -s "$standalone_archive"
test -s "$repo_root/Snakefile"
test -s "$repo_root/profiles/allinone/config.yaml"
test ! -e "$sif"
test ! -e "$sif_partial"
command -v docker >/dev/null
command -v apptainer >/dev/null

source_stage=$(mktemp -d "$output_dir/geneinit_source_stage.XXXXXX")
for item in Snakefile geneinit workflow config policies profiles envs containers README.md pyproject.toml; do
    cp -a -- "$repo_root/$item" "$source_stage/"
done

DOCKER_BUILDKIT=1 docker build \
    --build-arg "PREDICTOR_IMAGE=$predictor_image" \
    --build-context "geneinit_src=$source_stage" \
    --file "$dockerfile" \
    --tag "$docker_image" \
    "$script_dir"

apptainer build "$sif_partial" "docker-daemon://$docker_image"
mv "$sif_partial" "$sif"

apptainer exec --cleanenv "$sif" geneinit-software-versions \
    > "$version_report" 2>&1
grep -F 'VERSION_ID="24.04"' "$version_report"
grep -F $'snakemake\t9.23.1' "$version_report"
grep -F $'ingenannot\t0.0.20' "$version_report"
grep -F '3.5.0' "$version_report"
grep -F '1.02' "$version_report"

apptainer exec --cleanenv "$sif" /bin/bash -o pipefail -c '
set -euo pipefail
if grep -RIl "^#!/usr/bin/env perl$" /opt/ETP/bin; then
    echo "GeneMark-ETP still contains env-based Perl shebangs" >&2
    exit 1
fi
for module in Hash::Merge MCE::Mutex Math::Utils Parallel::ForkManager YAML YAML::XS; do
    /usr/bin/perl -M"${module}" -e 1
done
while IFS= read -r -d "" script; do
    PERL5LIB=/opt/ETP/bin/gmes/lib /usr/bin/perl -c "${script}"
done < <(find /opt/ETP/bin -type f -name "*.pl" -print0)
'

smoke_dir=$(mktemp -d "$output_dir/gmcore_allinone_smoke.XXXXXX")
smoke_home=$(mktemp -d "$output_dir/gmcore_allinone_home.XXXXXX")
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
printf 'allinone_sif\t%s\n' "$sif"
printf 'version_report\t%s\n' "$version_report"
printf 'licensed_core_smoke\tPASS\n'
printf 'etp_perl_smoke\tPASS\n'

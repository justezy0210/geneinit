# syntax=docker/dockerfile:1
ARG PREDICTOR_IMAGE=geneinit-predictor:ubuntu24.04-etp1.02-augustus3.5.0

FROM mambaorg/micromamba:2.9.0 AS micromamba

FROM ${PREDICTOR_IMAGE}

LABEL org.opencontainers.image.title="GeneInIt all-in-one stack" \
      org.opencontainers.image.description="Ubuntu 24.04 GeneInIt runtime with RNA evidence, GeneMark-ETP, AUGUSTUS, and InGenAnnot" \
      org.geneinit.execution-mode="in_container" \
      org.geneinit.genemark-etp.version="1.02" \
      org.geneinit.augustus.version="3.5.0" \
      org.geneinit.snakemake.version="9.23.1" \
      org.geneinit.ingenannot.version="0.0.20" \
      org.geneinit.license-key.embedded="false"

USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

COPY --from=micromamba /bin/micromamba /usr/local/bin/micromamba
COPY environment-main.yaml /tmp/environment-main.yaml
COPY environment-alignment.yaml /tmp/environment-alignment.yaml
COPY environment-stringtie.yaml /tmp/environment-stringtie.yaml
COPY environment-utr.yaml /tmp/environment-utr.yaml

ENV MAMBA_ROOT_PREFIX=/opt/micromamba

RUN set -eux; \
    micromamba create --yes --prefix /opt/geneinit-main \
        --file /tmp/environment-main.yaml; \
    micromamba create --yes --prefix /opt/geneinit-alignment \
        --file /tmp/environment-alignment.yaml; \
    micromamba create --yes --prefix /opt/geneinit-stringtie \
        --file /tmp/environment-stringtie.yaml; \
    micromamba create --yes --prefix /opt/geneinit-utr \
        --file /tmp/environment-utr.yaml; \
    micromamba clean --all --yes; \
    rm -rf /opt/micromamba/pkgs /root/.cache/pip \
        /tmp/environment-main.yaml /tmp/environment-alignment.yaml \
        /tmp/environment-stringtie.yaml /tmp/environment-utr.yaml

COPY --from=geneinit_src . /opt/GeneInIt
RUN /opt/geneinit-utr/bin/python /opt/GeneInIt/containers/allinone/patch_ingenannot.py
COPY geneinit-run /usr/local/bin/geneinit-run
COPY geneinit-software-versions /usr/local/bin/geneinit-software-versions

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONNOUSERSITE=1 \
    PYTHONPATH=/opt/GeneInIt \
    PATH=/opt/geneinit-main/bin:/opt/geneinit-alignment/bin:/opt/geneinit-stringtie/bin:/opt/geneinit-utr/bin:/opt/ETP/bin:/opt/ETP/tools:/opt/ETP/bin/gmes:/opt/ETP/bin/gmes/ProtHint/bin:/opt/ETP/bin/gmes/ProtHint/dependencies:/opt/Augustus/bin:/opt/Augustus/scripts:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Conda environments precede /usr/bin in PATH and may contain an incompatible
# Perl. Pin every GeneMark-ETP Perl entry point to the system interpreter whose
# CPAN modules are installed in the predictor base image.
RUN set -eux; \
    find /opt/ETP/bin -type f -name '*.pl' \
        -exec sed -i '1s|^#!/usr/bin/env perl$|#!/usr/bin/perl|' {} +; \
    ! grep -RIl '^#!/usr/bin/env perl$' /opt/ETP/bin; \
    for module in Hash::Merge MCE::Mutex Math::Utils Parallel::ForkManager YAML YAML::XS; do \
        /usr/bin/perl -M"${module}" -e 1; \
    done; \
    while IFS= read -r -d '' script; do \
        PERL5LIB=/opt/ETP/bin/gmes/lib /usr/bin/perl -c "${script}"; \
    done < <(find /opt/ETP/bin -type f -name '*.pl' -print0)

RUN set -eux; \
    chmod 0755 /usr/local/bin/geneinit-run /usr/local/bin/geneinit-software-versions; \
    test -s /opt/GeneInIt/Snakefile; \
    test -s /opt/GeneInIt/profiles/allinone/config.yaml; \
    cd /tmp; \
    python -c 'import geneinit'; \
    ingenannot --help >/dev/null; \
    geneinit-software-versions > /opt/GeneInIt/software_versions.build.txt; \
    grep -F 'VERSION_ID="24.04"' /opt/GeneInIt/software_versions.build.txt; \
    grep -F $'snakemake\t9.23.1' /opt/GeneInIt/software_versions.build.txt; \
    grep -F $'ingenannot\t0.0.20' /opt/GeneInIt/software_versions.build.txt; \
    grep -F '3.5.0' /opt/GeneInIt/software_versions.build.txt; \
    grep -F '1.02' /opt/GeneInIt/software_versions.build.txt

WORKDIR /work

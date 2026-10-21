# syntax=docker/dockerfile:1
FROM ubuntu:24.04 AS builder

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

COPY --from=sources GeneMark-ETP.tar.gz /opt/src/GeneMark-ETP.tar.gz
COPY --from=sources Augustus-v3.5.0.tar.gz /opt/src/Augustus-v3.5.0.tar.gz

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        build-essential \
        libboost-iostreams-dev \
        libhash-merge-perl \
        libmath-utils-perl \
        libmce-perl \
        libparallel-forkmanager-perl \
        libthread-queue-perl \
        libyaml-libyaml-perl \
        libyaml-perl \
        perl \
        python3 \
        zlib1g-dev; \
    mkdir -p /opt/ETP /opt/Augustus; \
    tar -xzf /opt/src/GeneMark-ETP.tar.gz -C /opt/ETP --strip-components=1; \
    tar -xzf /opt/src/Augustus-v3.5.0.tar.gz -C /opt/Augustus --strip-components=1; \
    grep -Fx '#define VERSION "3.5.0"' /opt/Augustus/include/types.hh; \
    sed -i 's/^COMPGENEPRED = true$/COMPGENEPRED = false/' /opt/Augustus/common.mk; \
    grep -Fx 'COMPGENEPRED = false' /opt/Augustus/common.mk; \
    make -C /opt/Augustus -j"$(nproc)" augustus; \
    find /opt/ETP/bin -type f -name '*.pl' \
        -exec sed -i '1s|^#!/usr/bin/env perl$|#!/usr/bin/perl|' {} +; \
    ! grep -RIl '^#!/usr/bin/env perl$' /opt/ETP/bin; \
    find /opt/ETP/bin -type f \( -name '*.pl' -o -name '*.py' -o -name '*.sh' \) \
        -exec chmod a+rx {} +; \
    find /opt/ETP/tools -maxdepth 1 -type f -exec chmod a+rx {} +; \
    chmod a+rx \
        /opt/ETP/bin/bam2hints \
        /opt/ETP/bin/probuild \
        /opt/ETP/bin/gmes/gmhmme3_original \
        /opt/ETP/bin/gmes/probuild \
        /opt/ETP/bin/gmst/gmhmmp \
        /opt/ETP/bin/gmst/probuild; \
    chmod a+rx /opt/Augustus/bin/* /opt/Augustus/scripts/*.pl; \
    PERL5LIB=/opt/ETP/bin/gmes/lib /usr/bin/perl -MYAML::XS -MMCE::Mutex \
        -MThread::Queue -MMath::Utils -MHash::Merge \
        -MParallel::ForkManager -MLineFit -e 1; \
    while IFS= read -r -d '' script; do \
        PERL5LIB=/opt/ETP/bin/gmes/lib /usr/bin/perl -c "${script}"; \
    done < <(find /opt/ETP/bin -type f -name '*.pl' -print0); \
    test -s /opt/ETP/bin/gmetp.pl; \
    test -s /opt/ETP/bin/format_back.pl; \
    test -x /opt/Augustus/bin/augustus; \
    test -x /opt/Augustus/bin/etraining; \
    test -s /opt/Augustus/scripts/new_species.pl; \
    test -s /opt/Augustus/scripts/gff2gbSmallDNA.pl; \
    test -s /opt/Augustus/scripts/optimize_augustus.pl; \
    test -s /opt/Augustus/scripts/join_mult_hints.pl; \
    test -s /opt/Augustus/config/species/generic/generic_parameters.cfg; \
    /opt/Augustus/bin/augustus --version 2>&1 | grep -F '3.5.0'; \
    grep -F 'my $version = "1.02";' /opt/ETP/bin/gmetp.pl

FROM ubuntu:24.04

LABEL org.opencontainers.image.title="GeneInIt predictor stack" \
      org.opencontainers.image.description="Clean Ubuntu 24.04 runtime with GeneMark-ETP 1.02 and AUGUSTUS 3.5.0" \
      org.opencontainers.image.base.name="docker.io/library/ubuntu:24.04" \
      org.opencontainers.image.source.augustus="https://github.com/Gaius-Augustus/Augustus/releases/tag/v3.5.0" \
      org.geneinit.genemark-etp.version="1.02" \
      org.geneinit.augustus.version="3.5.0" \
      org.geneinit.augustus.compgenepred="false" \
      org.geneinit.license-key.embedded="false"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        coreutils \
        findutils \
        gawk \
        grep \
        gzip \
        libboost-iostreams1.83.0 \
        libhash-merge-perl \
        libmath-utils-perl \
        libmce-perl \
        libparallel-forkmanager-perl \
        libthread-queue-perl \
        libyaml-libyaml-perl \
        libyaml-perl \
        perl \
        python3 \
        sed \
        tar; \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/ETP /opt/ETP
COPY --from=builder /opt/Augustus/bin /opt/Augustus/bin
COPY --from=builder /opt/Augustus/config /opt/Augustus/config
COPY --from=builder /opt/Augustus/scripts /opt/Augustus/scripts

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    GENEMARK_PATH=/opt/ETP/bin \
    AUGUSTUS_BIN_PATH=/opt/Augustus/bin \
    AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts \
    AUGUSTUS_CONFIG_PATH=/opt/Augustus/config \
    PERL5LIB=/opt/ETP/bin/gmes/lib \
    PATH=/opt/ETP/bin:/opt/ETP/tools:/opt/ETP/bin/gmes:/opt/ETP/bin/gmes/ProtHint/bin:/opt/ETP/bin/gmes/ProtHint/dependencies:/opt/Augustus/bin:/opt/Augustus/scripts:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN set -eux; \
    /usr/bin/perl -MYAML::XS -MMCE::Mutex -MThread::Queue -MMath::Utils \
        -MHash::Merge -MParallel::ForkManager -MLineFit -e 1; \
    ! grep -RIl '^#!/usr/bin/env perl$' /opt/ETP/bin; \
    while IFS= read -r -d '' script; do \
        /usr/bin/perl -c "${script}"; \
    done < <(find /opt/ETP/bin -type f -name '*.pl' -print0); \
    ! ldd /opt/Augustus/bin/augustus | grep -F 'not found'; \
    /opt/ETP/tools/bedtools --version; \
    /opt/ETP/tools/diamond version; \
    /opt/ETP/tools/samtools --version | head -n 1; \
    /opt/Augustus/bin/augustus --version 2>&1 | grep -F '3.5.0'; \
    grep -F 'my $version = "1.02";' /opt/ETP/bin/gmetp.pl; \
    test -s /opt/Augustus/scripts/new_species.pl; \
    test -s /opt/Augustus/scripts/gff2gbSmallDNA.pl; \
    test -s /opt/Augustus/scripts/optimize_augustus.pl; \
    test -s /opt/Augustus/scripts/join_mult_hints.pl; \
    test -s /opt/Augustus/config/species/generic/generic_parameters.cfg

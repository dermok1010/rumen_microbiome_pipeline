# git, java 17, wget, curl all already ship in this base image's conda env
# (confirmed 2026-09-14); unzip isn't used anywhere in this pipeline. The
# apt-get layer that used to install these was also broken -- bullseye's
# live security mirror no longer serves several of the exact package
# builds it pinned (openjdk-17-jre-headless, unzip, curl, libtiff5, ...
# 404 on deb.debian.org), since this base image is old enough now that
# those point releases have aged off the live mirror. Removing the layer
# fixes the build and removes an unnecessary dependency on that mirror
# staying available at all.
FROM quay.io/qiime2/amplicon:2024.2

WORKDIR /pipeline

CMD ["/bin/bash"]

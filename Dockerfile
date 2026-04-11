FROM quay.io/qiime2/amplicon:2024.2

RUN apt-get update && apt-get install -y \
    git \
    openjdk-17-jre \
    wget \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /pipeline

CMD ["/bin/bash"]

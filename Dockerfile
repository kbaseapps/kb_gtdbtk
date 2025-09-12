FROM kbase/sdkpython:3.8.10
LABEL maintainer="Dylan Chivian"

# -------------------------------
# System deps used by the toolchain & installers
# -------------------------------
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        libgomp1 unzip curl git xvfb ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Make sure conda/python is first on PATH
ENV PATH="/opt/conda3/bin:${PATH}"

# -------------------------------
# Python tooling & libraries used in this repo
# -------------------------------
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
        pytest pytest-cov mypy coveralls flake8 \
        jsonrpcbase requests pandas

# Prodigal barfs on NumPy >= 1.24 due to np.bool removal; keep 1.23.x
RUN pip install --no-cache-dir "numpy==1.23.1"

# -------------------------------
# GTDB-Tk 2.4.1
# -------------------------------
ENV GTDBTK_VERSION="2.4.1"
RUN pip install --no-cache-dir "gtdbtk==${GTDBTK_VERSION}"

# Where the GTDB database will be mounted/staged
ENV GTDBTK_DATA_PATH="/data/gtdbtk"

# -------------------------------
# GTDB-Tk runtime dependencies
#   - FastANI
#   - hmmer, prodigal, pplacer, fasttree
#   - Mash 2.3
# -------------------------------
# FastANI v1.33
ENV FASTANI_VERSION="v1.33"
RUN curl -LJO "https://github.com/ParBLiSS/FastANI/releases/download/${FASTANI_VERSION}/fastANI-Linux64-${FASTANI_VERSION}.zip" && \
    unzip -q "fastANI-Linux64-${FASTANI_VERSION}.zip" && \
    mv fastANI /usr/local/bin/ && \
    rm -f "fastANI-Linux64-${FASTANI_VERSION}.zip"

# HMMER / Prodigal / pplacer / FastTree from bioconda
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -y -c bioconda hmmer prodigal pplacer fasttree && \
    conda clean -afy

# Mash v2.3 (binary tarball)
ENV MASH_VERSION="v2.3"
RUN curl -LJO "https://github.com/marbl/Mash/releases/download/${MASH_VERSION}/mash-Linux64-${MASH_VERSION}.tar" && \
    tar xf "mash-Linux64-${MASH_VERSION}.tar" && \
    mv "mash-Linux64-${MASH_VERSION}/mash" /usr/local/bin/ && \
    rm -rf "mash-Linux64-${MASH_VERSION}" "mash-Linux64-${MASH_VERSION}.tar"

# -------------------------------
# Krona
# -------------------------------
WORKDIR /kb
RUN git clone https://github.com/marbl/Krona.git && \
    cd Krona/KronaTools && \
    ./install.pl

# -------------------------------
# ETE3 + PyQt5 (headless via xvfb when needed)
# Note: original pin kept
# -------------------------------
RUN pip install --no-cache-dir "ete3==3.1.3" "PyQt5==5.11.3" || true

# -------------------------------
# tree_utils (as before)
# -------------------------------
WORKDIR /kb/module
RUN git clone https://github.com/dcchivian/tree_utils && \
    mkdir -p bin && \
    mv tree_utils/gtdb/* bin/ && \
    chmod +x bin/*

# -------------------------------
# Module wiring
# -------------------------------
COPY ./ /kb/module
RUN mkdir -p /kb/module/work && \
    chmod -R a+rw /kb/module

WORKDIR /kb/module
RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]

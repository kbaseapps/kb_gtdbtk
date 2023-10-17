FROM kbase/sdkpython:3.8.10
MAINTAINER Dylan Chivian
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update -y
RUN apt-get install -y libgomp1 unzip

RUN pip install pip --upgrade
RUN pip install pytest pytest-cov mypy coveralls flake8 --upgrade \
    && pip install jsonrpcbase requests pandas --upgrade

# Prodigal barfing on newer numpy from np.bool
RUN pip install --upgrade numpy==1.23.1

# GTDB-Tk install
ENV GTDBTK_VERSION='2.3.2'
RUN pip install gtdbtk==${GTDBTK_VERSION}

# GTDB-Tk dependencies: FastANI
ENV FASTANI_VERSION='v1.33'
RUN curl -LJO https://github.com/ParBLiSS/FastANI/releases/download/${FASTANI_VERSION}/fastANI-Linux64-${FASTANI_VERSION}.zip \
    && unzip fastANI-Linux64-${FASTANI_VERSION}.zip \
    && mv fastANI /usr/local/bin/

# GTDB-Tk dependencies: hmmer, prodigal, pplacer, fasttree
RUN conda install -c bioconda hmmer prodigal pplacer fasttree --yes

# GTDB-Tk dependencies: update mash to v2.3 (from v1.1)
ENV MASH_VERSION='v2.3'
#RUN conda install -c bioconda "mash>=2.3" --yes
RUN curl -LJO https://github.com/marbl/Mash/releases/download/${MASH_VERSION}/mash-Linux64-${MASH_VERSION}.tar \
    && tar xf mash-Linux64-${MASH_VERSION}.tar \
    && mv mash-Linux64-${MASH_VERSION}/mash /usr/local/bin \
    && rm mash-Linux64-${MASH_VERSION}.tar

# Krona install
WORKDIR /kb
RUN git clone https://github.com/marbl/Krona.git && \
    cd Krona/KronaTools/ && \
    ./install.pl

# Install ETE3
RUN apt-get -y install xvfb
# Note: You must use PyQt5==5.11.3 on debian
RUN pip install ete3==3.1.3 PyQt5==5.11.3

# tree utils install
WORKDIR /kb/module
RUN git clone https://github.com/dcchivian/tree_utils && \
    mkdir bin && \
    mv tree_utils/gtdb/* bin/ && \
    chmod +x bin/*

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

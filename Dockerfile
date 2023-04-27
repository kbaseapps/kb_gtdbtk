FROM kbase/sdkpython:3.8.10
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update
RUN apt-get install -y libgomp1 unzip

RUN pip install pip --upgrade
RUN pip install pytest pytest-cov mypy coveralls flake8 --upgrade \
    && pip install jsonrpcbase requests pandas --upgrade

# Prodigal barfing on newer numpy from np.bool
RUN pip uninstall -y numpy
RUN yes | pip install numpy==1.23.1

# GTDB-Tk install
ENV GTDBTK_VERSION='2.2.1'
RUN pip install gtdbtk==${GTDBTK_VERSION}

# GTDB-Tk dependencies
ENV FASTANI_VERSION='v1.33'
RUN curl -LJO https://github.com/ParBLiSS/FastANI/releases/download/${FASTANI_VERSION}/fastANI-Linux64-${FASTANI_VERSION}.zip \
    && unzip fastANI-Linux64-${FASTANI_VERSION}.zip \
    && mv fastANI /usr/local/bin/

RUN conda install -c bioconda hmmer prodigal pplacer fasttree mash --yes

# Krona install
WORKDIR /kb
RUN git clone https://github.com/marbl/Krona.git && \
    cd Krona/KronaTools/ && \
    ./install.pl
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENV GTDBTK_DATA_PATH=/data
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

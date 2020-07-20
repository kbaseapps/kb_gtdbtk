FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update
RUN apt-get install libgomp1  

RUN curl -LJO https://github.com/ParBLiSS/FastANI/releases/download/v1.31/fastANI-Linux64-v1.31.zip \
&& unzip fastANI-Linux64-v1.31.zip \
&& mv fastANI /miniconda/bin/

RUN pip install pipenv==2018.11.26

ENV GTDBTK_DATA_PATH=/data
# conda updates to py 3.8 and everything breaks
RUN conda install -c bioconda hmmer prodigal pplacer fasttree --yes
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

# This stuff has to come after the COPY since it uses the pipfile in the repo
# really need a test build and a prod build. Not sure that's possible via sdk.
RUN pipenv install --system --deploy --ignore-pipfile --dev

RUN make all


ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update
RUN apt-get install libgomp1  
RUN pip install --upgrade pip

# update security libraries in the base image (deprecated approach)
# Does this update for python3 or just 2.7?  what about conda python?
RUN pip install cffi --upgrade \
    && pip install pyOpenSSL --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

RUN pip install pytest --upgrade \
    && pip install pytest-cov --upgrade \
    && pip install mypy --upgrade \
    && pip install coveralls --upgrade \
    && pip install flake8 --upgrade

RUN pip install jsonrpcbase --upgrade \
    && pip install requests --upgrade \
    && pip install pandas --upgrade \
    && pip install gtdbtk --upgrade

# update security libraries in the base image (deprecated approach)
# Does this update for python3 or just 2.7?  what about conda python?
#RUN pip install cffi --upgrade \
#    && pip install pyopenssl --upgrade \
#    && pip install ndg-httpsclient --upgrade \
#    && pip install pyasn1 --upgrade \
#    && pip install requests --upgrade \
#    && pip install 'requests[security]' --upgrade


# GTDB-Tk specific installs
ENV FASTANI_VERSION='v1.33'

RUN curl -LJO https://github.com/ParBLiSS/FastANI/releases/download/${FASTANI_VERSION}/fastANI-Linux64-${FASTANI_VERSION}.zip \
&& unzip fastANI-Linux64-${FASTANI_VERSION}.zip \
&& mv fastANI /miniconda/bin/

RUN pip install --upgrade pipenv
#RUN pip install pipenv==2018.11.26

ENV GTDBTK_DATA_PATH=/data
# conda updates to py 3.8 and everything breaks
RUN conda install -c bioconda hmmer prodigal pplacer fasttree mash --yes
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

# This stuff has to come after the COPY since it uses the pipfile in the repo
# really need a test build and a prod build. Not sure that's possible via sdk.
#RUN pipenv install --system --deploy --ignore-pipfile --dev

RUN make all


ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

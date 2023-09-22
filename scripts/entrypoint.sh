#!/bin/bash

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  #sh ./scripts/run_async.sh
  xvfb-run bash ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"

  # r207 refdata
  cd /data
  echo CWD `pwd`
  
  if [ ! -d "/data/r207" ] ; then
    export GTDB_VER_INT=207
    export GTDB_VER_FLT=207.0

    mkdir r${GTDB_VER_INT}
    cd r${GTDB_VER_INT}
  
    echo "Getting GTDB-Tk databases for r" ${GTDB_VER_INT}
    export GTDB_TK_DATA_DB=gtdbtk_r${GTDB_VER_INT}_v2_data.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/auxillary_files/${GTDB_TK_DATA_DB}
    tar xzf ${GTDB_TK_DATA_DB} --strip 1
    rm ${GTDB_TK_DATA_DB}
    if [ ! -d mash ] ; then
	mkdir mash
    fi
    chmod 777 mash

    echo "Getting GTDB Archaea metadata"
    export ARC_METADATA=ar53_metadata_r${GTDB_VER_INT}.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${ARC_METADATA}
    tar xzf ${ARC_METADATA}
    rm ${ARC_METADATA}

    echo "Getting GTDB Bacteria metadata"
    export BAC_METADATA=bac120_metadata_r${GTDB_VER_INT}.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${BAC_METADATA}
    tar xzf ${BAC_METADATA}
    rm ${BAC_METADATA}
  fi

  # r214 refdata
  cd /data
  if [ ! -d "/data/r214" ] ; then
    export GTDB_VER_INT=214
    export GTDB_VER_FLT=214.1

    mkdir r${GTDB_VER_INT}
    cd r${GTDB_VER_INT}

    echo "Getting GTDB-Tk databases for r" ${GTDB_VER_INT}
    export GTDB_TK_DATA_DB=gtdbtk_r${GTDB_VER_INT}_data.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/auxillary_files/${GTDB_TK_DATA_DB}
    tar xzf ${GTDB_TK_DATA_DB} --strip 1
    rm ${GTDB_TK_DATA_DB}
    if [ ! -d mash ] ; then
	mkdir mash
    fi
    chmod 777 mash
    
    echo "Getting GTDB Archaea metadata"
    export ARC_METADATA=ar53_metadata_r${GTDB_VER_INT}.tsv.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${ARC_METADATA}
    gunzip ${ARC_METADATA}
    rm ${ARC_METADATA}

    echo "Getting GTDB Bacteria metadata"
    export BAC_METADATA=bac120_metadata_r${GTDB_VER_INT}.tsv.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${BAC_METADATA}
    gunzip ${BAC_METADATA}
    rm ${BAC_METADATA}
  fi

  # don't repeat if refdata prepared
  cd /data
  if [[ -d "r207/taxonomy" && -d "r207/fastani" && -d "r207/markers" && -s "r207/ar53_metadata_r207.tsv" && -s "r207/bac120_metadata_r207.tsv" && -d "r214/taxonomy" && -d "r214/fastani" && -d "r214/markers" && -s "r214/ar53_metadata_r214.tsv" && -s "r214/bac120_metadata_r214.tsv" ]] ;  then
    touch __READY__
  else
    echo "init failed"
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi

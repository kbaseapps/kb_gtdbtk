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
  cd /data
  if [! -d /data/r207] ; then
    mkdir r207
    cd r207

    export GTDB_VER_INT=207
    export GTDB_VER_FLT=207.0
  
    echo "Getting GTDB-Tk databases"
    export GTDB_TK_DATA_DB=gtdbtk_r${GTDB_VER_INT}_v2_data.tar.gz
    curl -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/auxillary_files/${GTDB_TK_DATA_DB}
    tar xvzf ${GTDB_TK_DATA_DB} --strip 1
    rm ${GTDB_TK_DATA_DB}

    echo "Getting GTDB Archaea metadata"
    export ARC_METADATA=ar53_metadata_r${GTDB_VER_INT}.tar.gz
    curl -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${ARC_METADATA}
    tar xvzf ${ARC_METADATA}
    rm ${ARC_METADATA}

    echo "Getting GTDB Bacteria metadata"
    export BAC_METADATA=bac120_metadata_r${GTDB_VER_INT}.tar.gz
    curl -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${BAC_METADATA}
    tar xvzf ${BAC_METADATA}
    rm ${BAC_METADATA}
  fi

  cd /data
  if [[ -d "r207/taxonomy" && -d "r207/fastani" && -d "r207/markers" && -s "r207/ar53_metadata_r207.tsv" && -s "bac120_metadata_r207.tsv" ]] ; then
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

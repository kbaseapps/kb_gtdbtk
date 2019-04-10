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
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  cd /data
  echo "Getting GTDB-Tk database"
  curl -O https://data.ace.uq.edu.au/public/gtdbtk/release_86/gtdbtk.r86_v2_data.tar.gz
  tar xvzf gtdbtk.r86_v2_data.tar.gz --strip 1
  rm gtdbtk.r86_v2_data.tar.gz
  if [[ -d "taxonomy" && -d "fastani" && -d "markers"]] ; then
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

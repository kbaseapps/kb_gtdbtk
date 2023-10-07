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

  ################  
  # r207 refdata #
  ################  
  export GTDB_VER_INT=207
  export GTDB_VER_FLT=207.0
  cd /data
  
  if [ ! -d "/data/r${GTDB_VER_INT}" ] ; then

    mkdir r${GTDB_VER_INT}
    cd /data/r${GTDB_VER_INT}

    # Download gtdb-tk refdata
    echo "Getting GTDB-Tk databases for r"${GTDB_VER_INT}
    export GTDB_TK_DATA_DB=gtdbtk_r${GTDB_VER_INT}_v2_data.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/auxillary_files/${GTDB_TK_DATA_DB}
    tar xzf ${GTDB_TK_DATA_DB} --strip 1
    rm ${GTDB_TK_DATA_DB}

    # make mash db by running gtdbtk (can't do later because /data read-omly)
    if [ ! -d mash ] ; then
	mkdir mash
    fi
    export local_fna_dir="fna_dir"
    export mash_db_file="gtdb_ref_sketch.msh"
    export mash_db_dir="/data/r${GTDB_VER_INT}/mash"
    export mash_db_path="${mash_db_dir}/${mash_db_file}"
    if [ ! -s ${mash_db_path} ] ; then
	echo "Generating GTDB-Tk mash db for r"${GTDB_VER_INT}
	cd /kb/module
	mkdir -p ${local_fna_dir}
	export fna_file="GCF_000008665.1_assembly.fa.gz"
	cp /kb/module/test/data/${fna_file} ${local_fna_dir}
	export GTDBTK_DATA_PATH="/data/r${GTDB_VER_INT}"
	gtdbtk classify_wf --out_dir ./ --genome_dir ${local_fna_dir} --extension gz --cpus 8 --min_perc_aa 10 --mash_db ${mash_db_path} &> mash.log
	cd /data/r${GTDB_VER_INT}
    fi

    # Download archaeal metadata
    echo "Getting GTDB Archaea metadata for r"${GTDB_VER_INT}
    export ARC_METADATA=ar53_metadata_r${GTDB_VER_INT}.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${ARC_METADATA}
    tar xzf ${ARC_METADATA}
    rm ${ARC_METADATA}

    # Download bacterial metadata
    echo "Getting GTDB Bacteria metadata for r"${GTDB_VER_INT}
    export BAC_METADATA=bac120_metadata_r${GTDB_VER_INT}.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${BAC_METADATA}
    tar xzf ${BAC_METADATA}
    rm ${BAC_METADATA}
  fi

    
  ################  
  # r214 refdata #
  ################
  export GTDB_VER_INT=214
  export GTDB_VER_FLT=214.1

  cd /data
  if [ ! -d "/data/r${GTDB_VER_INT}" ] ; then

    mkdir r${GTDB_VER_INT}
    cd r${GTDB_VER_INT}

    # Download gtdb-tk refdata
    echo "Getting GTDB-Tk databases for r"${GTDB_VER_INT}
    export GTDB_TK_DATA_DB=gtdbtk_r${GTDB_VER_INT}_data.tar.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/auxillary_files/${GTDB_TK_DATA_DB}
    tar xzf ${GTDB_TK_DATA_DB} --strip 1
    rm ${GTDB_TK_DATA_DB}

    # make mash db by running gtdbtk (can't do later because /data read-omly)
    if [ ! -d mash ] ; then
	mkdir mash
    fi
    export local_fna_dir="fna_dir"
    export mash_db_file="gtdb_ref_sketch.msh"
    export mash_db_dir="/data/r${GTDB_VER_INT}/mash"
    export mash_db_path="${mash_db_dir}/${mash_db_file}"
    if [ ! -s ${mash_db_path} ] ; then
	echo "Generating GTDB-Tk mash db for r"${GTDB_VER_INT}
	cd /kb/module
	mkdir -p ${local_fna_dir}
	export fna_file="GCF_000008665.1_assembly.fa.gz"
	cp /kb/module/test/data/${fna_file} ${local_fna_dir}
	export GTDBTK_DATA_PATH="/data/r${GTDB_VER_INT}"
	gtdbtk classify_wf --out_dir ./ --genome_dir ${local_fna_dir} --extension gz --cpus 8 --min_perc_aa 10 --mash_db ${mash_db_path} &> mash.log
	cd /data/r${GTDB_VER_INT}
    fi
    
    # Download archaeal metadata
    echo "Getting GTDB Archaea metadata for r"${GTDB_VER_INT}
    export ARC_METADATA=ar53_metadata_r${GTDB_VER_INT}.tsv.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${ARC_METADATA}
    gunzip ${ARC_METADATA}

    # Download bacterial metadata
    echo "Getting GTDB Bacteria metadata for r"${GTDB_VER_INT}
    export BAC_METADATA=bac120_metadata_r${GTDB_VER_INT}.tsv.gz
    curl -s -O https://data.gtdb.ecogenomic.org/releases/release${GTDB_VER_INT}/${GTDB_VER_FLT}/${BAC_METADATA}
    gunzip ${BAC_METADATA}
  fi

  # don't repeat if refdata prepared
  cd /data
  if [[ -d "r207/taxonomy" && -d "r207/fastani" && -d "r207/markers" && -s "r207/mash/gtdb_ref_sketch.msh" && -s "r207/ar53_metadata_r207.tsv" && -s "r207/bac120_metadata_r207.tsv" && -d "r214/taxonomy" && -d "r214/fastani" && -d "r214/markers" && -s "r214/mash/gtdb_ref_sketch.msh" && -s "r214/ar53_metadata_r214.tsv" && -s "r214/bac120_metadata_r214.tsv" ]] ;  then
    touch __READY__
  else
    echo "init failed"

    if [ ! -d "r207/taxonomy" ] ; then
	echo "missing r207/taxonomy"
    fi
    if [ ! -d "r207/fastani" ] ; then
	echo "missing r207/fastani"
    fi
    if [ ! -d "r207/markers" ] ; then
	echo "missing r207/markers"
    fi
    if [ ! -s "r207/mash/gtdb_ref_sketch.msh" ] ; then
	echo "missing r207/mash/gtdb_ref_sketch.msh"
    fi
    if [ ! -s "r207/ar53_metadata_r207.tsv" ] ; then
	echo "missing r207/ar53_metadata_r207.tsv"
    fi
    if [ ! -s "r207/bac120_metadata_r207.tsv" ] ; then
	echo "missing r207/bac120_metadata_r207.tsv"
    fi
    if [ ! -d "r214/taxonomy" ] ; then
	echo "missing r214/taxonomy"
    fi
    if [ ! -d "r214/fastani" ] ; then
	echo "missing r214/fastani"
    fi
    if [ ! -d "r214/markers" ] ; then
	echo "missing r214/markers"
    fi
    if [ ! -s "r214/mash/gtdb_ref_sketch.msh" ] ; then
	echo "missing r214/mash/gtdb_ref_sketch.msh"
    fi
    if [ ! -s "r214/ar53_metadata_r214.tsv" ] ; then
	echo "missing r214/ar53_metadata_r214.tsv"
    fi
    if [ ! -s "r214/bac120_metadata_r214.tsv" ] ; then
	echo "missing r214/bac120_metadata_r214.tsv"
    fi

  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi

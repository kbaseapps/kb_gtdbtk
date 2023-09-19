SERVICE = kb_gtdbtk
SERVICE_CAPS = kb_gtdbtk
SPEC_FILE = kb_gtdbtk.spec
URL = https://kbase.us/services/kb_gtdbtk
DIR = $(shell pwd)
LIB_DIR = lib
SCRIPTS_DIR = scripts
TEST_DIR = test
LBIN_DIR = bin
WORK_DIR = /kb/module/work/tmp
EXECUTABLE_SCRIPT_NAME = run_$(SERVICE_CAPS)_async_job.sh
STARTUP_SCRIPT_NAME = start_server.sh
TEST_SCRIPT_NAME = run_tests.sh
# see https://stackoverflow.com/a/23324703/643675
MAKEFILE_DIR:=$(strip $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST)))))
PYPATH=$(MAKEFILE_DIR)/$(LIB_DIR):$(MAKEFILE_DIR)/$(TEST_DIR)

.PHONY: test

default: compile

all: compile build build-startup-script build-executable-script build-test-script

compile:
	kb-sdk compile $(SPEC_FILE) \
		--out $(LIB_DIR) \
		--pysrvname $(SERVICE_CAPS).$(SERVICE_CAPS)Server \
		--pyimplname $(SERVICE_CAPS).$(SERVICE_CAPS)Impl;
	
	kb-sdk compile $(SPEC_FILE) \
		--out . \
		--html \

build:
	chmod +x $(SCRIPTS_DIR)/entrypoint.sh

build-executable-script:
	mkdir -p $(LBIN_DIR)
	echo '#!/bin/bash' > $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'script_dir=$$(dirname "$$(readlink -f "$$0")")' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export PYTHONPATH=$$script_dir/../$(LIB_DIR):$$PATH:$$PYTHONPATH' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'python -u $$script_dir/../$(LIB_DIR)/$(SERVICE_CAPS)/$(SERVICE_CAPS)Server.py $$1 $$2 $$3' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	chmod +x $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)

build-startup-script:
	mkdir -p $(LBIN_DIR)
	echo '#!/bin/bash' > $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'script_dir=$$(dirname "$$(readlink -f "$$0")")' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'export KB_DEPLOYMENT_CONFIG=$$script_dir/../deploy.cfg' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'export PYTHONPATH=$$script_dir/../$(LIB_DIR):$$PATH:$$PYTHONPATH' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'uwsgi --master --processes 5 --threads 5 --http :5000 --wsgi-file $$script_dir/../$(LIB_DIR)/$(SERVICE_CAPS)/$(SERVICE_CAPS)Server.py' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	chmod +x $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)

build-test-script:
	echo '#!/bin/bash' > $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'script_dir=$$(dirname "$$(readlink -f "$$0")")' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'export KB_DEPLOYMENT_CONFIG=$$script_dir/../deploy.cfg' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'export KB_AUTH_TOKEN=`cat /kb/module/work/token`' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'echo "Removing temp files..."' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'rm -rf $(WORK_DIR)/*' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'echo "...done removing temp files."' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	# TODO flake8 and bandit
	echo 'MYPYPATH=$(MAKEFILE_DIR)/$(LIB_DIR) mypy --namespace-packages $(LIB_DIR)/$(SERVICE_CAPS)/core' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'PYTHONPATH=$(PYPATH) pytest --verbose --cov $(LIB_DIR)/$(SERVICE_CAPS) --cov-config=$(TEST_DIR)/coveragerc $(TEST_DIR) --cov-report=html' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'mv .coverage /kb/module/work/' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'rm -r /kb/module/work/test_coverage' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'mv htmlcov /kb/module/work/test_coverage' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	chmod +x $(TEST_DIR)/$(TEST_SCRIPT_NAME)

test:
	if [ ! -f /kb/module/work/token ]; then echo -e '\nOutside a docker container please run "kb-sdk test" rather than "make test"\n' && exit 1; fi
	xvfb-run bash $(TEST_DIR)/$(TEST_SCRIPT_NAME)

test-sdkless:
	# TODO flake8 and bandit
	MYPYPATH=$(MAKEFILE_DIR)/$(LIB_DIR) mypy --namespace-packages $(LIB_DIR)/$(SERVICE_CAPS)/core
	PYTHONPATH=$(PYPATH) pytest --verbose --cov $(LIB_DIR)/$(SERVICE_CAPS) --cov-config=$(TEST_DIR)/coveragerc $(TEST_DIR) --ignore=$(TEST_DIR)/kb_gtdbtk_server_test.py

clean:
	rm -rfv $(LBIN_DIR)

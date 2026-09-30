SHELL := /bin/bash

SCRIPT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
PARTICIPANTS := tissue_41_diffusion tissue_42_decay tissue_43_cell_response tissue_44_twoway
BDM_PATH ?=

.SUFFIXES:
.PHONY: help build $(PARTICIPANTS) clean

help:
	@echo "make build                         Build the cells participant"
	@echo "make run-<case>                    Run one of: $(PARTICIPANTS)"
	@echo "make clean                         Remove local build/output artifacts"
	@echo "Set BDM_PATH=/path/to/thisbdm.sh unless BioDynaMo is already active."

build:
	@bash -c 'if [[ -z "$${BDMSYS:-}" ]]; then \
		[[ -n "$(BDM_PATH)" && -f "$(BDM_PATH)" ]] || { \
			echo "Set BDM_PATH to BioDynaMo thisbdm.sh" >&2; exit 1; }; \
		source "$(BDM_PATH)"; fi; \
		cd "$(SCRIPT_DIR)" && bdm build'

$(PARTICIPANTS): build
	@echo "Build complete for $@"

# Explicit rules, not a run-% pattern. Listing these names in .PHONY disables
# pattern-rule search, and Make then reports "Nothing to be done".
define RUN_CASE
.PHONY: run-$(1)
run-$(1):
	@cd "$(SCRIPT_DIR)" && BDM_PATH="$(BDM_PATH)" ./run.sh "$(1)"
endef
$(foreach case,$(PARTICIPANTS),$(eval $(call RUN_CASE,$(case))))

clean:
	@cd "$(SCRIPT_DIR)" && ./clean.sh

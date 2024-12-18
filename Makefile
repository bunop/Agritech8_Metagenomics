
# Start by specifying the project names and their corresponding scripts and stores.
PROJECTS := technical_replicates plot_iNEXT reactor_vs_algae duckweed_and_box

# Create a generic rule that sets the TAR_PROJECT environment variable and runs
# tar_make() for the specified project. This approach eliminates redundancy and
# keeps the Makefile concise.
.PHONY: $(PROJECTS) restore list

$(PROJECTS):
	@echo "Building project: $@"
	@Rscript -e 'Sys.setenv(TAR_PROJECT = "$@"); targets::tar_make()'

restore:
	@echo "Restoring R environment..."
	@Rscript -e 'renv::restore()'

# Rule to list all available projects
list:
	@echo "Available projects:"
	@for project in $(PROJECTS); do echo "  $$project"; done

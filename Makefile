
# Start by specifying the project names and their corresponding scripts and stores.
PROJECTS := technical_replicates plot_iNEXT reactor_vs_algae duckweed_and_box

# Create a generic rule that sets the TAR_PROJECT environment variable and runs
# tar_make() for the specified project. This approach eliminates redundancy and
# keeps the Makefile concise.
.PHONY: $(PROJECTS) restore list krona_bacteria krona_fungi

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

# krona rules
krona_bacteria:
	@echo "Calling Krona on bacteria samples"
	@set -e; Rscript scripts/run_krona.R -i results-bacteria/phyloseq/dada2_phyloseq.rds \
		-m data/metadata_bacteria_fix.tsv -o results-bacteria-krona

krona_fungi:
	@echo "Calling Krona on fungi samples"
	@set -e; Rscript scripts/run_krona.R -i results-fungi/phyloseq/dada2_phyloseq.rds \
		-m data/metadata_fungi_fix.tsv -o results-fungi-krona

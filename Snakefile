

rule all:
	input:
		"results/figures/figure-001.pdf",
		"docs/001-introductory-linkage-sims.html",
		"docs/003-permutation-methods-figure.html",





rule introductory_linkage_simulations:
	input:
		Rmd = "Rmd/admix-fracts-of-hybrids.Rmd"
	output:
		html = "docs/001-introductory-linkage-sims.html",
		fig = "results/figures/figure-001.pdf"
	log:
		"results/logs/introductory_linkage_simulations/log.txt"
	threads: 8
	conda:
		"envs/rendering.yaml"
	script:
		"R/render-rmd-for-snakemake.R"


rule figure_of_permutation_methods:
	input:
		Rmd = "Rmd/permutation-methods-figure.Rmd"
	output:
		html = "docs/003-permutation-methods-figure.html",
		#fig = "results/figures/figure-003.pdf"
	log:
		"results/logs/figure_of_permutation_methods/log.txt"
	threads: 1
	conda:
		"envs/rendering.yaml"
	script:
		"R/render-rmd-for-snakemake.R"

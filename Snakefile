
freqs=["A"]
Qs=["1"]
Ls=["100", "1000", "10000", "100000"]
Ns=["25", "50", "100", "250"]
ns=["1", "10", "50"]



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


rule simulate_bias_sim_data_sets:
	params:
		freq="{freq}",
		Q="{Q}",
		L="{nloc}",
		N="{N}",
		n="{n}",
		rep="{rep}",
		prefix="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink"
	log:
		"results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/R-sim.log"
	envmodules:
		"R/4.0.3" 
	output:
		multiext("results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink", ".ped", ".map", ".pop")
	script:
		"R/bias-sim-script.R"


rule run_admixture:
	input:
		ped="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink.ped"
	params:
		plink_in="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink",
		bed_out="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink",
		seed_file="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/R.seed",
		bed_in="plink.bed",
		admix_out="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink.2.Q",
		key="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/key.tsv",
		ad_dir="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}", # admixture writes to current working directory.  Lame!
	output:
		sup="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/supervised_Q.tsv",
		uns="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/unsupervised_Q.tsv"
	log:
		plink="results/bias_sims/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/log_of_plink.txt",
	conda:
		"envs/plink2_admixture.yaml"
	shell:
		" plink2 --make-bed --pedmap {params.plink_in} --out {params.bed_out} > {log.plink} 2>&1;  "
		" RSEED=$(cat {params.seed_file});  "
		" CUR_DIR=$(pwd);  "
		" cd {params.ad_dir}; admixture  -s $(($RSEED + 1)) {params.bed_in}  2   > log_of_admixture_uns.txt 2>&1; cd $CUR_DIR; "
		" paste {params.key} {params.admix_out} > {output.uns}; "
		" cd {params.ad_dir}; admixture  -s $(($RSEED + 2)) --supervised {params.bed_in}  2  > log_of_admixture_sup.txt 2>&1; cd $CUR_DIR; "
		" paste {params.key} {params.admix_out} > {output.sup}; "

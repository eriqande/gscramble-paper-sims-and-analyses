
# Here are the different values for the simulations.

freqs=["A"]  # in the R script, this translates to: list(a = 1, b = 8, t = 0.01)
Qs=["1", "1_sim_the_knowns"]     # in the R script, these translate to: Qs <- c(0, 0.125, 0.25, 0.375, 0.5)
Ls=["100", "1000", "10000", "100000"]  # number of loci
Ns=["25", "50", "100", "250"]          # size of each reference population
ns=["3", "12", "24"]                   # number of individuals simulated for each Q value

# Note, for the gscramble simulations, most of the hybrid categories
# get two individuals simulated.  So, we will want 240 simulations
# of those, and we will just set their n to 2.  

# we want to devise these simulations so that for each
# value of L, N, and n, and Q, we have a total of 480
# simulated individuals.  That means, for:
#  - n = 3: reps = 160
#  - n = 12: reps = 40
#  - n = 24: reps = 20

# so, for testing, we could do 16, 4, and 2.  Or even smaller value.
# I can set those value here:
R3=160
R12=40
R24=20

# Here, we make a list of all the different Q-output files that
# we want to have as input to make one giant, gzipped text
# file of the output for tidyverse processing:
threes = expand("results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/{cond}vised_Q.tsv",
	bias_sims="bias_sims", freq=freqs, Q=Qs, nloc=Ls, N=Ns, n=3, rep=range(1,R3+1), cond=['super', 'unsuper'])
twelves = expand("results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/{cond}vised_Q.tsv",
	bias_sims="bias_sims", freq=freqs, Q=Qs, nloc=Ls, N=Ns, n=12, rep=range(1,R12+1), cond=['super', 'unsuper'])
twenty_fours = expand("results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/{cond}vised_Q.tsv",
	bias_sims="bias_sims", freq=freqs, Q=Qs, nloc=Ls, N=Ns, n=24, rep=range(1,R24+1), cond=['super', 'unsuper'])

# and here we request the gscramble sims
gscrambles = expand("results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/{cond}vised_Q.tsv",
	bias_sims="gscramble", freq=freqs, Q=1, nloc=[1000], N=50, n=2, rep=range(1,24+1), cond=['super', 'unsuper'])
	
# and these are the same thing, but with strongly diverged markers
# so we can verify that our simulations and gscramble are working correctly
gscramble_diffs = expand("results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/{cond}vised_Q.tsv",
	bias_sims="gscramble_diff", freq=freqs, Q=1, nloc=[1000], N=50, n=2, rep=range(1,24+1), cond=['super', 'unsuper'])



# this little function doubles the memory requested when
# N=250, L=100000
def mem_func(wildcards):
	if(wildcards.N=="250" and wildcards.nloc=="100000"):
		return 9600
	else:
		return 4800


rule all:
	input:
		#"results/figures/figure-001.pdf",
		"results/compiled/Q-values-from-sims.tsv.gz",
		#"docs/001-introductory-linkage-sims.html",
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
		prefix="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink",
		sim_type="{bias_sims}"
	log:
		log="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/R-sim.log",
		snake_obj="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/snake_obj.rds"
	envmodules:
		"R/4.0.3"
	resources:
		mem_mb=mem_func 
	output:
		multiext("results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink", ".ped", ".map", ".pop")
	script:
		"R/bias-sim-script.R"


rule run_admixture:
	input:
		ped="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink.ped"
	params:
		plink_in="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink",
		bed_out="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink",
		seed_file="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/R.seed",
		bed_in="plink.bed",
		admix_out="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/plink.2.Q",
		key="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/key.tsv",
		ad_dir="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}", # admixture writes to current working directory.  Lame!
		sim_type="{bias_sims}"
	output:
		sup="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/supervised_Q.tsv",
		uns="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/unsupervised_Q.tsv"
	log:
		plink="results/{bias_sims}/freq_{freq}/Qs_{Q}/L_{nloc}/N_{N}/n_{n}/Rep_{rep}/log_of_plink.txt",
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


rule compile_Qs:
	input:
		Qfiles=gscrambles + gscramble_diffs # + threes + twelves + twenty_fours
	output:
		"results/compiled/Q-values-from-sims.tsv.gz"
	log:
		"results/logs/compile_Qs.txt"
	shell:
		" (for i in {input.Qfiles}; do "
		" awk -v file=$i 'BEGIN {{OFS=\"\t\"}} {{print file, $0}}' $i; "
		" done | gzip -c > {output}) 2> {log} "

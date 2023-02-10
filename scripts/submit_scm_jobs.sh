#-- Launch slurm job array
echo -e "Launching analyses... "
sbatch -p lemmium --array="1-$(wc -l scripts/arg_list_scm.txt | cut -d ' ' -f1)" ./scripts/job.sh
echo "done"

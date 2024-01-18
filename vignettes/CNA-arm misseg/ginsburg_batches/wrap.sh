for i in {1..30}; do
    sbatch "batch_program_fitting_arm_step3_stats_${i}.sh"
    sleep 1
done
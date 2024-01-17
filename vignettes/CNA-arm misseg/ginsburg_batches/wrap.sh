for i in {1..100}; do
    sbatch "batch_program_fitting_arm_step2_library_${i}.sh"
    sleep 1
done
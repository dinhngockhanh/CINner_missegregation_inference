list <- c(2:30)
for (i in list) {
    #---Copy and reconfigure batch file
    newfilename <- paste0("batch_program_fitting_arm_step3_stats_", i, ".sh")
    file.copy("batch_program_fitting_arm_step3_stats_1.sh", newfilename)
    tx <- readLines(newfilename)
    tx <- gsub(pattern = "program_fitting_arm_step3_stats_1.r", replace = paste0("program_fitting_arm_step3_stats_", i, ".r"), x = tx)
    tx <- gsub(pattern = "routput_stats_1", replace = paste0("routput_stats_", i), x = tx)
    writeLines(tx, con = newfilename)
    #---Copy and reconfigure R file
    newfilename <- paste0("program_fitting_arm_step3_stats_", i, ".r")
    file.copy("program_fitting_arm_step3_stats_1.r", newfilename)
    tx <- readLines(newfilename)
    tx <- gsub(pattern = "ABC_simcount_start = 0,", replace = paste0("ABC_simcount_start = ", (i - 1) * 500, ","), x = tx)
    writeLines(tx, con = newfilename)
}

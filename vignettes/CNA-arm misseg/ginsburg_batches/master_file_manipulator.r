list <- c(2:100)
for (i in list) {
    #---Copy and reconfigure batch file
    newfilename <- paste0("batch_program_fitting_arm_step2_library_", i, ".sh")
    file.copy("batch_program_fitting_arm_step2_library_1.sh", newfilename)
    tx <- readLines(newfilename)
    tx <- gsub(pattern = "program_fitting_arm_step2_library_1.r", replace = paste0("program_fitting_arm_step2_library_", i, ".r"), x = tx)
    tx <- gsub(pattern = "routput_lib_1", replace = paste0("routput_lib_", i), x = tx)
    writeLines(tx, con = newfilename)
    #---Copy and reconfigure R file
    newfilename <- paste0("program_fitting_arm_step2_library_", i, ".r")
    file.copy("program_fitting_arm_step2_library_1.r", newfilename)
    tx <- readLines(newfilename)
    tx <- gsub(pattern = "ABC_simcount_start = 0,", replace = paste0("ABC_simcount_start = ", (i - 1) * 1000, ","), x = tx)
    writeLines(tx, con = newfilename)
}

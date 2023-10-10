# sensitivity_values_bulk <- seq(10, 100, by = 10)
# for (i in sensitivity_values_bulk) {
#     setwd(paste0("/Users/xiangzijin/Downloads/Fitting_whole_chroms_N_data_bulk_", i))
#     filename <- paste0("Fitting_whole_chroms_N_data_bulk_", i, "_ABC_input_1.rda")
#     load(filename)
#     ABC_input_all <- ABC_input
#     for (batch in 2:5) {
#         filename <- paste0("Fitting_whole_chroms_N_data_bulk_", i, "_ABC_input_", batch, ".rda")
#         load(filename)
#         ABC_input_all$sim_param <- rbind(ABC_input_all$sim_param, ABC_input$sim_param)
#         for (j in 1:length(ABC_input$sim_stat)) {
#             ABC_input_all$sim_stat[[j]] <- rbind(ABC_input_all$sim_stat[[j]], ABC_input$sim_stat[[j]])
#         }
#     }
#     filename <- paste0("Fitting_whole_chroms_N_data_bulk_", i, "_ABC_input.rda")
#     ABC_input <- ABC_input_all
#     setwd("/Users/xiangzijin/Downloads")
#     save(ABC_input, file = filename)
# }

sensitivity_values_sc <- seq(5, 50, by = 5)
for (i in sensitivity_values_sc) {
    setwd(paste0("/Users/xiangzijin/Downloads/Fitting_whole_chroms_N_data_sc_", i))
    filename <- paste0("Fitting_whole_chroms_N_data_sc_", i, "_ABC_input_1.rda")
    load(filename)
    ABC_input_all <- ABC_input
    for (batch in 2:250) {
        filename <- paste0("Fitting_whole_chroms_N_data_sc_", i, "_ABC_input_", batch, ".rda")
        load(filename)
        ABC_input_all$sim_param <- rbind(ABC_input_all$sim_param, ABC_input$sim_param)
        for (j in 1:length(ABC_input$sim_stat)) {
            ABC_input_all$sim_stat[[j]] <- rbind(ABC_input_all$sim_stat[[j]], ABC_input$sim_stat[[j]])
        }
    }
    filename <- paste0("Fitting_whole_chroms_N_data_sc_", i, "_ABC_input.rda")
    ABC_input <- ABC_input_all
    setwd("/Users/xiangzijin/Downloads")
    save(ABC_input, file = filename)
}

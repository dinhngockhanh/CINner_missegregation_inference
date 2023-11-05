sensitivity_values_bulk <- seq(10, 100, by = 10)
original <- c(2:500)
subset <- c(
    19, 20, 24, 27, 38, 39, 47, 54, 61, 66, 68, 73, 76, 81, 97, 122, 132, 135, 146, 151,
    155, 159, 162, 190, 193, 205, 209, 223, 248, 279, 288, 289, 310, 317, 325, 379, 396, 397, 409,
    412, 418, 421, 429, 436, 446, 463, 483, 496
)
new_bulk <- original[!original %in% subset]
for (i in sensitivity_values_bulk) {
    setwd(paste0("/Users/xiangzijin/Downloads/Fitting_whole_chroms_N_data_bulk_", i))
    filename <- paste0("Fitting_whole_chroms_N_data_bulk_", i, "_ABC_input_1.rda")
    load(filename)
    ABC_input_all <- ABC_input
    for (batch in new_bulk) {
        filename <- paste0("Fitting_whole_chroms_N_data_bulk_", i, "_ABC_input_", batch, ".rda")
        load(filename)
        ABC_input_all$sim_param <- rbind(ABC_input_all$sim_param, ABC_input$sim_param)
        for (j in 1:length(ABC_input$sim_stat)) {
            ABC_input_all$sim_stat[[j]] <- rbind(ABC_input_all$sim_stat[[j]], ABC_input$sim_stat[[j]])
        }
    }
    filename <- paste0("Fitting_whole_chroms_N_data_bulk_", i, "_ABC_input_partial.rda")
    ABC_input <- ABC_input_all
    setwd("/Users/xiangzijin/Downloads")
    save(ABC_input, file = filename)
}

# sensitivity_values_sc <- seq(5, 50, by = 5)
# for (i in sensitivity_values_sc) {
#     setwd(paste0("/Users/xiangzijin/Downloads/Fitting_whole_chroms_N_data_sc_", i))
#     filename <- paste0("Fitting_whole_chroms_N_data_sc_", i, "_ABC_input_1.rda")
#     load(filename)
#     ABC_input_all <- ABC_input
#     for (batch in 2:250) {
#         filename <- paste0("Fitting_whole_chroms_N_data_sc_", i, "_ABC_input_", batch, ".rda")
#         load(filename)
#         ABC_input_all$sim_param <- rbind(ABC_input_all$sim_param, ABC_input$sim_param)
#         for (j in 1:length(ABC_input$sim_stat)) {
#             ABC_input_all$sim_stat[[j]] <- rbind(ABC_input_all$sim_stat[[j]], ABC_input$sim_stat[[j]])
#         }
#     }
#     filename <- paste0("Fitting_whole_chroms_N_data_sc_", i, "_ABC_input.rda")
#     ABC_input <- ABC_input_all
#     setwd("/Users/xiangzijin/Downloads")
#     save(ABC_input, file = filename)
# }



# new_bulk

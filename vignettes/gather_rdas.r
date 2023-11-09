sensitivity_values_bulk <- seq(10, 100, by = 10)
original <- c(2:500)
subset <- c(
    19, 20, 24, 27, 37:40, 47:48, 53, 54, 61, 66, 68, 73, 76, 81, 97, 121, 122, 131, 132, 135, 136, 145, 146, 151, 152,
    155, 159, 161, 162, 190, 193, 194, 205, 209, 223, 248, 279, 288, 289, 209, 309, 310, 317, 318, 325, 379, 380, 396, 397, 409,
    410, 412, 417, 418, 421, 429, 436, 445, 446, 463, 483, 495, 496
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

#---Gather simulation librarys of 1000&500
filename <- paste0("Simpler_DLP&BULK_DNA", "_ABC_input1.rda")
load(filename)
ABC_input_all <- ABC_input
for (batch in 2:19) {
    filename <- paste0("Simpler_DLP&BULK_DNA", "_ABC_input", batch, ".rda")
    load(filename)
    ABC_input_all$sim_param <- rbind(ABC_input_all$sim_param, ABC_input$sim_param)
    for (i in 1:length(sim_stat)) {
        ABC_input_all$sim_stat[[i]] <- rbind(ABC_input_all$sim_stat[[i]], ABC_input$sim_stat[[i]])
    }
}
filename <- "Simpler_DLP&BULK_DNA_ABC_input.rda"
save(ABC_input_all, file = filename)

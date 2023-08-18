#---Gather simulation librarys of 1000&500
ABC_input_model_variables <- list()
ABC_input_n_samples <- 0
ABC_input_parameter_IDs <- c()
ABC_input_sim_param <- c()
ABC_input_sim_stat <- c()
ABC_input_sim_sample_stat <- list()
for (batch in 1:20) {
    ####################################################################
    filename <- paste0("Simpler_DLP&BULK_DNA", "_ABC_input", batch, ".rda")
    ####################################################################
    load(filename)
    ABC_input_model_variables <- ABC_input$model_variables
    ABC_input_n_samples <- ABC_input$n_samples
    ABC_input_parameter_IDs <- ABC_input$parameter_IDs
    ABC_input_sim_param <- rbind(ABC_input_sim_param, ABC_input$sim_param)
    ABC_input_sim_stat <- rbind(ABC_input_sim_stat, ABC_input$sim_stat)
    ABC_input_sim_sample_stat <- c(ABC_input_sim_sample_stat, ABC_input$sim_sample_stat)
    load(filename)
}

ABC_input <- list()
ABC_input$model_variables <- ABC_input_model_variables
ABC_input$n_samples <- ABC_input_n_samples
ABC_input$parameter_IDs <- ABC_input_parameter_IDs
ABC_input$sim_param <- ABC_input_sim_param
ABC_input$sim_stat <- ABC_input_sim_stat
ABC_input$sim_sample_stat <- ABC_input_sim_sample_stat
filename <- "Simpler_DLP&BULK_DNA_ABC_input.rda"
save(ABC_input, file = filename)

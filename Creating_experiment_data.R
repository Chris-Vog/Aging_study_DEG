# Setting working directory
setwd("Timur/")
list.files()

# Loading required packages
suppressMessages(library(readxl))
suppressMessages(library(dplyr))

# Creating vector with column names
col_names <- c("ID" , "ERR" , "ERS" , "donor" , "age" , "year" , "sex" , "site_of_biopsis" , "file")

# Loading experiment data
experiment <- read_xlsx(path = "Zuordnung samples.xlsx", col_names = FALSE)
colnames(experiment) <- col_names

# Removing unnecessary columns
experiment$ERR <- NULL
experiment$ERS <- NULL
experiment$year <- NULL
experiment$file <- NULL
experiment$donor <- NULL

# Categorize donors into age groups
experiment <- mutate(experiment, age_groups = if_else(age < 26, "young", 
                                                      if_else(age > 59, "old", "middle")))

experiment$age_groups <- factor(experiment$age_groups, labels = c("young" , "middle" , "old"), ordered = TRUE)

# Remove donors without expression data
processed_data <- list.files("Analysis/STAR/")
missing_data <- base::intersect(processed_data, experiment$ID)
experiment <- experiment[experiment$ID %in% missing_data,]

# Saving experiment data
write.table(experiment, file = "experiment.txt", sep = "," , col.names = TRUE)

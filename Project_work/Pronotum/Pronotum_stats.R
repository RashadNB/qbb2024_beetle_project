# Step 1: Add File IDs to TPS Data
add_file_id <- function(tps_data, file_name) {
  modified_data <- c()
  entry_buffer <- c()  # Buffer to store data until ID line
  
  for (line in tps_data) {
    entry_buffer <- c(entry_buffer, line)
    
    # Insert File ID when ID line is reached
    if (grepl("^ID=", line)) {
      entry_buffer <- c(entry_buffer, paste("Original_File_ID=", file_name, sep = ""))
      modified_data <- c(modified_data, entry_buffer)
      entry_buffer <- c()  # Clear the buffer for next entry
    }
  }
  return(modified_data)
}

# List of Pronotum TPS files
tps_files <- c("Aesalini.TPS", "Ceratognathini.TPS", "Chiasognathini.TPS", 
               "Lamprimini.TPS", "Lucanini.TPS", "Nicagini.TPS", "Platycerini.TPS",
               "Platyceroidini.TPS", "Streptocerini.TPS", "Syndesinae.TPS")

# Combine all TPS data with file IDs
combined_tps_data <- c()

for (tps_file in tps_files) {
  tps_data <- readLines(tps_file)
  file_name <- gsub("\\.TPS$", "", tps_file)
  modified_data <- add_file_id(tps_data, file_name)
  combined_tps_data <- c(combined_tps_data, modified_data)
}
writeLines(combined_tps_data, "Pronotum.TPS")  # Save combined data
cat("Modified TPS data saved to: Pronotum.TPS\n")

# Step 2: Load and Prepare Data
library(geomorph)
library(ggplot2)
library(dplyr)
library(RRPP)

# Load the TPS data
tps_file <- "Pronotum.TPS"
landmark_data <- readland.tps(tps_file, specID = "ID")

# Extract tribe from the Original_File_ID
lines <- readLines(tps_file)
file_origin <- factor(gsub("^Original_File_ID=", "", lines[grep("^Original_File_ID=", lines)]))
tribe <- file_origin  # Assuming Tribe corresponds to file origin

# Step 3: Perform Procrustes Alignment
gpa_result <- gpagen(landmark_data)

# Step 4: Procrustes ANOVA
group_test <- procD.lm(gpa_result$coords ~ tribe, iter = 999)
anova_summary <- summary(group_test)
print(anova_summary)  # Display overall ANOVA results

# Save ANOVA results
anova_table <- as.data.frame(anova_summary$table)
write.csv(anova_table, "Procrustes_ANOVA_results_pronotum.csv", row.names = FALSE)
cat("ANOVA results saved to 'Procrustes_ANOVA_results_pronotum.csv'\n")

# Step 5: Pairwise Procrustes ANOVA
pairwise_results <- pairwise(group_test, groups = tribe)
pairwise_summary <- summary(pairwise_results)
# Step 5: Pairwise Procrustes ANOVA
pairwise_results <- pairwise(group_test, groups = tribe)

# Check if pairwise_results contains a valid pairwise table
if (!is.null(pairwise_results$pairwise.table) && nrow(pairwise_results$pairwise.table) > 0) {
  pairwise_table <- as.data.frame(pairwise_results$pairwise.table)
  
  # Adjust p-values and add significance levels
  pairwise_table$Adjusted_P <- p.adjust(pairwise_table$p.value, method = "BH")
  pairwise_table$Significance <- case_when(
    pairwise_table$Adjusted_P <= 0.001 ~ "***",
    pairwise_table$Adjusted_P <= 0.01 ~ "**",
    pairwise_table$Adjusted_P <= 0.05 ~ "*",
    TRUE ~ ""
  )
  
  # Save pairwise results to CSV
  write.csv(pairwise_table, "Pairwise_Procrustes_results_pronotum.csv", row.names = FALSE)
  cat("Pairwise results saved to 'Pairwise_Procrustes_results_pronotum.csv'\n")
} else {
  cat("No valid pairwise comparisons were found.\n")
}

# Extract pairwise results
pairwise_table <- as.data.frame(pairwise_summary$pairwise.table)
pairwise_table$Adjusted_P <- p.adjust(pairwise_table$p.value, method = "BH")  # Adjust p-values
pairwise_table$Significance <- case_when(
  pairwise_table$Adjusted_P <= 0.001 ~ "***",
  pairwise_table$Adjusted_P <= 0.01 ~ "**",
  pairwise_table$Adjusted_P <= 0.05 ~ "*",
  TRUE ~ ""
)
write.csv(pairwise_table, "Pairwise_Procrustes_results_pronotum.csv", row.names = FALSE)
cat("Pairwise results saved to 'Pairwise_Procrustes_results_pronotum.csv'\n")

# Step 6: Heatmap for Pairwise Comparisons
# Prepare heatmap data
tribes <- levels(tribe)
heatmap_data <- expand.grid(Tribe1 = tribes, Tribe2 = tribes)
heatmap_data$PValue <- NA

for (i in 1:(length(tribes) - 1)) {
  for (j in (i + 1):length(tribes)) {
    pair <- paste(tribes[i], tribes[j], sep = " - ")
    if (pair %in% rownames(pairwise_table)) {
      p_val <- pairwise_table$p.value[rownames(pairwise_table) == pair]
      heatmap_data$PValue[heatmap_data$Tribe1 == tribes[i] & heatmap_data$Tribe2 == tribes[j]] <- p_val
      heatmap_data$PValue[heatmap_data$Tribe1 == tribes[j] & heatmap_data$Tribe2 == tribes[i]] <- p_val
    }
  }
}

# Add significance stars
heatmap_data$Significance <- case_when(
  heatmap_data$PValue <= 0.001 ~ "***",
  heatmap_data$PValue <= 0.01 ~ "**",
  heatmap_data$PValue <= 0.05 ~ "*",
  TRUE ~ ""
)

# Plot the heatmap
heatmap_data$PValue[is.na(heatmap_data$PValue)] <- 1  # Replace NA with non-significant
ggplot(heatmap_data, aes(x = Tribe1, y = Tribe2, fill = PValue)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Significance), color = "black", size = 3) +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 1), name = "P-Value") +
  labs(title = "Significant Pairwise Procrustes ANOVA Heatmap for Pronotum",
       x = "Tribe 1", y = "Tribe 2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save filtered data for heatmap
write.csv(heatmap_data, "Significant_Heatmap_Data_Pronotum.csv", row.names = FALSE)
cat("Significant heatmap data saved to 'Significant_Heatmap_Data_Pronotum.csv'\n")








# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Step 1: Load the pairwise results
pairwise_data <- read.csv("Pairwise_Procrustes_results_pronotum.csv", check.names = FALSE, row.names = 1)

# Step 2: Extract p-value columns
# Assuming columns with "P." prefix contain p-values
p_value_columns <- grep("^P\\.", colnames(pairwise_data), value = TRUE)

# Convert the p-value matrix into a long format
heatmap_data <- pairwise_data %>%
  select(all_of(p_value_columns)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Tribe2",
    values_to = "PValue"
  )

# Step 3: Add `Tribe1` column using row names from `pairwise_data`
heatmap_data <- heatmap_data %>%
  mutate(
    Tribe1 = rep(rownames(pairwise_data), each = length(p_value_columns)),
    Tribe2 = sub("^P\\.", "", Tribe2)  # Remove "P." prefix from column names
  )

# Step 4: Add significance stars based on p-values
heatmap_data <- heatmap_data %>%
  mutate(Significance = case_when(
    PValue <= 0.001 ~ "***",
    PValue <= 0.01 ~ "**",
    PValue <= 0.05 ~ "*",
    TRUE ~ ""
  ))

# Step 5: Plot the heatmap
ggplot(heatmap_data, aes(x = Tribe1, y = Tribe2, fill = PValue)) +
  geom_tile(color = "white") +  # Add white tile borders
  geom_text(aes(label = Significance), color = "black", size = 3) +  # Add significance stars
  scale_fill_gradient(low = "white", high = "red", name = "P-Value", na.value = "gray") +  # Color gradient
  labs(title = "Pairwise Procrustes ANOVA P-Value Heatmap",
       x = "Tribe 1",
       y = "Tribe 2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(1.5, "lines")
  )

# Save the heatmap data to a CSV for further inspection if needed
write.csv(heatmap_data, "Procrustes_PValue_Heatmap_Data.csv", row.names = FALSE)



# Ensure Tribe1 gets the correct tribe names
tribe_names <- unique(heatmap_data$Tribe2)  # Use Tribe2 to get the unique names
heatmap_data$Tribe1 <- rep(tribe_names, each = length(tribe_names))  # Repeat tribe names for each row

# Convert Tribe1 and Tribe2 to factors
heatmap_data$Tribe1 <- factor(heatmap_data$Tribe1, levels = tribe_names)
heatmap_data$Tribe2 <- factor(heatmap_data$Tribe2, levels = tribe_names)

# Plot the heatmap
ggplot(heatmap_data, aes(x = Tribe1, y = Tribe2, fill = PValue)) +
  geom_tile(color = "white") +  # Add white borders to tiles
  geom_text(aes(label = Significance), color = "black", size = 3) +  # Add significance stars
  scale_fill_gradient(
    low = "red", high = "white",
    name = "P-Value",
    na.value = "gray",
    trans = "reverse"  # Reverse scale so red represents more significant p-values
  ) +
  labs(
    title = "Pairwise Procrustes ANOVA P-Value Heatmap",
    x = "Tribe 1",
    y = "Tribe 2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  )


# Add a numeric significance level for mapping to colors
heatmap_data <- heatmap_data %>%
  mutate(SignificanceLevel = case_when(
    Significance == "***" ~ 3,
    Significance == "**" ~ 2,
    Significance == "*" ~ 1,
    TRUE ~ 0  # No stars
  ))

# Plot the heatmap
ggplot(heatmap_data, aes(x = Tribe1, y = Tribe2, fill = SignificanceLevel)) +
  geom_tile(color = "white") +  # Add white borders to tiles
  geom_text(aes(label = Significance), color = "black", size = 3) +  # Add significance stars
  scale_fill_gradient(
    low = "lightblue", high = "darkblue",
    name = "Significance",
    breaks = c(0, 1, 2, 3),
    labels = c("No significance", "*", "**", "***")
  ) +
  labs(
    title = "Pairwise Procrustes ANOVA Heatmap (Significance Level)",
    x = "Tribe 1",
    y = "Tribe 2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  )





# Add a numeric significance level for mapping to colors
heatmap_data <- heatmap_data %>%
  mutate(SignificanceLevel = case_when(
    Significance == "***" ~ 3,
    Significance == "**" ~ 2,
    Significance == "*" ~ 1,
    TRUE ~ 0  # No stars
  ))

# Plot the heatmap
ggplot(heatmap_data, aes(x = Tribe1, y = Tribe2, fill = SignificanceLevel)) +
  geom_tile(color = "white") +  # Add white borders to tiles
  geom_text(aes(label = Significance), color = "black", size = 3) +  # Add significance stars
  scale_fill_gradient(
    low = "lightblue", high = "darkblue",
    name = "P-Value Cutoff",
    breaks = c(0, 1, 2, 3),
    labels = c(
      "P > 0.05 (ns)",  # No significance
      "0.01 < P ≤ 0.05 (*)",  # Single star
      "0.001 < P ≤ 0.01 (**)",  # Double stars
      "P ≤ 0.001 (***)"  # Triple stars
    )
  ) +
  labs(
    title = "Pairwise Procrustes ANOVA Heatmap (P-Value Cutoffs)",
    x = "Tribe 1",
    y = "Tribe 2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  )

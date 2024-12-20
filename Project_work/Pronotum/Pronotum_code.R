#Organizing the data

# Function to add a file ID with global entry order to each TPS entry
add_file_id <- function(tps_data, file_name, start_counter) {
  modified_data <- c()
  
  # Initialize an entry buffer to store each block of data until ID line
  entry_buffer <- c()
  
  # Initialize a local counter for the current file
  local_counter <- start_counter
  
  for (line in tps_data) {
    entry_buffer <- c(entry_buffer, line)
    
    # If we reach an "ID=" line, we assume the entry is complete
    if (grepl("^ID=", line)) {
      # Create a new ID combining the original ID and the global order
      new_id <- paste(gsub("^ID=", "", line), local_counter, sep = "#")
      
      # Replace the original ID line with the new ID
      entry_buffer <- gsub("^ID=.*", paste("ID=", new_id, sep = ""), entry_buffer)
      
      # Insert the file ID after the ID line
      entry_buffer <- c(entry_buffer, paste("Original_File_ID=", file_name, sep = ""))
      
      # Append the current entry to the modified data
      modified_data <- c(modified_data, entry_buffer)
      
      # Clear the buffer for the next entry
      entry_buffer <- c()
      
      # Increment the local counter
      local_counter <- local_counter + 1
    }
  }
  
  # Return the modified data and the new global entry count
  return(list(modified_data = modified_data, new_count = local_counter))
}

# List of TPS files
tps_files <- c("Aesalini.TPS", "Ceratognathini.TPS", "Chiasognathini.TPS", 
               "Lamprimini.TPS", "Lucanini.TPS", "Nicagini.TPS", "Platycerini.TPS",
               "Platyceroidini.TPS", "Streptocerini.TPS", "Syndesinae.TPS")

# Initialize an empty vector to store the modified data from all files
combined_tps_data <- c()

# Initialize a global counter for the rank in the combined TPS file
global_counter <- 1

# Loop through each TPS file, read, modify, and append the modified data
for (tps_file in tps_files) {
  # Read the TPS file
  tps_data <- readLines(tps_file)
  
  # Extract the base file name (without ".TPS") to use as the file ID
  file_name <- gsub("\\.TPS$", "", tps_file)
  
  # Modify the TPS data by adding the file ID and global entry order
  result <- add_file_id(tps_data, file_name, global_counter)
  
  # Update the combined data
  combined_tps_data <- c(combined_tps_data, result$modified_data)
  
  # Update the global counter to the new value after processing the current file
  global_counter <- result$new_count
}

# Write the combined modified data to a new TPS file
writeLines(combined_tps_data, "Pronotum.TPS")

cat("Modified TPS data saved to: Pronotum.TPS")



###############################
# Load necessary packages
library(geomorph)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(viridis)
library(plotly)

tps_file <- "Pronotum.TPS"
landmark_data <- readland.tps(tps_file, specID = "ID")  # Assuming "ID" contains the specimen IDs

# Extract the file of origin from the TPS file
lines <- readLines(tps_file)
file_origin <- gsub("^Original_File_ID=", "", lines[grep("^Original_File_ID=", lines)])  # Extract file origin

# Create a new variable for Tribe (just the file origin)
tribe <- factor(file_origin)  # Convert to factor for grouping

# Perform Procrustes alignment (GPA)
gpa_result <- gpagen(landmark_data)

# Perform PCA on the aligned shapes using geomorph's gm.prcomp
pca_result <- gm.prcomp(gpa_result$coords)

# Extract the variance explained by PC1 and PC2
explained_variance <- 100 * (pca_result$sdev^2 / sum(pca_result$sdev^2))

# Create informative axis labels
pc1_label <- paste0("PC1 (", round(explained_variance[1], 2), "% variance)")
pc2_label <- paste0("PC2 (", round(explained_variance[2], 2), "% variance)")

# Extract PCA scores and prepare data for plotting
pca_scores <- as.data.frame(pca_result$x)  # PCA scores from gm.prcomp
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))  # Ensure column names are PC1, PC2, etc.
pca_scores$Tribe <- tribe  # Add Tribe as a factor

# Create convex hulls for each tribe
hulls <- pca_scores %>%
  group_by(Tribe) %>%
  filter(n() > 2) %>%  # Only keep tribes with more than two points
  slice(chull(PC1, PC2))  # Find the convex hull points

# PCA plot with ggplot
pca_plot <- ggplot() +
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, fill = Tribe, group = Tribe), alpha = 0.2, show.legend = FALSE) + # Draw hulls first
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = Tribe, text = rownames(pca_scores)), size = 0.5) + # Then draw points on top
  labs(title = "Shape of Stag Beetle Pronota",
       x = pc1_label,
       y = pc2_label) +
  theme_minimal() +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(title = "Tribe", override.aes = list(size = 3)))  # Adjust title only

# Convert to interactive plot
interactive_plot <- ggplotly(pca_plot, tooltip = "text") 

# Hide hoverinfo for hulls by adjusting the traces after ggplotly
for (i in seq_along(interactive_plot$x$data)) {
  if (interactive_plot$x$data[[i]]$type == "scatter" && interactive_plot$x$data[[i]]$mode == "lines") {
    interactive_plot$x$data[[i]]$hoverinfo <- "none"  # Remove tooltips for hulls
  }
}

# Show the interactive plot
interactive_plot

#################################################
#Plotting individual sample shapes 

# Function to extract landmarks for a specific sample ID from the TPS file
extract_landmarks <- function(tps_file, sample_id) {
  # Read the TPS file
  tps_data <- readLines(tps_file)
  
  # Initialize variables
  sample_landmarks <- NULL
  current_entry <- NULL
  is_target_sample <- FALSE
  
  for (line in tps_data) {
    # Check if the line starts with "ID="
    if (grepl(paste0("^ID=", sample_id), line)) {
      is_target_sample <- TRUE  # Start capturing landmarks for this sample
      current_entry <- c(line)  # Include the ID line
    } else if (is_target_sample) {
      if (grepl("^LM=", line) || grepl("^[0-9]", line)) {
        # If it's a landmark line, append it to the current entry
        current_entry <- c(current_entry, line)
      } else if (grepl("^ID=", line) && !grepl(paste0("^ID=", sample_id), line)) {
        # If we reach a new ID and it's not our target sample, stop capturing
        break
      }
    }
  }
  
  # If we captured landmarks, create a matrix
  if (!is.null(current_entry) && length(current_entry) > 1) {
    sample_landmarks <- matrix(as.numeric(unlist(strsplit(current_entry[-1], " "))), ncol = 2, byrow = TRUE)
  }
  
  return(sample_landmarks)
}

# Specify the TPS file and the sample ID you want to plot
tps_file <- "Pronotum.TPS"
sample_id <- "1223"  # Change this to the ID of the sample you want to plot

# Extract landmarks for the specific sample
sample_landmarks <- extract_landmarks(tps_file, sample_id)

# Check if landmarks were successfully extracted
if (is.null(sample_landmarks)) {
  stop("No landmarks found for the specified sample ID.")
}

# Create a data frame for plotting
landmark_df <- data.frame(x = sample_landmarks[, 1], y = sample_landmarks[, 2])

# Plot the shape using ggplot2
ggplot(landmark_df, aes(x = x, y = y)) +
  geom_point(size = 1) +  # Plot the landmarks
  geom_path(size = 1, color = "blue") +  # Connect the landmarks
  labs(title = paste("Shape of Sample ID:", sample_id),
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()

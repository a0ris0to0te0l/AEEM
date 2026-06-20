# Clear wortkspace
rm(list = ls())

# Clear plots
check_dev <- dev.list()
if(!is.null(check_dev)){
  dev.off(dev.list()["RStudioGD"])  
}

# Get location of current script
fileloc <- dirname(rstudioapi::getSourceEditorContext()$path)

# Set working directory to script location
setwd(fileloc)

# Remove fileloc variable
rm(fileloc, check_dev)

library(tidyverse)
library(moments)
library(scales)

#Import the base .csv file
base_data <- read_csv("7.4-Sofia properties.csv")

#Import the fixed encoding .csv 
base_fix_enconding <- read_csv("Sofia_properties_encoding.csv")

#Import the fixed encoding .csv
median_per_quarted <- read_csv("median_per_quartid.csv")

#type Descriptives
# Create a frequency table for property types
tbl_type <- table(base_fix_enconding$type)

# Convert the table to a data frame
type_df <- as.data.frame(tbl_type)

# Rename columns for clarity
colnames(type_df) <- c("Property_Type", "Count")

# Calculate percentages
type_df$Percentage <- (type_df$Count / sum(type_df$Count)) * 100

# Create the histogram for property types with percentages
type_hist <- ggplot(type_df, aes(x = Property_Type, y = Percentage)) +
  geom_bar(stat = "identity", color = "lightgray", fill = "black", alpha = 0.5) +
  labs(title = "Distribution of Property Types in Sofia", 
       x = "Property Type", 
       y = "Percentage (%)") +  # Updated y-label
  theme_minimal() +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            vjust = -0.5,  # Adjusts the position of the text above the bars
            color = "black")  # Text color

# Save the histogram to a PDF file
pdf("property_type_histogram.pdf", width = 8, height = 6)  # Open the PDF device

# Display the histogram
print(type_hist)

# Close the PDF device
dev.off()

#quartid Descriptives
# Step 1: Create frequency table for quartid
tbl_quartid <- table(base_fix_enconding$quartid)
quartid_df <- as.data.frame(tbl_quartid)

# Step 2: Rename columns
colnames(quartid_df) <- c("Quartid", "Count")

# Step 3: Calculate total count and threshold for "Other"
total_count <- sum(quartid_df$Count)
threshold <- total_count * 0.01  # 1%

# Step 4: Group low-count quartids into "Other"
quartid_df$Quartid <- ifelse(quartid_df$Count < threshold,
                             "Other", 
                             as.character(quartid_df$Quartid))

# Step 5: Summarize counts for the new grouping
final_quartid_df <- aggregate(Count ~ Quartid, data = quartid_df, sum)

# Step 6: Sort the data frame by Count in descending order
final_quartid_df <- final_quartid_df[order(-final_quartid_df$Count), ]

# Step 7: Calculate percentages
final_quartid_df$Percentage <- (final_quartid_df$Count / total_count) * 100

# Step 8: Create a label for the legend
final_quartid_df$Label <- paste(final_quartid_df$Quartid, " (", round(final_quartid_df$Percentage, 1), "%)", sep = "")

# Step 9: Create a color scale
colors <- colorRampPalette(c("black", "white"))(nrow(final_quartid_df))

# Step 10: Create a pie chart with sorted colors and labels in the legend
quartid_pie_chart <- ggplot(final_quartid_df, aes(x = "", y = Count, fill = factor(Label, levels = final_quartid_df$Label))) + 
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colors) +  
  labs(title = "Distribution of Quartid in Sofia") +
  theme_void() +  
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1))  # Adjusting the title position

# Save the histogram to a PDF file
pdf("quartid_pie_chart.pdf",width = 8, height = 6)  # Open the PDF device

# Display the pie chart
print(quartid_pie_chart)

# Close the PDF device
dev.off()

#price Descriptives
# Create frequency table for price
tbl_price <- table(base_fix_enconding$`price (Euro)`)

# Convert the table to a data frame
price_df <- as.data.frame(tbl_price)
colnames(price_df) <- c("Price", "Count")

# Convert Price to numeric for proper binning (if it's not already)
price_df$Price <- as.numeric(as.character(price_df$Price))  # Ensure Price is numeric

# Create the histogram with custom bin width
price_hist <- ggplot(base_fix_enconding, aes(x = `price (Euro)`)) +
  geom_histogram(aes(y = ..count../sum(..count..) * 100), 
                 binwidth = 25000, color = "lightgray", fill = "black", alpha = 0.5) +  # Set bin width to 50,000
  labs(title = "Distribution of Property Prices in Sofia", 
       x = "Price (Euro)", 
       y = "Percentage (%)") +  # Updated y-label
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(base_fix_enconding$`price (Euro)`), by = 25000), labels = scales::comma)  # Absolute values & increase by 50,000

# Save the histogram to a PDF file
pdf("price_histogram.pdf",width = 8, height = 6)  # Open the PDF device

# Display the price histogram
print(price_hist)

# Close the PDF device
dev.off()

summary(price_df)
price_max <- max(base_fix_enconding$`price (Euro)`)
price_min <- min(base_fix_enconding$`price (Euro)`)
price_mean <- mean(base_fix_enconding$`price (Euro)`)
price_mad <- mad(base_fix_enconding$`price (Euro)`)
price_sd <- sd(base_fix_enconding$`price (Euro)`)
price_skewness <- skewness(base_fix_enconding$`price (Euro)`)
price_kurtosis <- kurtosis(base_fix_enconding$`price (Euro)`)

#area Descriptives
# Create frequency table for area
tbl_area <- table(base_fix_enconding$area)

# Convert the table to a data frame
area_df <- as.data.frame(tbl_area)
colnames(area_df) <- c("Area", "Count")

# Convert Area to numeric for proper binning (if it's not already)
area_df$Area <- as.numeric(as.character(area_df$Area))  # Ensure Area is numeric

# Create the histogram with custom bin width
area_hist <- ggplot(base_fix_enconding, aes(x = area)) +
  geom_histogram(aes(y = ..count../sum(..count..) * 100), 
                 binwidth = 25, color = "lightgray", fill = "black", alpha = 0.5) +  # Set bin width to 10
  labs(title = "Distribution of Property Areas in Sofia", 
       x = "Area (sq m)", 
       y = "Percentage (%)") +  # Updated y-label to reflect percentages
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(base_fix_enconding$area, na.rm = TRUE), by = 25), 
                     labels = scales::comma)  # Absolute values & increase by 10

# Save the histogram to a PDF file
pdf("area_histogram.pdf",width = 8, height = 6)  # Open the PDF device

# Display the area histogram
print(area_hist)

# Close the PDF device
dev.off()

summary(area_df)
area_max <- max(base_fix_enconding$area)
area_min <- min(base_fix_enconding$area)
area_mean <- mean(base_fix_enconding$area)
area_mad <- mad(base_fix_enconding$area)
area_sd <- sd(base_fix_enconding$area)
area_skewness <- skewness(base_fix_enconding$area)
area_kurtosis <- kurtosis(base_fix_enconding$area)


#sq_m Descriptives
# Create frequency table for euro/sq.m.
tbl_sq_m <- table(base_fix_enconding$`euro/sq.m.`)

# Convert the table to a data frame
sq_m_df <- as.data.frame(tbl_sq_m)
colnames(sq_m_df) <- c("Euro_per_sq_m", "Count")

# Convert Euro_per_sq_m to numeric for proper binning (if it's not already)
sq_m_df$Euro_per_sq_m <- as.numeric(as.character(sq_m_df$Euro_per_sq_m))  # Ensure values are numeric

# Create the histogram with a bin width of 100
sq_m_hist <- ggplot(base_fix_enconding, aes(x = `euro/sq.m.`)) +
  geom_histogram(aes(y = ..count../sum(..count..) * 100), 
                 binwidth = 250, color = "lightgray", fill = "black", alpha = 0.5) +  # Set bin width to 100
  labs(title = "Distribution of Property Prices per Square Meter in Sofia", 
       x = "Price per sq m (Euro)", 
       y = "Percentage (%)") +  # Updated y-label to reflect percentages
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(base_fix_enconding$`euro/sq.m.`, na.rm = TRUE), by = 250), 
                     labels = scales::comma)  # Absolute values & increase by 100

# Save the histogram to a PDF file
pdf("sq_m_histogram.pdf",width = 8, height = 6)  # Open the PDF device

# Display the price per square meter histogram
print(sq_m_hist)

# Close the PDF device
dev.off()

summary(sq_m_df)
sq_m_max <- max(base_fix_enconding$`euro/sq.m.`)
sq_m_min <- min(base_fix_enconding$`euro/sq.m.`)
sq_m_mean <- mean(base_fix_enconding$`euro/sq.m.`)
sq_m_mad <- mad(base_fix_enconding$`euro/sq.m.`)
sq_m_sd <- sd(base_fix_enconding$`euro/sq.m.`)
sq_m_skewness <- skewness(base_fix_enconding$`euro/sq.m.`)
sq_m_kurtosis <- kurtosis(base_fix_enconding$`euro/sq.m.`)

#constr. Descriptives
# Step 1: Create frequency table for construction types
tbl_constr <- table(base_fix_enconding$constr.)
constr_df <- as.data.frame(tbl_constr)
colnames(constr_df) <- c("Construction_Type", "Count")

# Step 2: Calculate percentages
total_count <- sum(constr_df$Count)
constr_df$Percentage <- (constr_df$Count / total_count) * 100

# Step 3: Create the bar chart using percentages
constr_chart <- ggplot(constr_df, aes(x = Construction_Type, y = Percentage)) +
  geom_bar(stat = "identity", color = "lightgray", fill = "black", alpha = 0.5) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            vjust = -0.5,  # Adjust the vertical position of the text
            color = "black") +  # Text color
  labs(title = "Distribution of Property Construction Types in Sofia", 
       x = "Construction Type", 
       y = "Percentage") +
  theme_minimal()

# Save the histogram to a PDF file
pdf("const_histogram.pdf",width = 8, height = 6)  # Open the PDF device

# Step 4: Display the construction types chart
print(constr_chart)

# Close the PDF device
dev.off()


#year Descriptives
# Step 1: Ensure the year column is numeric for proper binning
base_fix_enconding$year <- as.numeric(as.character(base_fix_enconding$year))

# Step 1.1: Change any 0 values to NA
base_fix_enconding$year[base_fix_enconding$year == 0] <- NA

year_df <- na.omit(as.data.frame(base_fix_enconding$year_filtered))

# Step 2: Remove rows with NA values in the year column (if any)
base_fix_enconding <- na.omit(base_fix_enconding)

# Step 3: Create custom bin ranges by decade
breaks <- seq(1860, ceiling(max(base_fix_enconding$year)), by = 10)

# Create the histogram data
year_hist <- ggplot(base_fix_enconding, aes(x = year)) +
  geom_histogram(aes(y = ..count../sum(..count..)), 
                 breaks = breaks, color = "lightgray", alpha = 0.5) +  # Specify custom breaks
  labs(title = "Distribution of Property Years in Sofia", 
       x = "Year", 
       y = "Percentage") +  
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(breaks = breaks) +  # Use the same breaks for x-axis labels
  theme_minimal()

# Save the histogram to a PDF file
pdf("year_histogram.pdf",width = 8, height = 6)  # Open the PDF device

# Step 4: Display the histogram
print(year_hist)

# Close the PDF device
dev.off()


summary(year_df)
year_max <- max(as.numeric(base_fix_enconding$year_filtered),na.rm = T)
year_min <- min(base_fix_enconding$year_filtered,na.rm = T)
year_mean <- mean(base_fix_enconding$year_filtered,na.rm = T)
year_mad <- mad(base_fix_enconding$year_filtered,na.rm = T)
year_sd <- sd(base_fix_enconding$year_filtered,na.rm = T)
year_skewness <- skewness(base_fix_enconding$year_filtered,na.rm = T)
year_kurtosis <- kurtosis(base_fix_enconding$year_filtered,na.rm = T)

count_non_na_year <- sum(!is.na(as.numeric(year_df$`base_fix_enconding$year`)))

#floor Descriptives
# Create a frequency table for floor levels
tbl_floor <- table(na.omit(base_fix_enconding$floor_filtered))

# Convert the table to a data frame
floor_df <- as.data.frame(tbl_floor)

# Rename columns for clarity
colnames(floor_df) <- c("Floor", "Count")

# Calculate percentages
floor_df$Percentage <- (floor_df$Count / sum(floor_df$Count)) * 100

# Create the histogram for floors
floor_hist <- ggplot(floor_df, aes(x = Floor, y = Percentage)) +
  geom_bar(stat = "identity", color = "lightgray", fill = "black", alpha = 0.5) +
  labs(title = "Distribution of Floor Levels in Sofia", 
       x = "Floor Level", 
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal()

# Save the histogram to a PDF file
pdf("floor_histogram.pdf", width = 8, height = 6)  # Open the PDF device

# Display the histogram
print(floor_hist)

# Close the PDF device
dev.off()

summary(base_fix_enconding$floor_filtered,na.rm=T)
floor_max <- max(base_fix_enconding$floor_filtered,na.rm=T)
floor_min <- min(base_fix_enconding$floor_filtered,na.rm=T)
floor_mean <- mean(base_fix_enconding$floor_filtered,na.rm=T)
floor_mad <- mad(base_fix_enconding$floor_filtered,na.rm=T)
floor_sd <- sd(base_fix_enconding$floor_filtered,na.rm=T)
floor_skewness <- skewness(base_fix_enconding$floor_filtered,na.rm=T)
floor_kurtosis <- kurtosis(base_fix_enconding$floor_filtered,na.rm=T)



  #heating Descriptives
  # Create a frequency table for heating types
  tbl_heating <- table(base_fix_enconding$TEC)
  
  # Convert the table to a data frame
  heating_df <- as.data.frame(tbl_heating)
  
  # Ensure there are exactly two levels and rename columns for clarity
  if (ncol(heating_df) == 2) {
    colnames(heating_df) <- c("Heating_Type", "Count")
  } else {
    stop("Unexpected number of columns in heating_df. Please check the input data.")
  }
  
  # Calculate percentages for the labels
  heating_df$Percentage <- (heating_df$Count / sum(heating_df$Count)) * 100
  
  # Create the pie chart
  heating_pie_chart <- ggplot(heating_df, aes(x = "", y = Count, fill = Heating_Type)) +
    geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.5) +  # Change border color to black
    coord_polar(theta = "y") +  # Convert to pie chart
    scale_fill_manual(values = c("YES" = "black", "NO" = "black")) +  # Set both colors to black
    labs(title = "Distribution of Heating Types in Sofia", fill = "Heating Type") +
    theme_void() +  # Use a void theme for a cleaner look
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
              position = position_stack(vjust = 0.5), 
              color = "white")  # Add percentage labels
  
  
  # Save the histogram to a PDF file
  pdf("heating_pie_chart.pdf", width = 8, height = 6)  # Open the PDF device
  
  # Display the pie chart
  print(heating_pie_chart)
  
  # Close the PDF device
  dev.off()
  
  #seller Descriptives
  #Create a frequency table for sellers
  tbl_seller <- table(base_fix_enconding$seller)
  
  # Convert the table to a data frame
  seller_df <- as.data.frame(tbl_seller)
  
  # Rename columns for clarity
  colnames(seller_df) <- c("Seller_Type", "Count")
  
  # Calculate percentages for the labels
  seller_df$Percentage <- (seller_df$Count / sum(seller_df$Count)) * 100
  
  # Create the pie chart
  seller_pie_chart <- ggplot(seller_df, aes(x = "", y = Count, fill = Seller_Type)) +
    geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.5) +  # Border color
    coord_polar(theta = "y") +  # Convert to pie chart
    scale_fill_manual(values = c("Agency" = "black", "Private" = "lightgray")) +  # Define colors for each Seller_Type
    labs(title = "Distribution of Sellers in Sofia", fill = "Seller Type") +  # Title and fill label
    theme_void() +  # Clean look
    theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 8)) +  # Legend adjustments
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
              position = position_stack(vjust = 0.5), 
              color = "white")  # Add percentage labels
  
  # Save the histogram to a PDF file
  pdf("seller_pie_chart.pdf", width = 8, height = 6)  # Open the PDF device
  
  # Display the pie chart
  print(seller_pie_chart)
  
  # Close the PDF device
  dev.off()
  
  
  
  # Create a frequency table for the property types based on the Difference from Median
  tbl_diff_median <- table(base_fix_enconding$Diff_median)
  
  # Convert the frequency table to a data frame for easier manipulation
  diff_median_df <- as.data.frame(tbl_diff_median)
  
  # Rename the columns for clarity: first column is 'Difference_from_Median', second is 'Count'
  colnames(diff_median_df) <- c("Difference_from_Median", "Count")
  
  # Calculate the percentage of each count relative to the total
  diff_median_df$Percentage <- (diff_median_df$Count / sum(diff_median_df$Count)) * 100
  
  # Sort the data frame by 'Difference_from_Median' in ascending order for better visualization
  diff_median_df <- diff_median_df[order(diff_median_df$Difference_from_Median), ]
  
  # Create a histogram to visualize the distribution of premium values
  diff_median_hist <- ggplot(base_fix_enconding, aes(x = Diff_median)) +
    geom_histogram(aes(y = (..count.. / sum(..count..)) * 100),  # Calculate percentage for each bin
                   bins = 30, color = "lightgray", fill = "black", alpha = 0.5) +
    labs(title = "Distribution of Premium", 
         x = "Premium", 
         y = "Percentage (%)") +  # Set the y-axis label to Percentage
    scale_x_continuous(breaks = seq(-1200, max(base_fix_enconding$Diff_median, na.rm = TRUE), by = 250)) +  # Set x-axis labels with rounded scale
    theme_minimal()  # Use a minimal theme for better aesthetics
  
  # Save the histogram to a PDF file for later viewing
  pdf("difference_from_median_histogram.pdf", width = 8, height = 6)  # Open the PDF device
  
  # Display the histogram
  print(diff_median_hist)
  
  # Close the PDF device to finalize the file
  dev.off()  # Close the PDF device
  
  diff_max <- max(as.numeric(base_fix_enconding$Diff_median), na.rm = TRUE)
  diff_min <- min(as.numeric(base_fix_enconding$Diff_median), na.rm = TRUE)
  diff_mean <- mean(as.numeric(base_fix_enconding$Diff_median), na.rm = TRUE)
  diff_mad <- mad(as.numeric(base_fix_enconding$Diff_median), na.rm = TRUE)
  diff_sd <- sd(as.numeric(base_fix_enconding$Diff_median), na.rm = TRUE)
  diff_skewness <- skewness(as.numeric(base_fix_enconding$Diff_median), na.rm = TRUE)
  diff_kurtosis <- kurtosis(as.numeric(base_fix_enconding$Diff_median), na.rm = TRUE)
  
  # Create the percentage histogram
  median_sq_df_hist <- ggplot(median_per_quarted, aes(x = MEDIAN_SQ_M)) +
    geom_histogram(aes(y = (..count.. / sum(..count..)) * 100),  # Calculate percentage
                   bins = 20, color = "lightgray", fill = "black", alpha = 0.5) +
    labs(title = "Percentage Distribution of Median Square Prices of Neighborhood Medians", 
         x = "Median Square Price (Euro)", 
         y = "Percentage (%)") +
    theme_minimal() +
    scale_x_continuous(
      breaks = seq(0, max(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE), by = 100),  # Adjust increments as needed
      labels = function(x) scales::comma(x, accuracy = 1)  # Format labels with commas
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
      plot.title = element_text(hjust = 0.5)  # Center title for added aesthetic appeal
    )
  
  # Save the histogram to a PDF file
  pdf("median_percentage_histogram.pdf", width = 8, height = 6)  # Open the PDF device
  
  # Display the histogram
  print(median_sq_df_hist)
  
  # Close the PDF device
  dev.off()
  
  median_max <- max(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE)
  median_min <- min(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE)
  median_mean <- mean(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE)
  median_mad <- mad(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE)
  median_sd <- sd(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE)
  median_skewness <- skewness(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE)
  median_kurtosis <- kurtosis(median_per_quarted$MEDIAN_SQ_M, na.rm = TRUE)
  
  
  
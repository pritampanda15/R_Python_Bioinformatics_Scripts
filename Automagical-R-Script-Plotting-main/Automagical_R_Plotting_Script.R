###########################################################################################################
###########################################################################################################
################################                                       ####################################
################################   The Automagical R Plotting Script   ####################################
################################                                       ####################################
################################          Bioinformatics Copilot       ####################################
################################                         	             ####################################
################################                                       ####################################
###########################################################################################################
###########################################################################################################


# This script automatically creates a bunch of relevant plots based on the type of data you enter
# First it automatically detects your file's delimiter (tab or comma) and whether it has a header
# Then it categorizes each column by the type of data contained in the column
# And finally it makes a bunch of plots fitting each column and combinations of different columns


###############    Step 1:  Change this filename    ############### 
filename <- "~/Automagical/variants.bed"

  # other examples:
  # filename <- "~/Desktop/my_folder/myfile.txt"
  # filename <- "~/Desktop/my_folder/my_other_file.csv"
  
  # Files can be tab or comma separated, and can have a header or not, though headers help to label the plots correctly
  
###############    Step 2: Hold down Command and Shift together while you press Enter. This will run the whole script    ############### 



# Define styles and plot sizes (optional)

png_width = 2000
png_height = 1000
png_res = 300


# Load necessary packages
library(ggplot2)


# Load a file, check for the best delimiter to use, and determine if it has a header
# Note: This is a really useful function you can include in other scripts to load any .csv or tab-delimited file into a data.frame
load_any_file <- function(filename) {
    data <- read.csv(filename, sep="\t")
    sep <- "\t"
    header <- TRUE
    
    csv <- read.csv(filename)
    if (length(names(csv)) > length(names(data))) {
      data <- csv
      sep <- ","
    }
    remove(csv)
    
    no_header_data <- read.csv(filename, sep=sep, header = FALSE)
    types_without_header <- sapply(no_header_data, class)
    types_with_header <- sapply(data, class)
    if (identical(unname(types_with_header),unname(types_without_header))) {
      header <- FALSE
      data <- read.csv(filename, sep=sep, header = header)
    }
    remove(no_header_data)
    return(data)
}


data <- load_any_file(filename)

# For outputting messages to the user (that's you!) neatly
pretty_print <- function(my_string) {
  cat(paste(my_string, "\n"))
}

pretty_print("If you experience issues with the plots, check the table below to make sure the top of the file got imported correctly:")
head(data)

# Inspect the types of each column
pretty_print("These are the types that each column is getting interpreted as by R itself:")
sapply(data, class)

# Set column identities to classes as default
column_names <- names(data)
column_identities <- sapply(data,class)

# This function returns the set of input strings that match known human chromosome names
valid_chromosome_names <- function(chrom_names) {
  possible_chrom_names <- c(seq(1,40), paste("chr",seq(1,40)), "X","chrX","Y","chrY","M","MT","chrM","chrMT")
  return(intersect(chrom_names, possible_chrom_names))
}

# This function replaces NA's with a prettier default value
rename_NAs <- function(column,replacement_value) {
  column = factor(column, levels=c(levels(column), replacement_value))
  column[is.na(column)] <- replacement_value
  return(column)
}

# Replace bad header characters like spaces with underscores for better file names
underscorify <- function(my_string) {
  return(gsub(" ", "_", my_string))
}

# Loop through each column, clean up data, and assign column identities
for (column_index in seq(1,length(column_names))) {
  pretty_print(column_names[column_index])
  
  column_name <- column_names[column_index]
  
  if (class(data[,column_name]) == "factor") {
    # Check if this factor may be chromosome names, if so filter to main chromosomes
    if (length(unique(data[,column_name])) == 1) {
      column_identities[column_index] <- "same"
    } else if (max(table(data[,column_name])) == 1 && length(table(data[,column_name])) > 100) {
      column_identities[column_index] <- "ID"
    } else {
      categories <- table(data[,column_name])
      matching_chrom_names <- valid_chromosome_names(names(categories))
      if (sum(categories[names(categories) %in% matching_chrom_names])/nrow(data) > 0.90) {
        # Apply natural sort to factor data to deal with possible chromosome names and standardize ordering of categories
        library(naturalsort)
        column_identities[column_index] <- "chromosome"
        data[,column_name] <- factor(data[,column_name], levels = naturalsort(matching_chrom_names))
        data[,column_name] <- rename_NAs(data[,column_name], "other")
      } else if (length(table(data[,column_name])) <= 10 ) {
        column_identities[column_index] <- "category"
      }
    }
  } else if (class(data[,column_name]) == "integer") {
    column_identities[column_index] <- "number"
    
    # This prevents position columns from being the secondary variable in scatterplots:
    position_indicators <- c("pos","start","stop","end","position","coord")
    if (length(grep(paste(position_indicators,collapse="|"), column_name, value=FALSE)) > 0) {
      column_identities[column_index] <- "position"
    }
    # /end of position classification
  }
  pretty_print(column_identities[column_index])
}

categorical_variables <- column_names[column_identities == "category"]
numerical_variables <- column_names[column_identities == "number"]


theme_set(theme_gray()+theme(
  strip.text.y = element_blank() # remove labels from facet_grid since they often don't fit and we use colors shown in the legend instead
  ))

# Here you can replace the color palette if you want:
large_palette <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#ff9896", "#c5b0d5", "#8c564b", "#e377c2", "#bcbd22", "#9edae5", "#c7c7c7", "#d62728", "#ffbb78", "#98df8a", "#ff7f0e", "#f7b6d2", "#c49c94", "#dbdb8d", "#aec7e8", "#17becf", "#2ca02c", "#7f7f7f", "#1f77b4", "#9467bd")


for (column_index in seq(1,length(column_names))) {
  column_name <- column_names[column_index]
  pretty_print(paste("Plotting",column_name))
  
  if (column_identities[column_index] == "ID") {
    # don't plot
    pretty_print("    Ignoring factor where all rows have different values")
  } else if (column_identities[column_index] == "same") {
    # don't plot
    pretty_print(paste("    Ignoring factor where all rows have the same value:", data[1,column_name]))
  } else if (column_identities[column_index] == "factor" || column_identities[column_index] == "category" || column_identities[column_index] == "chromosome") {
    
    ############# CATEGORICAL ##############
    # Plot counts split on x-axis by this categorical variable
    if (column_identities[column_index] == "chromosome") {
      # Don't add colors if there are many types
      png(file=paste(filename,".plots.individual_variables.",underscorify(column_name),".bar_chart",".png", sep=''),width=png_width,height=png_height,res=png_res)
      print(ggplot(data, aes_string(x=column_name)) + geom_bar())
      dev.off()
      
      png(file=paste(filename,".plots.individual_variables.",underscorify(column_name),".bar_chart.flipped",".png", sep=''),width=png_width,height=png_height,res=png_res)
      print(ggplot(data, aes_string(x=column_name)) + geom_bar() + coord_flip() + scale_x_discrete(limits = rev(levels(data[,column_name]))))
      dev.off()
    } else {
      png(file=paste(filename,".plots.individual_variables.",underscorify(column_name),".bar_chart",".png", sep=''),width=png_width,height=png_height,res=png_res)
      print(ggplot(data, aes_string(x=column_name, fill=column_name)) + geom_bar() + scale_fill_manual(values=large_palette, drop=FALSE))
      dev.off()
      
      png(file=paste(filename,".plots.individual_variables.",underscorify(column_name),".bar_chart.flipped",".png", sep=''),width=png_width,height=png_height,res=png_res)
      print(ggplot(data, aes_string(x=column_name, fill=column_name)) + geom_bar() + scale_fill_manual(values=large_palette, drop=FALSE) + coord_flip() + scale_x_discrete(limits = rev(levels(data[,column_name]))))
      dev.off()
    }
    
    # Split by each of the other categorical variables
    for (fill_variable in categorical_variables) {
      if (fill_variable != column_name) {
        pretty_print(paste("-> splitting by",fill_variable))
        
        if (column_identities[column_index] == "chromosome" || length(table(data[,column_name]))>6) {
          # faceted
          png(file=paste(filename,".plots.combinations.",underscorify(column_name),".split_by.",underscorify(fill_variable),".bar_chart",".faceted",".png", sep=''),width=png_width,height=png_height,res=png_res)
          print(ggplot(data, aes_string(x=column_name, fill=fill_variable)) + geom_bar() + facet_grid(reformulate(".",fill_variable)) + scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values=large_palette, drop=FALSE))
          dev.off()
          # stacked
          png(file=paste(filename,".plots.combinations.",underscorify(column_name),".split_by.",underscorify(fill_variable),".bar_chart",".faceted",".png", sep=''),width=png_width,height=png_height,res=png_res)
          print(ggplot(data, aes_string(x=column_name, fill=fill_variable)) + geom_bar() + scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values=large_palette, drop=FALSE))
          dev.off()
        } else {
          # faceted
          png(file=paste(filename,".plots.combinations.",underscorify(column_name),".split_by.",underscorify(fill_variable),".bar_chart",".faceted",".png", sep=''),width=png_width,height=png_height,res=png_res)
          print(ggplot(data, aes_string(x=column_name, fill=fill_variable)) + geom_bar() + facet_grid(reformulate(".",fill_variable)) + coord_flip() + scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values=large_palette, drop=FALSE) + scale_x_discrete(limits = rev(levels(data[,column_name]))))
          dev.off()
        }
      }
    }
  } else if (column_identities[column_index] %in% c("number","position")) {
    binwidth <- (max(data[,column_name])-min(data[,column_name]))/50
    if (class(data[,column_name]) == "integer") {
      if (binwidth < 1) {
        binwidth <- 1
      }
    }
    
    png(file=paste(filename,".plots.individual_variables.",underscorify(column_name),".histogram",".png", sep=''),width=png_width,height=png_height,res=png_res)
    print(ggplot(data, aes_string(x=column_name)) + geom_histogram(binwidth=binwidth))
    dev.off()
    
    # Split by each of the other categorical variables
    for (fill_variable in categorical_variables) {
      pretty_print(paste("-> splitting by",fill_variable))
      
      # faceted
      png(file=paste(filename,".plots.combinations.",underscorify(column_name),".split_by.",underscorify(fill_variable),".histogram",".png", sep=''),width=png_width,height=png_height,res=png_res)
      print(ggplot(data, aes_string(x=column_name, fill=fill_variable))  + geom_histogram(binwidth=binwidth) + facet_grid(reformulate(".",fill_variable)) + scale_fill_manual(values=large_palette, drop=FALSE))
      dev.off()
      
      # boxplot
      png(file=paste(filename,".plots.combinations.",underscorify(column_name),".split_by.",underscorify(fill_variable),".boxplot",".png", sep=''),width=png_width,height=png_height,res=png_res)
      print(ggplot(data, aes_string(y=column_name, x=fill_variable,fill=fill_variable)) + geom_boxplot() + coord_flip() + scale_fill_manual(values=large_palette, drop=FALSE) + scale_x_discrete(limits = rev(levels(data[,column_name]))))
      dev.off()
    }
    
    # Against another numeric variable
    for (secondary_variable in numerical_variables) {
      if (secondary_variable != column_name) {
        # scatterplot
        png(file=paste(filename,".plots.combinations.",underscorify(column_name),".versus.",underscorify(secondary_variable),".scatterplot",".png", sep=''),width=png_width,height=png_height,res=png_res)
        print(ggplot(data, aes_string(x=column_name, y=secondary_variable)) + geom_point())
        dev.off()
        for (fill_variable in categorical_variables) {
          # scatterplot colored by categorical variables
          png(file=paste(filename,".plots.combinations.",underscorify(column_name),".versus.",underscorify(secondary_variable),".scatterplot",".colored_by.",underscorify(fill_variable),".png", sep=''),width=png_width,height=png_height,res=png_res)
          print(ggplot(data, aes_string(x=column_name, y=secondary_variable,color=fill_variable)) + geom_point() + scale_color_manual(values=large_palette, drop=FALSE))
          dev.off()
        }
      }
    }
  } else {
    pretty_print(paste("    Column identity not recognized:", column_identities[column_index]))
  }
}


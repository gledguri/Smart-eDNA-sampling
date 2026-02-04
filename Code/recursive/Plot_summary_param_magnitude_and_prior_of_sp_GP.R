# library(purrr)
# library(ggplot2)
# library(here)
# library(dplyr)

file <- here('Plots','Pred_maps','Round 8')

rds_files <- list.files(file, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

dfs <- imap(rds_files, function(path, i) {
	df <- readRDS(path)
	species <- basename(dirname(path))
	colnames(df) <- paste0(colnames(df), "_", i, "_", species)
	# Convert all columns to character first to avoid type conflicts
	df <- df %>% mutate(across(everything(), as.character))
	
	# Pivot longer within the function
	df_long <- df %>%
		pivot_longer(
			cols = everything(),
			names_to = "col_name",
			values_to = "value"
		) %>%
		# Extract draw number and statistic from column names
		separate(col_name, 
						 into = c("magnitude", "stat", "extra"), 
						 sep = "_", 
						 extra = "merge") %>%
		# Create parameter column
		mutate(
			param = case_when(
				stat == "param" ~ value,
				TRUE ~ NA_character_
			)
		) %>%
		# Fill down parameter names and filter
		fill(param, .direction = "down") %>%
		filter(stat != "param") %>%
		# Convert numeric columns back to numeric
		mutate(value = as.numeric(value)) %>%
		# Pivot wider to get final format
		pivot_wider(
			names_from = stat,
			values_from = value
		) %>%
		# Add species and file index info
		mutate(
			species = species,
			file_index = i
		) %>%
		select(file_index, species, magnitude, param, mean, `2.5%`, `97.5%`)
	
	return(df_long)
})
final_df <- reduce(dfs, rbind)

p1 <- final_df %>% 
	ggplot()+
	geom_point(aes(x=mean, y=param, color=species), 
						 position=position_dodge(width=0.4))+
	geom_errorbar(aes(y=param, xmin=`2.5%`, xmax=`97.5%`, color=species), 
								width=0.3, position=position_dodge(width=0.4))+
	facet_wrap(~magnitude)+
	xlim(-10,8)+
	theme_bw()
p1
ggsave(here('Plots','Pred_maps','Summary_round_8.jpg'),p1,width = 14,height = 6)

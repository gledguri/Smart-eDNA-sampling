estimate_costs <- function(
		n_collect,          # number of samples collected
		n_extract,          # number extracted
		n_qpcr_1,             # number run on qPCR
		n_qpcr_2,             # number run on qPCR
		n_meta,             # number run on metabarcoding
		boat_cost_day = 3000,
		samples_per_day = 50,
		collect_cost = 10,
		extract_cost = 5,
		qpcr_cost = 2,
		meta_cost = 4,
		eff_qpcr1 = 0.80,      # efficiency of assay 1
		eff_qpcr2 = 0.40,       # efficiency of assay 2
		eff_qpcr2_1 = 0.80     # efficiency of assay 2
) {
	
	# Boat days needed (always ceiling)
	boat_days <- ceiling(n_collect / samples_per_day)
	
	# Cost components
	boat_total   <- boat_days * boat_cost_day
	collect_total <- n_collect * collect_cost
	extract_total <- n_extract * extract_cost
	qpcr_total    <- (n_qpcr_1 + n_qpcr_2) * qpcr_cost
	meta_total    <- n_meta * meta_cost
	
	pos_qpcr1 <- eff_qpcr1 * n_qpcr_1
	pos_qpcr2 <- eff_qpcr2 * n_qpcr_2
	neg_qpcr1 <- min((n_qpcr_1-pos_qpcr1),n_qpcr_2)
	pos_qpcr2_added <- neg_qpcr1 * eff_qpcr2_1
	total_qpcr_pos <- pos_qpcr1 + pos_qpcr2_added
	
	# Total budget
	total_cost <- boat_total + collect_total + extract_total + qpcr_total + meta_total
	
	# Return a detailed breakdown
	list(
		boat_days      = boat_days,
		boat_total     = boat_total,
		collect_total  = collect_total,
		extract_total  = extract_total,
		qpcr_total     = qpcr_total,
		meta_total     = meta_total,
		total_cost     = total_cost,
		pos_qpcr1     = pos_qpcr1,
		pos_qpcr2     = pos_qpcr2,
		pos_qpcr2_added = pos_qpcr2_added,
		total_qpcr_pos     = total_qpcr_pos
	)
}



estimate_costs(
	n_collect = 1000,
	n_extract = 1000,
	n_qpcr_1  = 1000,
	n_qpcr_2  = 200,
	n_meta    = 1000,
	eff_qpcr1 = 0.80,      # efficiency of assay 1
	eff_qpcr2 = 0.40,       # efficiency of assay 2
	eff_qpcr2_1 = 0.60     # efficiency of assay 2 when assay 1 is negative
)




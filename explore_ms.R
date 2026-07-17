suppressMessages({library(macrosheds); library(dplyr)})
v <- ms_var_catalog
cat("== chem_category ==\n"); print(table(v$chem_category, useNA = "ifany"))
cat("\n== variable_code prefixes (first 3 chars) ==\n"); print(head(sort(table(substr(v$variable_code, 1, 3)), decreasing = TRUE), 12))
cat("\n== discharge-like variables ==\n")
print(v %>% filter(grepl("discharge|IS_discharge", variable_code, ignore.case = TRUE)) %>% head(3) %>% as.data.frame())
cat("\n== ms_vars_ts ==\n"); print(head(as.data.frame(ms_vars_ts), 8))
cat("\n== site_type ==\n"); print(table(ms_site_data$site_type))
cat("\n== ws_status ==\n"); print(table(ms_site_data$ws_status))

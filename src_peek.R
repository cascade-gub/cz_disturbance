suppressMessages(library(macrosheds))
print(macrosheds::ms_download_core_data)
cat("\n\n=== figshare id helpers ===\n")
print(head(macrosheds:::file_ids_for_r_package2, 20))

# Parse analytical code and repository utilities without executing analyses.
files <- unlist(lapply(c("figures", "scripts"), function(folder) {
  list.files(folder, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
}), use.names = FALSE)
stopifnot(length(files) > 0)
for (f in files) parse(f, keep.source = FALSE)
message(length(files), " R files parsed successfully.")

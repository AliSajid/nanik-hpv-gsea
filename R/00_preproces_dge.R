# 00_export_dge.R
# Reads the DESeq2 workbook, keeps every sheet with "All" in its name and writes
# each one to data/dge/<sheet>.csv (R-safe file names, column names and writing).

library(readxl)

# ---- Config ------------------------------------------------------------------
path          <- file.path("data", "nanik", "20251204_coding_total_results.xlsx")
out_dir       <- file.path("data", "dge")
sheet_pattern <- "(^|_)All($|_)"   # "All" as a whole word: HPV_Main_Effect_All, not "Small_x"

# ---- Helpers -----------------------------------------------------------------
# File-system and R-safe name: ASCII letters/digits/underscore, no leading digit
safe_name <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT", sub = "_")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- "sheet"
  ifelse(grepl("^[0-9]", x), paste0("x", x), x)
}

# Write to a temp file in the target folder, check it reads back identically in
# shape, then move it into place. A crash never leaves a half-written CSV behind.
write_csv_safe <- function(df, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = "dge_", tmpdir = dirname(file), fileext = ".tmp")
  on.exit(unlink(tmp), add = TRUE)

  write.csv(df, tmp, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

  check <- read.csv(tmp, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(check) != nrow(df) || !identical(names(check), names(df))) {
    stop("Read-back check failed for ", basename(file))
  }
  if (!file.rename(tmp, file)) stop("Could not move file into place: ", file)
  invisible(file)
}

# ---- Select sheets -----------------------------------------------------------
all_sheets <- excel_sheets(path)
sheets     <- grep(sheet_pattern, all_sheets, value = TRUE)
if (length(sheets) == 0) stop("No sheet matching '", sheet_pattern, "' in ", path)

file_names <- safe_name(sheets)
if (anyDuplicated(file_names)) {
  stop("Sheet names collide after sanitising: ",
       paste(sheets[file_names %in% file_names[duplicated(file_names)]], collapse = ", "))
}

# ---- Export ------------------------------------------------------------------
failed <- character()

for (i in seq_along(sheets)) {
  tryCatch({
    message("Sheet: ", sheets[i])

    df <- as.data.frame(read_excel(path, sheet = sheets[i], .name_repair = "minimal"))

    names(df)[1] <- "gene"                         # unnamed first column holds gene symbols
    names(df)    <- make.names(names(df), unique = TRUE)
    df$gene      <- trimws(as.character(df$gene))
    df           <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
    if (anyDuplicated(df$gene)) warning("Duplicate gene symbols in ", sheets[i])

    out_file <- file.path(out_dir, paste0(file_names[i], ".csv"))
    write_csv_safe(df, out_file)
    message("  -> ", out_file, " (", nrow(df), " genes)")
  }, error = function(e) {
    message("Sheet ", sheets[i], " failed: ", conditionMessage(e))
    failed <<- c(failed, sheets[i])
  })
}

if (length(failed) > 0) warning("Failed sheets: ", paste(failed, collapse = ", "))

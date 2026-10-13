# Readers of decimal text that the tests hold pilotr's output to.

# Read a decimal the way R does when built without long-double arithmetic, as on CRAN's
# no-long-double check and on Apple silicon. R accumulates the digits in a double and then
# applies the decimal exponent in one division or multiplication. Once the digits exceed 2^53,
# as every 17-digit number's do, the accumulation itself rounds, and the result can be one unit
# in the last place from the nearest double. The same holds for an exponent beyond +/-22, whose
# power of ten is not exact.
no_long_double <- function(s) {
  neg <- substr(s, 1, 1) == "-"
  if (neg) s <- substring(s, 2)
  expn <- 0
  if (grepl("[eE]", s)) {
    halves <- strsplit(s, "[eE]")[[1]]
    expn <- as.numeric(halves[2])
    s <- halves[1]
  }
  parts <- strsplit(s, ".", fixed = TRUE)[[1]]
  if (length(parts) > 1) expn <- expn - nchar(parts[2])
  ans <- 0
  for (d in strsplit(paste0(parts, collapse = ""), "")[[1]]) ans <- 10 * ans + as.numeric(d)
  fac <- 10^abs(expn)
  ans <- if (expn < 0) ans / fac else ans * fac
  if (neg) -ans else ans
}

# The decimal numbers R's parser reads in a piece of R code, as their source text. Strings,
# raw strings among them, and comments hold no such tokens, so the JSON text a generated script
# embeds is not among them, nor are integer constants such as 1L.
r_decimal_tokens <- function(code) {
  pd <- utils::getParseData(parse(text = code, keep.source = TRUE))
  toks <- pd$text[pd$token == "NUM_CONST"]
  toks[grepl("^[0-9.][0-9.eE+-]*$", toks)]
}

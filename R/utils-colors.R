#' Get Default Muted Colors
#'
#' Returns a named vector of colors based on the number of groups.
#' Uses a muted, pleasing palette (Okabe-Ito inspired or similar).
#'
#' @param groups A character vector of group names.
#'
#' @return A named character vector associated with the groups.
#' @keywords internal
get_default_colors <- function(groups) {
  n <- length(groups)
  
  # A pleasing muted palette (Okabe-Ito inspired + extras)
  # 1. Orange, 2. Sky Blue, 3. Bluish Green, 4. Yellow, 5. Blue, 6. Vermilion, 7. Reddish Purple, 8. Grey
  base_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#882255", "#AA4499", "#332288", "#DDCC77", "#117733", "#44AA99"
  )
  
  if (n > length(base_palette)) {
    warning("More groups than default colors available. Recycling colors.")
    colors <- rep(base_palette, length.out = n)
  } else {
    colors <- base_palette[1:n]
  }
  
  names(colors) <- groups
  return(colors)
}

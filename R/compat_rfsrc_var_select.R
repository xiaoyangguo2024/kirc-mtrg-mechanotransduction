# =============================================================================
# compat_rfsrc_var_select.R
# Restore randomForestSRC::var.select removed in randomForestSRC >= 3.4
# =============================================================================
local({
  ns <- asNamespace("randomForestSRC")
  if (!exists("var.select.rfsrc", envir = ns, inherits = FALSE)) {
    e <- new.env(parent = ns)
    sys.source("R/vendor/var.select.rfsrc_3.3.1.R", envir = e)
    for (nm in ls(e, all.names = TRUE)) {
      assign(nm, get(nm, envir = e), envir = ns)
    }
    if (!exists("var.select", envir = ns, inherits = FALSE)) {
      assign("var.select", get("var.select.rfsrc", envir = ns), envir = ns)
    }
    try(registerS3method("var.select", "rfsrc", get("var.select.rfsrc", envir = ns), envir = ns), silent = TRUE)
  }
  for (fn in c("var.select", "var.select.rfsrc")) {
    if (exists(fn, envir = ns, inherits = FALSE)) {
      assign(fn, get(fn, envir = ns), envir = .GlobalEnv)
    }
  }
})

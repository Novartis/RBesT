# ensure that an exisitng user .Rprofile is respected
if(exists("~/.Rprofile")) source("~/.Rprofile")

# sets up brms caching and use of cmdstanr if available
local({
    if(requireNamespace("cmdstanr", quietly=TRUE)) {
        # instruct brms to use cmdstanr as backend and cache all Stan binaries
        brms_cache_dir <- Sys.getenv("BRMS_CACHE_DIR", here::here(".brms_cache"))
        # create cache directory if not yet available
        dir.create(brms_cache_dir, FALSE)
        options(brms.backend="cmdstanr", cmdstanr_write_stan_file_dir=brms_cache_dir)
    }
})

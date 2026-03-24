refresh_shiny_app <- function(
  project_dir = ".",
  app_file = "app.r",
  app_name,
  account,
  branch = "main",
  mode = c("push_and_deploy", "pull_and_deploy"),
  commit_message = paste("Refresh app", Sys.time())
) {
  mode <- match.arg(mode)
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  
  git <- function(args) {
    status <- system2("git", c("-C", project_dir, args))
    if (!identical(status, 0L)) {
      stop("Git command failed: git ", paste(args, collapse = " "), call. = FALSE)
    }
    invisible(status)
  }
  
  if (mode == "push_and_deploy") {
    git(c("add", "-A"))
    
    commit_status <- system2(
      "git",
      c("-C", project_dir, "commit", "-m", commit_message)
    )
    
    if (!identical(commit_status, 0L)) {
      status_lines <- system2(
        "git",
        c("-C", project_dir, "status", "--porcelain"),
        stdout = TRUE
      )
      if (length(status_lines) > 0) {
        stop("Commit failed.", call. = FALSE)
      }
    }
    
    git(c("push", "origin", branch))
  }
  
  if (mode == "pull_and_deploy") {
    git(c("pull", "--ff-only", "origin", branch))
  }
  
  deploy_result <- rsconnect::deployApp(
    appDir = project_dir,
    appPrimaryDoc = app_file,
    appName = app_name,
    account = account,
    forceUpdate = TRUE
  )
  
  list(
    mode = mode,
    project_dir = project_dir,
    branch = branch,
    app_name = app_name,
    account = account,
    deploy_result = deploy_result
  )
}

# ##examples , refresh both github and shinysapp.io

# result <- refresh_shiny_app(
#   project_dir = "C:/Users/BolducF/Documents/ShinyApps/ShinyTransects",
#   app_file = "app.r",
#   app_name = "YOUR_APP_NAME",
#   account = "YOUR_ACCOUNT_NAME",
#   branch = "main",
#   mode = "push_and_deploy",
#   commit_message = "Update app"
# )

# result

# ## refresh only shinyapps.io

# result <- rsconnect::deployApp(
#   appDir = "C:/Users/BolducF/Documents/ShinyApps/ShinyTransects",
#   appPrimaryDoc = "app.r",
#   appName = "YOUR_APP_NAME",
#   account = "YOUR_ACCOUNT_NAME",
#   forceUpdate = TRUE
# )

# result

# #or if rsconnect on github is up to date

result <- rsconnect::deployApp(
  appDir = "C:/Users/BolducF/Documents/ShinyApps/ShinyTransects",
  appPrimaryDoc = "app.r",
  appName = "shiny-transects",
  account = "lreyi8-fran0ois-bolduc",
  forceUpdate = TRUE
)

result
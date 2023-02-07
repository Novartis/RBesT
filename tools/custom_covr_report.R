library("digest")
library("covr")
library("DT")
library("htmltools")
library("readr")

RBesT_covr_report <- function(
  x = package_coverage(),
  file = "covr/covr-RBesT.html",
  pkg_tar = Sys.getenv("PKG_FILE_NAME"),
  pkg_md5 = Sys.getenv("PKG_FILE_MD5"),
  browse = FALSE
)
{
  # Create any directories as needed
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  # Paths need to be absolute for save_html to work properly
  file <- file.path(normalizePath(dirname(file), mustWork = TRUE), basename(file))

  loadNamespace("htmltools")
  loadNamespace("DT")

  data <- covr:::to_report_data(x)

  # Color the td cells by coverage amount, like codecov.io does
  color_coverage_callback <- DT::JS(
    'function(td, cellData, rowData, row, col) {
  var percent = cellData.replace("%", "");
  if (percent > 90) {
    var grad = "linear-gradient(90deg, #edfde7 " + cellData + ", white " + cellData + ")";
  } else if (percent > 75) {
    var grad = "linear-gradient(90deg, #f9ffe5 " + cellData + ", white " + cellData + ")";
  } else {
    var grad = "linear-gradient(90deg, #fcece9 " + cellData + ", white " + cellData + ")";
  }
  $(td).css("background", grad);
}
')

  # Open a new file in the source tab and switch to it
  file_choice_callback <- DT::JS(
    "table.on('click.dt', 'a', function() {
  files = $('div#files div');
  files.not('div.hidden').addClass('hidden');
  id = $(this).text();
  files.filter('div[id=\\'' + id + '\\']').removeClass('hidden');
  $('ul.nav a[data-value=Source]').text(id).tab('show');
});")

  package_name <- attr(x, "package")$package
  package_version <- attr(x, "package")$version
  percentage <- sprintf("%02.2f%%", data$overall)

  table <- DT::datatable(
    data$file_stats,
    escape = FALSE,
    fillContainer = TRUE,
    options = list(
      searching = FALSE,
      dom = "t",
      paging = FALSE,
      columnDefs = list(
        list(targets = 6, createdCell = color_coverage_callback))),
    rownames = FALSE,
    class = "row-border",
    callback = file_choice_callback
  )
  table$sizingPolicy$defaultWidth <- "100%"
  table$sizingPolicy$defaultHeight <- NULL

  md5_pkg <- read_file(pkg_md5)
  md5_pkg <- strsplit(md5_pkg, " ")[[1]][1]

  ui <- covr:::fluid_page(
    htmltools::includeCSS(system.file("www/report.css", package = "covr")),
    covr:::column(8, offset = 2, size = "md",
           htmltools::HTML(
             paste0(
               "<h2>", package_name, " ", package_version, " coverage - ", percentage, "</h2>",
               "<h3>", pkg_tar, " MD5: ", md5_pkg, "</h3>"
              )
            ),
           covr:::tabset_panel(
             covr:::tab_panel("Files",
                       table
             ),
             covr:::tab_panel("Source", covr:::addHighlight(covr:::renderSourceTable(data$full)))
           )
    )
  )

  htmltools::save_html(ui, file)

  if (browse) {
    viewer <- getOption("viewer", utils::browseURL)
    viewer(file)
  }

  invisible(file)
}


RBesT_covr_report()

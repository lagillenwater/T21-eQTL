# Sources every scripts/lib helper once per test run.
#
# Kept OUT of the individual test files for two reasons: (1) the relative path
# only has to be resolved in one place, whether tests run from the repo root
# (Rscript tests/testthat.R) or from tests/testthat; (2) covr::file_coverage
# instruments the lib files itself, and a test file that re-sources them would
# silently replace the instrumented definitions and zero the coverage report.
lib_dir <- if (dir.exists("scripts/lib")) "scripts/lib" else file.path("..", "..", "scripts", "lib")
for (f in list.files(lib_dir, pattern = "[.]R$", full.names = TRUE)) source(f)

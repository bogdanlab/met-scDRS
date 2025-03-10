score.loader <- function(score.dir) {
    score.files <- list.files(score.dir, pattern = '\\.score.gz', full.names = TRUE);
    risk.score <- vector('list', length = length(score.files));
    names(risk.score) <- score.files;

    # read into the empty list:
    for (result in score.files) {
        risk.score[[result]] <- read.table(
            file = result,
            sep = '\t',
            header = TRUE,
            row.names = 1,
            stringsAsFactors = FALSE
            );
        }

    # simplify list names:
    list.names <- gsub(score.dir, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;
    return(risk.score)
}

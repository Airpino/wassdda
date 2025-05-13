
setMethod("[", "MatD",
          function(x, i, j, ..., drop = TRUE) {
            mat <- x@M

            # Default: tutte le righe/colonne
            if (missing(i)) i <- seq_len(nrow(mat))
            if (missing(j)) j <- seq_len(ncol(mat))

            # Estrazione
            sub_mat <- mat[i, j, drop = FALSE]

            # Se singolo elemento e drop=TRUE, ritorna distributionD
            if (drop && length(i) == 1 && length(j) == 1) {
              return(sub_mat[[1]])
            }

            # Altrimenti crea nuovo MatD
            new_obj <- new("MatD")
            new_obj@M <- sub_mat

            return(new_obj)
          }
)

# helpers


suppressMessages(require(config))

#' Process Config files
#' 
#' Function to process config file as in the config R library
#' 
#' @param file Character; The config file path
#' @param config character; The config of the config file (see config::get())
#' 
#' @return list with two elements: (1) config.file: List of the config parameters 
#' and (2) str.out a string summary of the parameters in the config filel
ProcessConfigFile <- function(file, config){
  
  default.configs <- config::get(file = file, config = "default")
  config.file <- config::get(file = file, config = config)
  
  str.out <- c()
  
  for (i1 in names(default.configs)){
    if (class(config.file[[i1]]) == "list"){
      str.out <- c(str.out, paste0("- ", i1, "\n"))
      for (i2 in names(config.file[[i1]])){
        if (class(config.file[[i1]][[i2]]) == "list"){
          str.out <- c(str.out, paste0("- ", i2, "\n"))
          for (i3 in names(config.file[[i1]][[i2]])){
            str.out <- c(str.out, paste0("\t\t- ", i3, ": ", paste(config.file[[i1]][[i2]][[i3]], collapse = ", "), "\n")) 
            if (!i3 %in% names(config.file[[i1]][[i2]])){
              config.file[[i1]][[i2]][[i3]] <- default.configs[[i1]][[i2]][[i3]]
            }
          }
        } else {
          str.out <- c(str.out, paste0("\t- ", i2, ": ", paste(config.file[[i1]][[i2]], collapse = ", "), "\n"))
          if (!i2 %in% names(config.file[[i1]])){
            config.file[[i1]][[i2]] <- default.configs[[i1]][[i2]]
          }
        }
      }
    } else {
      str.out <- c(str.out, paste0("- ", i1, ": ", paste(config.file[[i1]], collapse = ", "), "\n"))
      if (!i1 %in% names(config.file)){
        config.file[[i1]] <- default.configs[[i1]]
      }
    }
  }
  
  return(list("config.file" = config.file, "str.out" = paste(str.out, collapse = "")))
}


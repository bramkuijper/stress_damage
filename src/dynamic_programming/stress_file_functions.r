library("readr")
#' Function that reads parameters from 
#' the simulation output file 'stressLXXXXXX.txt
#' @param filename a string containing the name of the 'stressLXXXX' file
#' @return data frame with parameter names as column names and parameter values
#' as the corresponding row

read.parameters <- function(filename)
{
    # read the file as a list of lines
    f = readLines(filename)

    # search for this pattern in the file
    parameter.pattern <- "^PARAMETER"

    # aux variable to store the parameter data
    parameter.data <- NULL 

    # go through all lines & hunt for parameters
    for (line_i in seq(1,length(f)))
    {
        # see whether current line matches the starting
        # pattern of the parameter listing
        if (length(grep(
                       pattern=parameter.pattern,
                       x=f[[line_i]])
               ) > 0)
        {
            # get the remainder of the file
            # as a data.frame
            parameter.data <- read.table(
                    file=filename
                    ,skip=line_i,
                    ,sep="\t")

            break;
        }
    }

    # tranpose the parameter data, doing this with t()
    # is a nightmare. Rather use the values to fill a matrix
    # and then give the columns names later
    parameter.data.t = as.data.frame(matrix(parameter.data$V2,nrow=1))

    # give the parameter data set corresponding column names
    # also remove the handy colons from the column names
    names(parameter.data.t) <- gsub(
            pattern="[: ]",
            ,replacement=""
            ,x=parameter.data$V1)

    return(parameter.data.t)
} # end read.parameters()


#' Gets summary mortality statistics
#' from the fwdCalc file
read.fwdcalc <- function(filename)
{
    f = readLines(filename)

    death.data <- read_table(
            file=filename
            ,col_names=F
            ,n_max=3
            ,skip=1)

    # tranpose the fwdcalc data, doing this with t()
    # is a nightmare. Rather use the values to fill a matrix
    # and then give the columns names later
    death.data.t = as.data.frame(matrix(death.data$X2,nrow=1))

    # give the parameter data set corresponding column names
    # also remove the handy colons from the column names
    names(death.data.t) <- gsub(
            pattern="\\/i:",
            ,replacement=""
            ,x=death.data$X1)

    return(death.data.t)
} # end read.fwdCalc

#' Puts the stressLXXXX data file
#' in a data frame
read.stress.file <- function(filename)
{
    # read the file as a list of lines
    f = readLines(filename)

    data.start.pattern = "^t\t"
    data.end.pattern = "^\\D"

    data.range <- c(NA,NA)

    for (line_i in seq(1,length(f)))
    {
        if (is.na(data.range[[1]]) &&
                        length(grep(pattern =data.start.pattern, x= f[[line_i]])) > 0)
        {
            data.range[[1]] <- line_i
            next
        }

        if (!is.na(data.range[[1]]) &&
                length(grep(pattern=data.end.pattern, 
                                x=f[[line_i]])) > 0)
        {
            data.range[[2]] <- line_i
            break
        }
    } # end for line_i

    print(data.range)

    data <- read_table(file=filename
            ,skip= data.range[[1]] - 1
            ,n_max = data.range[[2]] - data.range[[1]] - 2)

    return(data)
} # end read.stress.file




#' Puts the simAttacks data file
#' in a data frame. 
#' simAttacks data is a bit problematic
#' as it contains two consecutive datasets:
#' one for acute stress, one for chronic
read.attack.file <- function(filename)
{
    # read in file line by line to obtain ranges
    f = readLines(filename)

    # whether we are at the acute part
    at.acute <- T

    # variables to store the line ranges
    range.acute <- c()
    range.chronic <- c()

    # loop through all lines in file 
    for (line_i in seq(1,length(f)))
    {
        if (length(grep("^time",f[[line_i]])) > 0)
        {
            if (at.acute)
            {
                range.acute <- c(line_i)
                at.acute <- F
            }
            else
            {
                range.acute <- c(range.acute, line_i - 3)
                range.chronic <- c(line_i)
                break
            }
        }
    } # end for line_i

    n_max = Inf
        
    if (length(range.acute) == 2)
    {
        n_max = range.acute[[2]] - range.acute[[1]]
    }        
    
    # read in the acute stress data
    data.acute <- read_table(
            file = filename
            ,skip = range.acute[[1]] - 1
            ,n_max = n_max
            ,col_names=T)

    data.acute["type"] <- "acute"
   
    # if there is a defined length of range acute
    # this means that there is a date set of the chronic
    # exposure below it
    if (length(range.acute) == 2)
    {
        # read in the chronic stress data
        data.chronic <- read_table(
                file = filename
                ,skip = range.chronic[[1]] - 1
                ,col_names=T)
    
        data.chronic["type"] <- "chronic"
        data.acute <- rbind(data.acute,data.chronic)
    }
    return(data.acute)
} # end read.attack.file()


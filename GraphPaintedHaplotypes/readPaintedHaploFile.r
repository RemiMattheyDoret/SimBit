readFirstPartOfHaplotype = function(str)
{
    ## Examples: "0-500 P434 I0 H0", "anc-end P434 I0 H0"
    generations = strsplit(str[1], '-')[[1]]
    paintedGeneration = generations[1]
    observedGeneration = generations[2]

    stopifnot(substring(str[2], 1, 1) == 'P')
    stopifnot(substring(str[3], 1, 1) == "I")
    stopifnot(substring(str[4], 1, 1) == "H")
    patch = as.integer(substring(str[2], 2))
    ind = as.integer(substring(str[3], 2))
    haplo = as.integer(substring(str[4], 2))

    return(c(paintedGeneration, observedGeneration, patch, paste(ind, haplo, sep="_"
        )))
}

readSegment = function(str)
{
    ss = strsplit(str, " ")[[1]]
    stopifnot(length(ss) == 4)
    stopifnot(substring(ss[3], 1, 1) == "P")
    stopifnot(substring(ss[4], 1, 1) == "I")
    left = as.integer(ss[1])
    right = as.integer(ss[2])
    patch = as.integer(substring(ss[3], 2))
    id = as.integer(substring(ss[4], 2))
 
    stopifnot(left < right)   

    return(c(left, right, patch, id))
}


readHaplotype = function(str)
{
    twoStrings = strsplit(str, ": ")[[1]]
    firstPart = strsplit(twoStrings[1], " ")[[1]]
    segmentsPart = strsplit(substring(twoStrings[2], 2, nchar(twoStrings[2])-1), "][", fixed=TRUE)[[1]]

    r = as.data.frame(matrix(NA, nrow=length(segmentsPart), ncol = 8))

    first = readFirstPartOfHaplotype(firstPart)

    for (seg_index in 1:length(segmentsPart))
    {
        r[seg_index,] = c(first, readSegment(segmentsPart[seg_index]))
    }

    return(r)
}

readPaintedHaploFile = function(path)
{
    ## read file
    lines = readLines(path)

    ## Initialize list that will be merged
    ldata = vector(mode = "list", length = length(lines))

    ## Read haplotypes and fill upp list
    for (line_index in 1:length(lines))
    {
        ldata[[length(ldata) + 1]] = readHaplotype(lines[line_index])
    }

    ## Merge list
    data = do.call("rbind", ldata)
    colnames(data) = c("paintedGeneration", "observedGeneration", "sampled_patch", "sampled_haploID", "left", "right", "patch", "ID")
    data$paintedGeneration = as.integer(data$paintedGeneration)
    data$observedGeneration = as.integer(data$observedGeneration)
    data$sampled_patch = as.integer(data$sampled_patch)
    data$left = as.integer(data$left)
    data$right = as.integer(data$right)
    data$patch = as.integer(data$patch)
    data$ID = as.integer(data$ID)
    return(data)
}
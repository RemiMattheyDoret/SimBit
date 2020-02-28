
GetParameterGridData = function(SGrid,Parameter_line)
{
	for (v in names(SGrid))
	{
		assign(v,SGrid[v][Parameter_line,], envir = .GlobalEnv)
	}
}

ParameterGrid = function(..., bigID = NULL)
{
	input = list(...)
	Names = rep("",length(input))
	for (big in 1:length(input))
	{
		name = names(input[[big]][1])
		if (length(name)<1) stop("error: Name missing")
		if (length(name)>1) stop("error: length(name)>1")
		Names[big] = name
	}
	FirstCol = lapply(input, "[[", 1)
	SGrid = expand.grid(FirstCol)
	Cols = FirstCol
	for (big in 1:length(input))
	{
		if (length(input[[big]]) > 1)
		{
			for (small in 2:length(input[[big]]))
			{
				Cols[[big]] = input[[big]][[small]]
				name = names(input[[big]][small])
				if (length(name)<1) stop("error: Name missing")
				if (length(name)>1) stop("error: length(name)>1")
				Names[length(Names)+1] = name
				SGrid = cbind(SGrid,expand.grid(Cols)[,big])
			}
		}
	}
	names(SGrid) = Names

	if (!is.null(bigID))
	{
		SGrid$bigID = bigID
		SGrid$smallID = 1:nrow(SGrid)
		SGrid$ID = paste0(SGrid$bigID, ".", SGrid$smallID)
	}

	return(SGrid)
}





###########################
##### Instal packages #####
###########################

needed_packages <- c("R6", "stringr")
need_to_install_packages = setdiff(needed_packages, rownames(installed.packages()))
if (length(need_to_install_packages) > 0) {
	cat(paste("I will try to install the following package(s):", paste(need_to_install_packages, collapse=", ")))
  install.packages(need_to_install_packages)  
}
for (needed_package in needed_packages)
{
	require(needed_package, character.only=TRUE)
}



#########################
##### ParameterFile #####
#########################

ParameterFile = R6Class(
	### Class name
	"ParameterFile",

	### private
	private = list(
		filePath = "",
		fileContent = list()
	),

	public = list(
		initialize = function(path)
		{
			filePath = path
			private$fileContent = readLines(path)
		},

		readParameter = function(parameterName)
		{
			stopifnot(class(parameterName) == "character")
			stopifnot(length(parameterName) == 1)

			index = grep(paste0("^",parameterName), private$fileContent)
			if (length(index) > 1)
			{
				stop(paste0("Found several line matching parameter name ", parameterName, "  in file '", private$filePath,"'. Matches found are at lines '",paste(index, collapse=" "),"'. That could be an internal error (either SimBit outputing several times the same variable name or this little readParameter function might not be working properly"))
			}
			if (length(index) == 0)
			{
				stop(paste0("Could not find parameter name ", parameterName, " in file '", private$filePath,"'"))
			}
			
			return(eval(parse(text = private$fileContent[[index+1]])))
		},

		readDispersalMatrixOverTime = function()
		{
			superlist = self$readParameter("fullFormForwardMigrationMatrix")

			r = list()
			for (i in 1:length(superlist))
			{
				d = data.frame(superlist[[i]])
				#rownames(d, prefix = "from_Patch_")
				#colnames(d, prefix = "to_Patch_")
				rownames(d) = paste0("from_Patch_", 1:nrow(d))
				colnames(d) = paste0("to_Patch_", 1:ncol(d))
				r[[i]] = d
			}
			return (r)
		}
	)
)

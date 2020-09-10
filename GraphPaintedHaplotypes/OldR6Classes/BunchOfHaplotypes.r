
BunchOfHaplotypes = R6::R6Class(
    ### Class name
    "BunchOfHaplotypes",

    ### Public members
    public = list(
        ### Attributes
        haplotypes = list(),
        
        ### Print
        print = function()
        {
            for (haplotype in self$haplotypes)
            {
                haplotype$print()
            }
        },

        ### Read file
        initialize = function(path)
		{
			## read file
			lines = readLines(path)

			## Initialize list to return
			self$haplotypes = rep(NA, length(lines))

			for (line_index in 1:length(lines))
			{
				self$haplotypes[line_index] = Haplotype$new(lines[line_index])
			}
		},

		### Copy haplotypes
        initialize = function(haplos)
		{
			self$haplotypes = haplos
		},


		### Remove part of the genome
		subsetGenome = function(from, to)
		{
			newHaplotypes = rep(NA, length(self$haplotypes))
			for (haplotype_index in 1:length(self$haplotypes))
			{
				newHaplotypes[haplotype_index] = self$haplotypes[haplotype_index]$subsetGenome(from, to)
			}

			return (BunchOfHaplotypes$new(newHaplotypes))
		},



		subsetHaplos = function(expr)
		{
			r = rep(NA, length(self$haplotypes))
			for (haplotype_index in 1:length(self$haplotypes))
			{
				if (self$haplotypes[haplotype_index].paintedGeneration )
			}
			return(r)
		}

		# getPaintedGenerations = function()
		# {
		# 	r = rep(NA, length(haplotypes))
		# 	for (haplotype_index in 1:length(haplotypes))
		# 	{
		# 		r[haplotype_index] = haplotypes[haplotype_index].paintedGeneration
		# 	}
		# 	return(r)
		# }
    )
)




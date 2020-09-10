Haplotype = R6::R6Class(
    ### Class name
    "Haplotype",

    ### Private members
    ### Add Segment
    private = list(
        readSegment = function(str)
        {
            ss = strsplit(str), " ")[[1]]
            stopifnot(lenght(ss) == 4)
            stopifnot(substring(ss[3], 1, 1) == "P")
            stopifnot(substring(ss[4], 1, 1) == "I")
            left = as.integer(ss[1])
            right = as.integer(ss[2])
            patch = as.integer(substring(ss[3], 2))
            id = as.integer(substring(ss[4], 2))

            return(c(left, right, patch, id))
        },

        readFirstPartOfHaplotype = function(str)
        {
            ## Examples: "0-500 P434 I0 H0", "anc-end P434 I0 H0"
            ss = strsplit(str, " ")

            generations = strsplit(ss[1], "-")[[1]]
            paintedGeneration = generations[1]
            observedGeneration = generations[2]

            stopifnot(substring(ss[2], 1, 1) == "P")
            stopifnot(substring(ss[3], 1, 1) == "I")
            stopifnot(substring(ss[4], 1, 1) == "H")
            patch = as.integer(substring(ss[2], 2))
            ind = as.integer(substring(ss[3], 2))
            haplo = as.integer(substring(ss[4], 2))
        }
    ),

    ### Public members
    public = list(
        ### Attributes
        segments = c(),
        paintedGeneration = -1,
        observedGeneration = -1,
        patch = -1,
        ind = -1,
        haplo = -1,

        ### Print
        print = function()
        {
            cat(paste0(self$paintedGeneration, "-", self$observedGeneration, " P", self$patch, " I", self$ind, " H", self$haplo, ": "))
            for (seg in self$segments)
            {
                seg$print()
            }
            cat("\n")
        },

        ### Read Haplotype
        initialize = function(str)
        {
            twoStrings = strsplit(str, ": ")[[1]]
            firstPart = strsplit(twoStrings[1], " ")[[1]]
            segmentsPart = strsplit(substr(twoStrings[2], 2, length(twoStrings[2])-1), " ")[[1]]

            readFirstPartOfHaplotype(firstPart)

            for (seg_str in segmentsPart)
            {
                readAndPushBackSegment(seg_str)
            }
        },

        ### Copy attributes
        initialize = function(segs, pG, oG, p, i, h)
        {
            self$segments = segs
            self$paintedGeneration = pG
            self$observedGeneration = oG
            self$patch = p
            self$ind = i
            self$haplo = h
        },

        ### Remove part of the genome
        subsetGenome = function(from, to)
        {
            newSegments = c()
            for (seg in self$segments)
            {
                if (seg.right < from) next
                if (seg.left >= to) break

                if (seg.left < from) seg.left = from
                if (seg.right > to) seg.right = to

                newSegments[length(newSegments) + 1] = seq
            }

            return (Haplotype$new(newSegments, paintedGeneration, observedGeneration, patch, ind, haplo))
        }
        
    )
)
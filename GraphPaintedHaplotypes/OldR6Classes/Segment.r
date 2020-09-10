Segment = R6::R6Class(
    ### Class name
    "Segment",

    ### Public members
    public = list(
        ### Attributes
        left = -1,
        right = -1,
        patch = -1,
        ID = -1,

        ### Initialize
        initialize = function(l, r, p, id)
        {
            self$left = l
            self$right = r
            self$patch = p
            self$ID = id

            stopifnot(self$left < self$right)
        },

        ### Print
        print = function()
        {
            cat(paste0("[",self$left, " -> ", self$right, " P", self$patch, " ", self$ID, "]")) 
        },


        ### Get width of segment
        width = function()
        {
            return(self$right - self$left)
        }
    )
)




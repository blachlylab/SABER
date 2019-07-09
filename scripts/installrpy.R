if (!require("rPython"))
{
    install.packages("rPython",repos='http://cran.us.r-project.org')
    if(!require("rPython")) stop("Package not found")
}


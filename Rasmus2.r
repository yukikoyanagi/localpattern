

source("StartClustSubsetVers4.txt")


# Rscript Rasmus2.r foobar/pattern
filnavn=commandArgs(TRUE)[1]
dirbase=dirname(filnavn)
filnavn=basename(filnavn)

inputdir=dirbase
outputdir=dirbase

inputfil=paste(inputdir, filnavn, sep="/")
outputfil=paste(outputdir, filnavn, sep="/")

rot=matrix(scan(inputfil),,4,byrow=T)
starttime=Sys.time()
StartNew(outputfil, rot)
stoptime=Sys.time()
difftime(stoptime, starttime, units="secs")


# Produces in labud:
# nummer, box-x, box-y, box-z, antal i mode box, samlede antal, 
# antal bokse, tæthed i modebox, sum af tætheder, volumen, peakness


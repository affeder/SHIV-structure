#Set up directory structure

#This script assumes it's in project directory = pd/code
#It assumes that there is already a pd/dat/ folder with appropriate files

#It creates the following:
createIfDoesntExist <- function(x){
    if(!dir.exists(x)){
        dir.create(x)
    }
}

createIfDoesntExist("../out")
#Main results (graphs and tables)
createIfDoesntExist("../out/graphs/")
createIfDoesntExist("../out/tables/")
#Compartmentalization test results
createIfDoesntExist("../out/amova/")
createIfDoesntExist("../out/kst/")
createIfDoesntExist("../out/sm/")
#BEAST results
createIfDoesntExist("../out/beast_out/")

createIfDoesntExist("../tmp")

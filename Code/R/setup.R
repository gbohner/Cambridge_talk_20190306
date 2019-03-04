setwd("~")
system("cd ~")
if (!file.exists("Cambridge_talk_20190306")){
  system("git clone https://github.com/gbohner/Cambridge_talk_20190306")
}
setwd("~/Cambridge_talk_20190306")
system("git pull")
setwd("~")
system("cp ~/Cambridge_talk_20190306/Code/R/* ~/")

# Set the libPaths to the already installed packages
curLibPaths = .libPaths()
curLibPaths[[2]] = paste0(curLibPaths[[2]],"/3.5")
.libPaths(curLibPaths)

# Therefore we do NOT need to reinstall packages on every user (takes very long)
#source("install_packages.R")

# Open the file that we want to edit
file.edit("scrnaseq_dimred.R")




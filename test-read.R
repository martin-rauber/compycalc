cooldown = read.csv(file.path(parent.folder, "cooldown_data.csv"), sep = ",", header = T)
cooldown= cooldown[,3:5]
cooldown
write.csv(cooldown, file =  paste("cooldown-test-file.csv",sep=""), row.names = T)
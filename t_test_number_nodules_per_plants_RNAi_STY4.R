set.seed(1407)

outfile <- "logs/t_test_number_nodules_per_plants_RNAi_STY4.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

control <- c(10,17,3,11,22,23,6,23,3,19,23,13,11,15,24,28,4,3,7,13,8,6,25,14,7,23,17,19,11,6,14,12,11,16,5,4,18,32)
STY4 <- c(3,7,10,1,5,10,3,2,16,18,11,8,6,14,16,17,22,5,0,12,16,2,2,15,0,5,10,11,15,11,9,3,5,16,4,12,3,25,0,18,0,9,6,23,8,11,15)

boxplot(control, STY4, names=c("Control", "STY4"))

t.test(control, STY4, var.equal=TRUE)
t.test(control, STY4, var.equal=FALSE)

closeAllConnections()

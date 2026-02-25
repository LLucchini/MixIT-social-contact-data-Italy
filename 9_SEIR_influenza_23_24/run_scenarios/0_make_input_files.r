rm(list=ls())

scen.school <- c("no_closure","uni","uni+upper_secondary",
                 "uni+secondary","uni+secondary+primary","all_closed")

scen.smart <- c("pre-lockdown","lifting lockdown", "lockdown")


# Creiamo tutti i possibili scenari combinati
combined_scenarios <- expand.grid(smart = scen.smart,school = scen.school)

combined_scenarios <- combined_scenarios[,2:1]
# Visualizziamo
print(combined_scenarios)

# Numero totale di scenari
num_scenarios <- nrow(combined_scenarios)
num_scenarios

nages <- 15

pop.scen.combined <- array(NA,c(length(scen.school),length(scen.smart),nages,2))
for(s in 1:length(scen.school)){
  for(t in 1:length(scen.smart)){
    pop.scen.combined[s,t,,] <- as.matrix(read.csv(paste("../../1_transmissibility_reductions/0_population_scenarios/pop_scenarios/scenario_smart_",t,"_sclose_",scen.school[s],".csv",sep="")))
  }
}

for(s in 1:num_scenarios){
  system(paste("mkdir scenario_",s,sep=""))
  system(paste("mkdir scenario_",s,"/input_files",sep=""))
  id.x <- which(scen.school==combined_scenarios[s,1])
  id.y <- which(scen.smart==combined_scenarios[s,2])
  pop.presence <- pop.scen.combined[id.x,id.y,,1]
  pop.not.presence <- pop.scen.combined[id.x,id.y,,2]
  write.table(pop.presence,file=paste("scenario_",s,"/input_files/pop_by_age_presence",sep=""),sep="\t",row.names = F,col.names = F)
  write.table(pop.not.presence,file=paste("scenario_",s,"/input_files/pop_by_age_not_presence",sep=""),sep="\t",row.names = F,col.names = F)
  system(paste("cp common_input/* scenario_",s,"/input_files/",sep=""))
  system(paste("cp ../sw/exec scenario_",s,sep=""))
  system(paste("cp common_input/run.sh scenario_",s,sep=""))
  
  
}

write.csv(combined_scenarios,file="readme_scenarios.csv")


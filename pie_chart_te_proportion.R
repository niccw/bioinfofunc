library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(see)
library(ggsci)

# repeatcraft output? raw repeatmasker output?
repeatcraft <- FALSE

# genome size
species_name <- "neoforv3"
genomeSize <- 37063664439  

if (repeatcraft){
  sp <- ""
}else{
  sp <- "neoforv3hm.repeatmasker.out.gff3.stat"
}


# colour map --------------------------------------------------------------
superfamily_ls <- c("Unknown","DNA","LTR","LINE","SINE","SINE?","RC","ARTEFACT","rRNA","snRNA","Simple_repeat","Satellite","Low_complexity","non_TE")
superfamily_color <- c("#B8B8B8B2","#D43F3AB2","#EEA236B2","#5CB85CB2","#46B8DAB2","#357EBDB2","#67464F","#428BC1","#9632B8B2","#9632B8B2","#D2E7D6","#D2E7D6","#D2E7D6","#FAFAEB")

names(superfamily_color) <- superfamily_ls


# read files --------------------------------------------------------------

# read as datafranme, add nonTE and species column, combine two df
sp_stat <- readr::read_tsv(sp)
sp_non_te <- genomeSize - sum(filter(sp_stat,level=="class")[["total_bp"]])



# plot non_te + all class -------------------------------------------------

#### plot non_TE + all class ####
#add percentage column to each species (for pie chart)
sp_stat <- rbind(sp_stat,c("class","non_TE",NA,sp_non_te,NA))

sp_stat$total_bp <- as.numeric(sp_stat$total_bp)
sp_stat <- mutate(sp_stat,perc=total_bp/genomeSize)

#### export cnt table ####
sp_stat <- mutate(sp_stat,perc100=perc*100)

write.csv(sp_stat,paste0(species_name,"_te_tbl.csv"),quote = F, row.names = F)

sp_stat$species <- species_name


# all fuper family is gone after factoring, be careful.
all_stat <- rbind(sp_stat)

# select superfamily manually
level_ls <- all_stat[all_stat["level"]=="class",c("name")] %>% pull()
print(level_ls)

cat("Remove unwanted superfamily by e.g level_ls[level_ls != `name`]\n")

if(match(level_ls, names(superfamily_color)) %>% is.na() %>% any()){
  cp <- superfamily_color[match(level_ls, names(superfamily_color))[!is.na(match(level_ls, names(superfamily_color)))]]
  # report what is missing
  cat("The colours of following superfamilies are not yet assigned:\n")
  cat(level_ls[match(level_ls, names(superfamily_color)) %>% is.na() %>% which()])
}else{
  cp <- superfamily_color[match(level_ls, names(superfamily_color))[!is.na(match(level_ls, names(superfamily_color)))]]
  cat("All good")
  }

if(length(cp) != length(level_ls)){
  cat("Not enough color! Check variable cp and level_ls \n")
  cat("level_ls")
  print(level_ls)
  cat("cp")
  print(cp)
}

all_stat$name <- factor(all_stat$name,levels = level_ls)


g <- ggplot(filter(all_stat,level=="class"),aes(x="",y=perc,fill=name))
g + geom_bar(stat = "identity",colour="#434744") + coord_polar("y",start=0) +
  scale_fill_manual(values=cp) +
  facet_grid(cols=vars(species)) +
  theme_modern() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave("te_composition_pie_w_nonTE_w0_annotation.svg",width=172,height = 192, units = "mm")

##### process df for all TE class ####
# start from begining, dont compute nonTE
sp_all_te <- sum(filter(sp_stat,level=="class")[["total_bp"]])
sp_stat <- mutate(sp_stat,perc=total_bp/sp_all_te)


sp_stat$species <- "sp"
all_stat <- rbind(sp_stat)
print(all_stat$name)

# select superfamily manually
level_ls <- all_stat[all_stat["level"]=="class",c("name")] %>% pull()
print(level_ls)

cat("Remove unwanted superfamily by e.g level_ls[level_ls != `name`]\n")

if(match(level_ls, names(superfamily_color)) %>% is.na() %>% any()){
  cp <- superfamily_color[match(level_ls, names(superfamily_color))[!is.na(match(level_ls, names(superfamily_color)))]]
  # report what is missing
  cat("The colours of following superfamilies are not yet assigned:\n")
  cat(level_ls[match(level_ls, names(superfamily_color)) %>% is.na() %>% which()])
}else{
  cp <- superfamily_color[match(level_ls, names(superfamily_color))[!is.na(match(level_ls, names(superfamily_color)))]]
  cat("All good")
}

if(length(cp) != length(level_ls)){
  cat("Not enough color! Check variable cp and level_ls \n")
  cat("level_ls")
  print(level_ls)
  cat("cp")
  print(cp)
}

g <- ggplot(filter(all_stat,level=="class"),aes(x="",y=perc,fill=name))
g + geom_bar(stat = "identity",colour="#434744") + coord_polar("y",start=0) +
  scale_fill_manual(values=cp) +
  facet_grid(cols=vars(species)) +
  theme_modern() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())



# ggsave
ggsave("te_composition_pie_w0_annotation.svg",width=268,height = 99, units = "mm")


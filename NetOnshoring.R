# This is a replication script for "Which Countries Have Offshored Carbon Dioxide Emissions in Net Terms?"
# The script imports the OECD ICIO data and creates the paper's tables and figures.
# Author: Aldy Darwili

# Load Libraries ----
packages = c("tidyverse",
             "openxlsx",
             "ggthemes",
             "magrittr",
             "diagonals",
             "rsdmx",
             "ggrepel")


# Now load or install&load all
package.check = lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Set timeout limit ---------
options(timeout = max(3000, getOption("timeout")))

# Load Raw Data -----
# We use ICIO 2021 data from OECD. We use the ICIO 2021 in R format.

## Intermediate input Z -------
load("ICIO2021econZ.RData") # this can be downloaded from the OECD ICIO database

## Final demand ---------
load("ICIO2021econFD.RData") # this can be downloaded from the OECD ICIO database

# Basic Information  -------
## Period of analysis ---------
year <- data.frame("year" = seq(from = 1995, to = 2018),
                   "n" = seq(from = 1, to = 24))

n_year <- seq(from = 1, to = 24)

## m x n row names (excluding CN1, CN2, MX1, MX2) -----
rnames <- data.frame("rowname" = rownames(ICIO2021econZ[1,,])) %>% 
  filter(!str_detect(rowname, "CN1|CN2|MX1|MX2")) %>% # excluding CN1, CN2, MX1, MX2
  separate(rowname, into=c("country", "sector"), sep="_")

# m x n Row names (Change CN1 and CN2 to CHN, and MX1 and MX2 to MEX)
rnames_alt <- bind_rows(rnames, rnames %>% filter(country == "MEX"), rnames %>% filter(country == "MEX"),
                        rnames %>% filter(country == "CHN"), rnames %>% filter(country == "CHN")) %>% 
  unite("rowname", country, sector, sep = "_") %>% 
  as_vector()

## Country -------
country <- rnames %>% 
  select(country) %>%
  unique()

country_long <- rnames %>% 
  select(country) %>% 
  as_vector()

## Sector ------------
sector <- read.xlsx("ICIO_sector_full.xlsx")

sector_long <- rnames %>% 
  select(sector) %>% 
  as_vector()

# Preparing Main Dataset (Z, FD, X) --------
list_Z <- list()
list_fd <- list()
list_x <- list()

for (i in n_year) {
  
  # Intermediate Input Matrix, Z
  Z <- ICIO2021econZ[i,,] # this line is to subset matrix from Z array
  rownames(Z) <- rnames_alt # naming the rows that contain CN1 and CN2 to CHN, and MX1 and MX2 to MEX
  colnames(Z) <- rnames_alt # naming the column that contain CN1 and CN2 to CHN, and MX1 and MX2 to MEX
  Z <- rowsum(Z, row.names(Z), reorder = FALSE) %>% # joining China and Mexico data: row
    t() %>% 
    rowsum(., row.names(Z), reorder = FALSE) %>% # joining China and Mexico data: column
    t()
  
  list_Z[[as.character(i)]] <- Z # storing the processed Z data into a list
  
  # Final Demand Matrix, Y
  FD <- ICIO2021econFD[i,,] # subsetting final demand matrix from final demand array
  rownames(FD) <- rnames_alt # naming the rows that contain CN1 and CN2 to CHN, and MX1 and MX2 to MEX
  FD <- rowsum(FD, row.names(FD), reorder = FALSE) # joining China and Mexico data: row
  
  list_fd[[as.character(i)]] <- FD[,-c(68:71)] # storing the processed FD data into a list; -c(68:71) is to remove CN1, CN2, MX1, MX2 (they all contain 0 values)
  
  # gross output X
  list_x[[as.character("")]] <- rowSums(Z) + rowSums(FD) # X (gross output) vector to the list object
  
  rm(Z, FD)
  gc()
  
}

rm(ICIO2021econFD, ICIO2021econZ)

# MRIO model ------
## Emission factor (emission intensity) data from OECD ----------
ef_oecd <- readRDS("EF_CO2_OECD.rds") # emission factor data can be downloaded from the OECD ICIO

## Unadjusted MRIO model ------
list_e <- list()
list_eeio <- list()

for (i in n_year) {
  
  # emission intensity
  q <- ef_oecd %>%
    filter(n == i) %>% 
    select(q) %>% 
    as_vector()
  
  # input coefficient matrix A
  A <- as.matrix(t((1/list_x[[i]]) * t(list_Z[[i]])))
  
  #removing NA value
  A[is.na(A)] <- 0
  
  # Leontief inverse matrix (merged China and Mexico)
  L <- solve(diag(NROW(A)) - A)
  
  # emission account
  eeio <- as.matrix(q * L %*% (list_fd[[i]])) %>% # q * L * y
    set_rownames(country_long) %>% # setting row names
    rowsum(., row.names(.), reorder = FALSE) # row sum by sector, results in m x m eeio matrix
  
  pbe <- eeio %>% # production-based emissions
    rowSums() %>% # sum over column to get PBE
    data.frame() %>% 
    rename(pbe_mrio = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(pbe_mrio = pbe_mrio / 1.0e+6)
  
  cbe <- eeio %>% # consumption-based emissions
    colSums() %>% 
    data.frame() %>% 
    rename(cbe_mrio = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(cbe_mrio = cbe_mrio / 1.0e+6)
  
  list_eeio[[as.character(i)]] <- pbe %>% 
    left_join(cbe, by = c("country", "n_year"))
  
  fatdiag(eeio, steps = 67) <- 0 # removing on-diagonal values (on-diagonal: domestic-domestic emissions)
  
  eex <- eeio %>% # emissions embodied in exports
    rowSums() %>% 
    data.frame() %>% 
    rename(eex_mrio = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(eex_mrio = eex_mrio / 1.0e+6)
  
  eem <- eeio %>% # emissions embodied in imports
    colSums() %>% 
    data.frame() %>% 
    rename(eem_mrio = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(eem_mrio = eem_mrio / 1.0e+6)
  
  list_e[[as.character("")]] <- eex %>% 
    left_join(eem, by = c("country", "n_year"))
  
  rm(q, A, L, eeio, pbe, cbe, eex, eem)
  gc()
  
}

saveRDS(list_e, "list_e.rds")
saveRDS(list_eeio, "list_eeio.rds")

## Technology-adjusted MRIO model -----

list_e_ctr_std <- list()
list_e_std <- list()

for (i in as_vector(country)) {
  for(j in c(seq(from = 1, to = 24))) {
    
    # adjusting emission intensity vector
    q_adj <- ef_oecd %>%
      filter(n == j) %>%
      filter(country == i) %>%
      select(q) %>%
      slice(rep(row_number(), 67)) %>% # replicating row m times
      as_vector()
    
    # input coefficient matrix A
    A <- as.matrix(t((1/list_x[[j]]) * t(list_Z[[j]])))
    
    A[is.na(A)] <- 0
    rownames(A) <- sector_long
    colnames(A) <- country_long
    
    # n x nm matrix of technological coefficient (as per Dietzenbacher)
    H <- A %>%
      rowsum(., row.names(.), reorder = FALSE) %>% # row sum by country; resulting in n x nm matrix
      rbind(., .[rep(1:45, 66),]) # replicating vertically from n x nm to nm x nm matrix
    
    # adjusted n x nm matrix of technological coefficient (e.g., all countries have a similar n x n matrix of technological coefficient)
    H_adj <- A[, colnames(A) %in% c(i)] %>% # subsetting the column for the country of interest; resulting in nm x n matrix
      rowsum(., row.names(.), reorder = FALSE) %>% # row sum by country; resulting in n x n matrix
      .[,rep(1:45, 67)] %>% # replicating horizontally; from n x n to n x nm matrix
      rbind(., .[rep(1:45, 66),]) # replicating vertically from n x nm to nm x nm matrix
    
    TS <- A / H # matric of intermediate trade structure
    
    is.na(TS) <- sapply(TS, is.infinite) #removing possible infinite value
    TS[is.na(TS)] = 0 # removing possible NA values
    
    # Adjusted technical coefficient matrix
    A_adj <- TS * H_adj
    
    # Adjusted Leontief inverse matrix
    L_adj <- solve(diag(NROW(A_adj)) - A_adj)
    
    # Adjusted Emission Account
    e_adj <- (q_adj * L_adj %*% (list_fd[[j]]))
    
    fatdiag(e_adj, steps = 67) <- 0 # this is to remove on-diagonal elements of the nm x m matrix of emission embodied input-output
    
    list_e_ctr_std[[as.character(j)]] <- e_adj %>%
      data.frame() %>%
      select(all_of(i))
    
    list_e_std[[as.character(i)]] <- list_e_ctr_std
    
    rm(q_adj, A, A_adj, H, H_adj, TS, L_adj, e_adj)
    gc()
    
  }
}

saveRDS(list_e_std, "list_e_std.rds")

rm(list_e_ctr_std)

list_eem_adj <- list()

for (i in as_vector(country)) {
  
  list_eem_adj[[i]] <- bind_cols(list_e_std[[i]]) %>% 
    colSums() %>% 
    data.frame() %>% 
    rename(eam_mrio = 1) %>% 
    mutate(country = c(i),
           n_year = c(seq(from = 1, to = 24))) %>% 
    select(country, n_year, eam_mrio)
  
}

EEM_adj <- bind_rows(list_eem_adj) %>% 
  as_tibble()

## Household Direct Emissions (from OECD) ------
DMHH_cb <- readRDS(paste0("DMHH_fd.rds")) %>% 
  rename("DMHH_cb" = "DMHH") %>% 
  mutate(year = as.numeric(year))

DMHH_pb <- readRDS(paste0("DMHH_prod.rds")) %>% 
  rename("DMHH_pb" = "DMHH") %>% 
  mutate(year = as.numeric(year))

mrio <- bind_rows(list_e) %>% # gathering all results from the MRIO model
  left_join(bind_rows(list_eeio), by = c("country", "n_year")) %>% 
  left_join(EEM_adj, by = c("country", "n_year")) %>% 
  mutate(eam_mrio = eam_mrio/1.0e+6) %>% 
  left_join(year %>% 
              rename(n_year = 2), by = "n_year") %>% 
  select(country, year, n_year, pbe_mrio, cbe_mrio, eex_mrio, eem_mrio, eam_mrio) %>% 
  left_join(bind_rows(DMHH_pb), by = c("country", "year")) %>% 
  left_join(bind_rows(DMHH_cb), by = c("country", "year")) %>% 
  mutate(pbe_mrio = pbe_mrio + DMHH_pb,
         cbe_mrio = cbe_mrio + DMHH_cb,
         beet_mrio = eex_mrio - eem_mrio,
         NEON_mrio = eex_mrio - eam_mrio,
         apbe_mrio = pbe_mrio - NEON_mrio) %>% 
  select(country, year, n_year, pbe_mrio, cbe_mrio, eex_mrio, eem_mrio, eam_mrio, beet_mrio, NEON_mrio, apbe_mrio)

saveRDS(mrio, "mrio.rds")

# EEBT Model ----
## Unadjusted EEBT Model ----
list_e_eebt <- list()
list_eebt <- list()

for (i in n_year) {
  
  # emission intensity
  q <- ef_oecd %>%
    filter(n == i) %>% 
    select(q) %>% 
    as_vector()
  
  # technical coefficient matrix A
  A <- as.matrix(t((1/list_x[[i]]) * t(list_Z[[i]])))
  A[is.na(A)] <- 0 #removing NA value
  A_rr <- matrix(0, 3015, 3015) 
  fatdiag(A_rr, steps = 67) <- fatdiag(A, steps = 67) # creating on-diagonal technical cofficient matrix
  
  # Leontief inverse matrix (merged China and Mexico)
  L <- solve(diag(NROW(A_rr)) - A_rr)
  
  # creating trade flows matrix
  Y <- list_fd[[i]]
  Z <- list_Z[[i]]
  
  fatdiag(Z, steps = 67) <- 0
  colnames(Z) <- country_long
  Z <- t(Z) %>% 
    rowsum(., row.names(.), reorder = FALSE) %>% 
    t()
  
  GO <- Y + Z # trade glows matrix (final demand + intermediate input trade)
  
  # emission account
  eebt <- as.matrix(q * L %*% GO) %>% # q * L * y
    set_rownames(country_long) %>% # setting row names
    rowsum(., row.names(.), reorder = FALSE) # row sum by sector, results in m x m eebt matrix
  
  pbe <- eebt %>% # production-based emissions
    rowSums() %>% # sum over column to get PBE
    data.frame() %>% 
    rename(pbe_eebt = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(pbe_eebt = pbe_eebt / 1.0e+6)
  
  cbe <- eebt %>% # consumption-based emissions
    colSums() %>% # sum over row to get CBE
    data.frame() %>% 
    rename(cbe_eebt = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(cbe_eebt = cbe_eebt / 1.0e+6)
  
  list_eebt[[as.character(i)]] <- pbe %>% 
    left_join(cbe, by = c("country", "n_year"))
  
  fatdiag(eebt, steps = 67) <- 0 # removing on-diagonal values (on-diagonal: domestic-domestic emissions)
  
  eex <- eebt %>% # emissions embodied in exports
    rowSums() %>% # sum over column to get emissions embodied in exports
    data.frame() %>% 
    rename(eex_eebt = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(eex_eebt = eex_eebt / 1.0e+6)
  
  eem <- eebt %>% # emissions embodied in imports
    colSums() %>% # sum over row to get emissions embodied in imports
    data.frame() %>% 
    rename(eem_eebt = 1) %>% 
    bind_cols(country, data.frame(i) %>% rename(n_year = 1), .) %>% 
    mutate(eem_eebt = eem_eebt / 1.0e+6)
  
  list_e_eebt[[as.character("")]] <- eex %>% 
    left_join(eem, by = c("country", "n_year"))
  
  rm(q, A, A_rr, L, eebt, pbe, cbe, eex, eem, GO, Y, Z)
  gc()
  
}

saveRDS(list_e_eebt, "list_e_eebt.rds")
saveRDS(list_eebt, "list_eebt.rds")

## Technology-adjusted EEBT model -----

list_eagm <- list()

for (i in n_year) {
  
  # emission intensity
  q <- ef_oecd %>%
    filter(n == i) %>% 
    select(q) %>% 
    as_vector()
  
  # input coefficient matrix A
  A <- as.matrix(t((1/list_x[[i]]) * t(list_Z[[i]])))
  A[is.na(A)] <- 0   # removing NA value
  A_rr <- matrix(0, 3015, 3015)
  fatdiag(A_rr, steps = 67) <- fatdiag(A, steps = 67) # creating on-diagonal technical cofficient matrix
  
  # Leontief inverse matrix (merged China and Mexico)
  L <- solve(diag(NROW(A_rr)) - A_rr)
  
  # final demand import
  Y <- list_fd[[i]]
  fatdiag(Y, steps = 67) <- 0
  
  Y_mod <- Y %>% 
    set_rownames(sector_long) %>% 
    rowsum(., row.names(.), reorder = FALSE) %>% 
    data.frame() %>% 
    rownames_to_column() %>% 
    rename("sector" = 1) %>% 
    pivot_longer(2:68, names_to = "country", values_to = "y") %>% 
    right_join(rnames, ., by = c("country", "sector")) %>% 
    select(y)
  
  # intermediate import
  Z <- list_Z[[i]]
  fatdiag(Z, steps = 67) <- 0
  
  Z_mod <- Z %>% 
    set_rownames(sector_long) %>% 
    rowsum(., row.names(.), reorder = FALSE) %>% 
    t() %>% 
    set_rownames(country_long) %>% 
    rowsum(., row.names(.), reorder = FALSE) %>% 
    t() %>% 
    data.frame() %>% 
    rownames_to_column() %>% 
    rename("sector" = 1) %>% 
    pivot_longer(2:68, names_to = "country", values_to = "Z") %>% 
    right_join(rnames, ., by = c("country", "sector")) %>% 
    select(Z)
  
  # gross import
  m <- as.matrix(Z_mod + Y_mod)
  
  # emission account
  eagm <- (q * L %*% m) %>% # if the gross imports were produced by domestic tech
    set_rownames(country_long) %>% 
    rowsum(., row.names(.), reorder = FALSE) %>% 
    data.frame() %>% 
    rename(eam_eebt = 1) %>% 
    bind_cols(country, .)  %>% 
    mutate(n_year = i) %>% 
    select(country, n_year, eam_eebt) %>% 
    mutate(eam_eebt = eam_eebt / 1.0e+6)
  
  list_eagm[[as.character("")]] <- eagm
  
  rm(q, A, A_rr, L, Y, Z, Y_mod, Z_mod, m, eagm)
  gc()
  
}

saveRDS(list_eagm, "list_eagm.rds")

eagm <- bind_rows(list_eagm)

pbe_cbe <- bind_rows(list_eebt) 

country_full <- read.xlsx("ICIO_country_full.xlsx") %>% 
  bind_cols(country, .)

eebt <- bind_rows(list_e_eebt) %>% # gathering all results from the EEBT model
  left_join(eagm, by = c("country", "n_year")) %>% 
  left_join(pbe_cbe, by = c("country", "n_year")) %>% 
  left_join(year %>% 
              rename(n_year = 2), by = "n_year") %>% 
  select(country, year, n_year, pbe_eebt, cbe_eebt, eex_eebt, eem_eebt, eam_eebt) %>% 
  left_join(bind_rows(DMHH_pb), by = c("country", "year")) %>% 
  left_join(bind_rows(DMHH_cb), by = c("country", "year")) %>% 
  mutate(pbe_eebt = pbe_eebt + DMHH_pb,
         cbe_eebt = cbe_eebt + DMHH_cb,
         beet_eebt = eex_eebt - eem_eebt,
         NEON_eebt = eex_eebt - eam_eebt,
         apbe_eebt = pbe_eebt - NEON_eebt) %>% 
  left_join(country_full, by = 'country') %>% 
  select(ICIO, country, year, n_year, pbe_eebt, cbe_eebt, eex_eebt, eem_eebt, eam_eebt, beet_eebt, NEON_eebt, apbe_eebt) 

saveRDS(eebt, "eebt.rds")

output <- eebt %>% 
  left_join(mrio, by = c("country", "year", "n_year")) %>% 
  rename(country = 1, isocode = 2)

saveRDS(output, "output.rds")

write.xlsx(output, "output.xlsx")


# Visualization ------
## Scatter Plot -------
library(pwt10)
data(pwt10.0)
pwt10.0 <- data.frame(as.list(pwt10.0), stringsAsFactors = FALSE) %>%
  as_tibble()

i <- sapply(pwt10.0 , is.factor)
pwt10.0 [i] <- lapply(pwt10.0 [i], as.character)
rm(i)

output <- readRDS("output.rds") %>% 
  left_join(pwt10.0 %>% 
              select(isocode, year, rgdpe, pop), by = c("isocode", "year")) %>% 
  mutate(gdp_pc = (rgdpe/pop)/1000) %>% 
  group_by(isocode) %>% 
  mutate(avg_pop = mean(pop, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(avg_pop >= 7.5)



plot_1 <- output %>% 
  filter(year %in% c(1995, 2007, 2018)) %>% 
  mutate(plot = NONS_mrio/((cbe_mrio + pbe_mrio)/2)*100) %>% 
  filter(!grepl("ROW", isocode)) %>% 
  select(isocode, year, gdp_pc, plot) %>%
  distinct_at(c("year", "isocode"), .keep_all = TRUE) %>% 
  mutate(Model = c("Net Onshoring"))

plot_2 <- output %>% 
  filter(year %in% c(1995, 2007, 2018)) %>% 
  mutate(plot = beet_mrio/((cbe_mrio + pbe_mrio)/2)*100) %>% 
  filter(!grepl("ROW", isocode)) %>% 
  select(isocode, year, gdp_pc, plot) %>%
  distinct_at(c("year", "isocode"), .keep_all = TRUE) %>% 
  mutate(Model = c("BEET"))

plot_facet <- bind_rows(plot_1, plot_2) %>% 
  mutate(Model = factor(Model, levels = c("Net Onshoring", "BEET"))) %>% 
  unite("Model", Model, year, sep = " ") %>% 
  mutate(Model = factor(Model, levels = c("Net Onshoring 1995", "BEET 1995", "Net Onshoring 2007", "BEET 2007", "Net Onshoring 2018", "BEET 2018")))


pdf("SctPlotFacet_1_new.pdf", width = 16, height = 24)
print(
  ggplot(plot_facet, aes(gdp_pc, plot)) +
    geom_hline(yintercept=0, size = 0.4, lty = 2) +
    geom_point(shape = 16, size = 3, alpha = 0.7) +
    stat_smooth(geom = "line", method = "lm", se = FALSE, size = 0.4, color = "#00008B", alpha = 0.7) +
    theme_pander(base_family = "Times") +
    geom_text_repel(aes(label = isocode), size = 6, alpha = 0.6, max.overlaps = getOption("ggrepel.max.overlaps", default = 12000)) +
    facet_wrap(~Model, scales = "free", ncol = 2) +
    scale_y_continuous() +
    theme(plot.background = element_rect(fill = "transparent"),
          legend.position = "none",
          plot.margin = margin(1.5, 1, 2, 0.5, "cm"),
          axis.title.x = element_text(size = 0),
          axis.title.y = element_text(size = 0),
          plot.caption = element_text(size = 8, vjust = -4),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          panel.spacing = unit(3, "lines"),
          strip.text.x = element_text(size = 26, vjust = 1.5),
          axis.line=element_line()))

dev.off()

## Line plots -----
# Emissions Embodied in Trade
Data_plot <- output %>% 
  group_by(isocode) %>% 
  mutate(avg_pop = mean(pop, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(avg_pop >= 7.5) %>% 
  select(country, year, eex_mrio,  eem_mrio) %>% 
  rename("EEX MRIO" = 3,
         "EEM MRIO" = 4) %>%
  gather(., "emission", "value", 3:4) %>% 
  separate(emission, into=c("flow", "method"), sep = " ") %>% 
  select(-method) %>% 
  mutate(country = as.character(country),
         flow = factor(flow, levels = c("EEX", "EEM")),
         value = value / 1000)


pdf("Line_EET.pdf", width = 17, height = 22)
print(
  ggplot(Data_plot, aes(year, value, color = flow, shape = flow)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_color_grey(start = 0.2, end = 0.7) +
    scale_fill_grey(start = 0.2, end = 0.7) +
    scale_shape_manual(values = c(1, 4)) +
    guides(colour=guide_legend(override.aes=list(shape=NA, size = 1))) +
    facet_wrap(~country, scales = "free_y", ncol = 5) +
    theme_pander(base_family = "Times") +
    theme(plot.title = element_text(size = 24, hjust = .5, vjust = 2.5),
          plot.subtitle = element_text(size = 20, hjust = .5, vjust = 2.5),
          plot.caption = element_text(vjust = 4.5),
          plot.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.border = element_rect(color = "white"),
          panel.spacing.y = unit(1, "lines"),
          plot.margin = margin(2, 1, 0, 0, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(y=expression(CO[2]* " Emissions (Gt)"))
)

dev.off()


# Emissions Avoided by Imports
Data_plot <- output %>% 
  group_by(isocode) %>% 
  mutate(avg_pop = mean(pop, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(avg_pop >= 7.5) %>% 
  select(country, year, eem_mrio, eam_mrio, eam_eebt) %>% 
  rename("EEM" = 3,
         "EAM-MRIO" = 4,
         "EAM-EEBT" = 5) %>% 
  gather(., "emission", "value", 3:5) %>% 
  mutate(country = as.character(country),
         emission = factor(emission, levels = c("EEM", "EAM-MRIO", "EAM-EEBT")),
         value = value / 1000)

#Custom Potrait
pdf("Line_EAM.pdf", width = 17, height = 22)
print(
  ggplot(Data_plot, aes(year, value, color = emission, shape = emission)) +
    geom_line() +
    geom_point() +
    scale_color_grey(start = 0.2, end = 0.7) +
    scale_fill_grey(start = 0.2, end = 0.7) +
    scale_shape_few() +
    facet_wrap(~country, scales = "free_y", ncol = 5) +
    theme_pander(base_family = "Times") +
    theme(plot.title = element_text(size = 24, hjust = .5, vjust = 2.5),
          plot.subtitle = element_text(size = 20, hjust = .5, vjust = 2.5),
          plot.caption = element_text(vjust = 4.5),
          plot.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.border = element_rect(color = "white"),
          panel.spacing.y = unit(1, "lines"),
          plot.margin = margin(2, 1, 0, 0, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(y=expression(CO[2]* " Emissions (Gt)"))
)

dev.off()

# Balance of Emissions Embodied in Trade
Data_plot <- output %>% 
  group_by(isocode) %>% 
  mutate(avg_pop = mean(pop, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(avg_pop >= 7.5) %>% 
  select(country, year, beet_mrio, NONS_mrio) %>% 
  rename("BEET MRIO" = 3,
         "Net-Onshoring MRIO" = 4) %>% 
  gather(., "emission", "value", 3:4) %>% 
  separate(emission, into=c("treatment", "method"), sep = " ") %>% 
  select(-method) %>% 
  mutate(country = as.character(country),
         treatment = factor(treatment, levels = c("Net-Onshoring", "BEET")),
         value = value / 1000)


#Custom Potrait
pdf("Line_NONS.pdf", width = 17, height = 22)
print(
  ggplot(Data_plot, aes(year, value, color = treatment, shape = treatment)) +
    geom_line() +
    geom_point(size = 1.5) +
    geom_hline(yintercept=0, size = 0.4, lty = 2) +
    scale_color_grey(start = 0.2, end = 0.7) +
    scale_fill_grey(start = 0.2, end = 0.7) +
    scale_shape_manual(values = c(1, 4)) +
    guides(colour=guide_legend(override.aes=list(shape=NA, size = 1))) +
    facet_wrap(~country, scales = "free_y", ncol = 5) +
    theme_pander(base_family = "Times") +
    theme(plot.title = element_text(size = 24, hjust = .5, vjust = 2.5),
          plot.subtitle = element_text(size = 20, hjust = .5, vjust = 2.5),
          plot.caption = element_text(vjust = 4.5),
          plot.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.border = element_rect(color = "white"),
          panel.spacing.y = unit(1, "lines"),
          plot.margin = margin(2, 1, 0, 0, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(y=expression(CO[2]* " Emissions (Gt)"))
)

dev.off()



pdf("Line_NONS_selected.pdf", width = 8, height = 8)
print(
  ggplot(Data_plot %>% 
           filter(country %in% c("China", "Germany", "India", "Japan", "United Kingdom", "United States")), aes(year, value, color = treatment, shape = treatment)) +
    geom_line() +
    geom_point(size = 1.5) +
    geom_hline(yintercept=0, size = 0.4, lty = 2) +
    scale_color_grey(start = 0.2, end = 0.7) +
    scale_fill_grey(start = 0.2, end = 0.7) +
    scale_shape_manual(values = c(1, 4)) +
    guides(colour=guide_legend(override.aes=list(shape = NA, size = 1))) +
    facet_wrap(~country, scales = "free_y", ncol = 2) +
    theme_pander(base_family = "Times") +
    theme(plot.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          strip.text.x = element_text(size = 16, margin = margin(0,0,0.5,0, "cm")),
          panel.border = element_rect(color = "white"),
          panel.spacing.x = unit(1, "lines"),
          panel.spacing.y = unit(2, "lines"),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))
)

dev.off()


Data_plot <- output %>%
  group_by(isocode) %>% 
  mutate(avg_pop = mean(pop, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(avg_pop >= 7.5) %>%
  select(country, year, pbe_mrio, apbe_mrio) %>% 
  rename("PBE_actual" = 3,
         "Domestic-technology CBE_actual" = 4) %>% 
  gather(., "emission", "value", 3:4) %>% 
  separate(emission, into=c("treatment", "method"), sep = "_") %>% 
  select(-method) %>% 
  mutate(country = as.character(country),
         treatment = factor(treatment, levels = c("Domestic-technology CBE", "PBE")),
         value = value / 1000)

#Custom Potrait
pdf("Line_PBE.pdf", width = 17, height = 22)
print(
  ggplot(Data_plot, aes(year, value, color = treatment, shape = treatment)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_color_grey(start = 0.2, end = 0.7) +
    scale_fill_grey(start = 0.2, end = 0.7) +
    scale_shape_manual(values = c(1, 4)) +
    guides(colour=guide_legend(override.aes=list(shape=NA, size = 1))) +
    facet_wrap(~country, scales = "free_y", ncol = 5) +
    theme_pander(base_family = "Times") +
    theme(plot.title = element_text(size = 24, hjust = .5, vjust = 2.5),
          plot.subtitle = element_text(size = 20, hjust = .5, vjust = 2.5),
          plot.caption = element_text(vjust = 4.5),
          plot.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.border = element_rect(color = "white"),
          panel.spacing.y = unit(1, "lines"),
          plot.margin = margin(2, 1, 0, 0, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(y=expression(CO[2]* " Emissions (Gt)"))
)

dev.off()



pdf("Line_PBE_selected.pdf", width = 8, height = 8)
print(
  ggplot(Data_plot %>% 
           filter(country %in% c("China", "Germany", "India", "Japan", "United Kingdom", "United States")), aes(year, value, color = treatment, shape = treatment)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_color_grey(start = 0.2, end = 0.7) +
    scale_fill_grey(start = 0.2, end = 0.7) +
    scale_shape_manual(values = c(1, 4)) +
    guides(colour=guide_legend(override.aes=list(shape = NA, size = 1))) +
    facet_wrap(~country, scales = "free_y", ncol = 2) +
    theme_pander(base_family = "Times") +
    theme(plot.background = element_rect(fill = "transparent"),
          legend.background = element_rect(fill = "transparent"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          strip.text.x = element_text(size = 16, margin = margin(0,0,0.5,0, "cm")),
          panel.border = element_rect(color = "white"),
          panel.spacing.x = unit(1, "lines"),
          panel.spacing.y = unit(2, "lines"),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))
)

dev.off()


## Tables -------------

library("kableExtra")

output <- readRDS("output.rds") %>% 
  left_join(pwt10.0 %>% 
              select(isocode, year, rgdpe, pop), by = c("isocode", "year")) %>% 
  mutate(gdp_pc = (rgdpe/pop)/1000)



# Balance of Emissions Embodied in Trade
dt95 <- output %>%
  filter(year == 1995) %>% 
  select(country, beet_mrio, NONS_eebt, NONS_mrio) %>% 
  set_colnames(c("country", "BEET", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame()

dt07 <- output %>%
  filter(year == 2007) %>% 
  select(beet_mrio, NONS_eebt, NONS_mrio) %>% 
  set_colnames(c("BEET", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

dt18 <- output %>%
  filter(year == 2018) %>% 
  select(beet_mrio, NONS_eebt, NONS_mrio) %>% 
  set_colnames(c("BEET", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

write_data_sel <- bind_cols(dt95 , dt07, dt18) %>%
  arrange(country) %>% 
  mutate_if(is.numeric, function(x, na.rm=TRUE) (x/1000)) %>% 
  set_colnames(c("Country", "BEET", "Adjusted EEBT", "Adjusted MRIO", "BEET", "Adjusted EEBT", "Adjusted MRIO", "BEET", "Adjusted EEBT", "Adjusted MRIO"))

# Saving table in LaTex kable format
print(
  kable(write_data_sel, 
        digits = 3,
        longtable = T,
        caption = "Balance of Emissions Embodied in Trade (Gt of CO2 Emissions)",
        label = "netons-gt",
        format="latex", 
        booktabs=TRUE) %>%
    column_spec(2:10, width = "1cm") %>% 
    kable_styling(font_size = 8) %>%
    add_header_above(c(" " = 1,"1995" = 3, "2007" = 3, "2018" = 3)),
  size = "\\fontsize{8pt}{10pt}\\selectfont",
  floating = TRUE, 
  latex.environments = "center",
  caption.placement = "top",
  math.style.exponents = TRUE,
  file = "Table_NONS_Long.tex"
)



# Balance of Emissions Embodied in Trade (% of PBE)
dt95 <- output %>%
  filter(year == 1995) %>% 
  mutate(NONS_eebt = NONS_eebt/pbe_eebt*100,
         NONS_mrio = NONS_mrio/pbe_mrio*100,
         beet_mrio = beet_mrio/pbe_mrio*100) %>% 
  select(country, beet_mrio, NONS_eebt, NONS_mrio) %>% 
  set_colnames(c("country", "BEET", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame()

dt07 <- output %>%
  filter(year == 2007) %>% 
  mutate(NONS_eebt = NONS_eebt/pbe_eebt*100,
         NONS_mrio = NONS_mrio/pbe_mrio*100,
         beet_mrio = beet_mrio/pbe_mrio*100) %>% 
  select(beet_mrio, NONS_eebt, NONS_mrio) %>% 
  set_colnames(c("BEET", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

dt18 <- output %>%
  filter(year == 2018) %>% 
  mutate(NONS_eebt = NONS_eebt/pbe_eebt*100,
         NONS_mrio = NONS_mrio/pbe_mrio*100,
         beet_mrio = beet_mrio/pbe_mrio*100) %>% 
  select(beet_mrio, NONS_eebt, NONS_mrio) %>% 
  set_colnames(c("BEET", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

write_data_sel2 <- bind_cols(dt95 , dt07, dt18) %>%
  arrange(country) %>% 
  #mutate_if(is.numeric, function(x, na.rm=TRUE) (x/1000)) %>% 
  set_colnames(c("Country", "BEET", "Adjusted EEBT", "Adjusted MRIO", "BEET", "Adjusted EEBT", "Adjusted MRIO", "BEET", "Adjusted EEBT", "Adjusted MRIO"))

# Saving table in LaTex kable format
print(
  kable(write_data_sel, 
        digits = 3,
        longtable = T,
        caption = "Balance of Emissions Embodied in Trade in Percentage of PBE",
        label = "netons-gt",
        format="latex", booktabs=TRUE) %>%
    column_spec(2:10, width = "1cm") %>% 
    kable_styling(font_size = 8) %>%
    add_header_above(c(" " = 1,"1995" = 3, "2007" = 3, "2018" = 3)),
  size = "\\fontsize{8pt}{10pt}\\selectfont",
  floating = TRUE, 
  latex.environments = "center",
  caption.placement = "top",
  math.style.exponents = TRUE,
  file = "Table_NONS_PBE_Long.tex"
)



# Production based emissions
dt95 <- output %>%
  filter(year == 1995) %>% 
  select(country, pbe_mrio, apbe_eebt, apbe_mrio) %>% 
  set_colnames(c("country", "PBE", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame()

dt07 <- output %>%
  filter(year == 2007) %>% 
  select(pbe_mrio, apbe_eebt, apbe_mrio) %>% 
  set_colnames(c("PBE", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

dt18 <- output %>%
  filter(year == 2018) %>% 
  select(pbe_mrio, apbe_eebt, apbe_mrio) %>% 
  set_colnames(c("PBE", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

write_data_sel <- bind_cols(dt95 , dt07, dt18) %>%
  arrange(country) %>% 
  mutate_if(is.numeric, function(x, na.rm=TRUE) (x/1000)) %>% 
  set_colnames(c("Country", "PBE", "Adjusted EEBT", "Adjusted MRIO", "PBE", "Adjusted EEBT", "Adjusted MRIO", "PBE", "Adjusted EEBT", "Adjusted MRIO"))

# Saving table in LaTex kable format
print(
  kable(write_data_sel, 
        digits = 3,
        longtable = T,
        caption = "Production-Based Emissions (Gt of CO2 Emissions)",
        label = "pbe-gt",
        format="latex", 
        booktabs=TRUE) %>%
    column_spec(2:10, width = "1cm") %>% 
    kable_styling(font_size = 8) %>%
    add_header_above(c(" " = 1,"1995" = 3, "2007" = 3, "2018" = 3)),
  size = "\\fontsize{8pt}{10pt}\\selectfont",
  floating = TRUE, 
  latex.environments = "center",
  caption.placement = "top",
  math.style.exponents = TRUE,
  file = "Table_APBE_Long.tex"
)



# Emissions Embodied in Import
dt95 <- output %>%
  filter(year == 1995) %>% 
  select(country, eem_mrio, eam_eebt, eam_mrio) %>% 
  set_colnames(c("country", "EEM", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame()

dt07 <- output %>%
  filter(year == 2007) %>% 
  select(eem_mrio, eam_eebt, eam_mrio) %>% 
  set_colnames(c("EEM", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

dt18 <- output %>%
  filter(year == 2018) %>% 
  select(eem_mrio, eam_eebt, eam_mrio) %>% 
  set_colnames(c("EEM", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame()  

write_data_sel <- bind_cols(dt95 , dt07, dt18) %>%
  arrange(country) %>% 
  mutate_if(is.numeric, function(x, na.rm=TRUE) (x/1000)) %>% 
  set_colnames(c("Country", "EEM", "Adjusted EEBT", "Adjusted MRIO", "EEM", "Adjusted EEBT", "Adjusted MRIO", "EEM", "Adjusted EEBT", "Adjusted MRIO"))

# Saving table in LaTex kable format
print(
  kable(write_data_sel, 
        digits = 3,
        longtable = T,
        caption = "Balance of Emissions Embodied in Import (Gt of CO2 Emissions)",
        label = "eem-gt",
        format="latex", 
        booktabs=TRUE) %>%
    column_spec(2:10, width = "1cm") %>% 
    kable_styling(font_size = 8) %>%
    add_header_above(c(" " = 1,"1995" = 3, "2007" = 3, "2018" = 3)),
  size = "\\fontsize{8pt}{10pt}\\selectfont",
  floating = TRUE, 
  latex.environments = "center",
  caption.placement = "top",
  math.style.exponents = TRUE,
  file = "Table_EAM_Long.tex"
)



# Emissions Embodied in Imports (% of PBE)
dt95 <- output %>%
  filter(year == 1995) %>% 
  mutate(eem_mrio = eem_mrio/pbe_mrio*100,
         eam_eebt = eam_eebt/pbe_eebt*100,
         eam_mrio = eam_mrio/pbe_mrio*100) %>% 
  select(country, eem_mrio, eam_eebt, eam_mrio) %>% 
  set_colnames(c("country", "EEM", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame()

dt07 <- output %>%
  filter(year == 2007) %>% 
  mutate(eem_mrio = eem_mrio/pbe_mrio*100,
         eam_eebt = eam_eebt/pbe_eebt*100,
         eam_mrio = eam_mrio/pbe_mrio*100) %>% 
  select(eem_mrio, eam_eebt, eam_mrio) %>% 
  set_colnames(c("EEM", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

dt18 <- output %>%
  filter(year == 2018) %>% 
  mutate(eem_mrio = eem_mrio/pbe_mrio*100,
         eam_eebt = eam_eebt/pbe_eebt*100,
         eam_mrio = eam_mrio/pbe_mrio*100) %>% 
  select(eem_mrio, eam_eebt, eam_mrio) %>% 
  set_colnames(c("EEM", "Adjusted EEBT", "Adjusted MRIO")) %>% 
  data.frame() 

write_data_sel <- bind_cols(dt95 , dt07, dt18) %>%
  arrange(country) %>% 
  #mutate_if(is.numeric, function(x, na.rm=TRUE) (x/1000)) %>% 
  set_colnames(c("Country", "EEM", "Adjusted EEBT", "Adjusted MRIO", "EEM", "Adjusted EEBT", "Adjusted MRIO", "EEM", "Adjusted EEBT", "Adjusted MRIO"))

# Saving table in LaTex kable format
print(
  kable(write_data_sel, 
        digits = 3,
        longtable = T,
        caption = "Emissions Embodied in Import in Percentage of PBE",
        label = "eem-ppb-gt",
        format="latex", booktabs=TRUE) %>%
    column_spec(2:10, width = "1cm") %>% 
    kable_styling(font_size = 8) %>%
    add_header_above(c(" " = 1,"1995" = 3, "2007" = 3, "2018" = 3)),
  size = "\\fontsize{8pt}{10pt}\\selectfont",
  floating = TRUE, 
  latex.environments = "center",
  caption.placement = "top",
  math.style.exponents = TRUE,
  file = "Table_EAM_PBE_Long.tex"
)


# Regression ----
library(stargazer)
library(plm)

output <- readRDS("output.rds") %>% 
  left_join(pwt10.0 %>% 
              select(isocode, year, rgdpe, pop), by = c("isocode", "year")) %>% 
  mutate(gdp_pc = (rgdpe/pop)/1000) %>% 
  group_by(isocode) %>% 
  mutate(avg_pop = mean(pop, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(avg_pop >= 7.5) %>% 
  filter(year %in% c(1995, 2007, 2018))

ols_data <- output %>% 
  select(isocode, year, pop, NONS_mrio, beet_mrio, pbe_mrio, cbe_mrio, gdp_pc) %>% 
  mutate(NONS_mrio = NONS_mrio/((cbe_mrio + pbe_mrio)/2)*100,
         beet_mrio = beet_mrio/((cbe_mrio + pbe_mrio)/2)*100)

NONS_1 <- lm(NONS_mrio ~ gdp_pc, 
             data = ols_data %>% # stargazer automatically applies heteroskedasticity-robust SE
               filter(year %in% c(1995)))

NONS_2 <- lm(NONS_mrio ~ gdp_pc, 
             data = ols_data %>% # stargazer automatically applies heteroskedasticity-robust SE
               filter(year %in% c(2007)))

NONS_3 <- lm(NONS_mrio ~ gdp_pc, 
             data = ols_data %>% # stargazer automatically applies heteroskedasticity-robust SE
               filter(year %in% c(2018)))

beet_1 <- lm(beet_mrio ~ gdp_pc, 
             data = ols_data %>% # stargazer automatically applies heteroskedasticity-robust SE
               filter(year %in% c(1995)))

beet_2 <- lm(beet_mrio ~ gdp_pc, 
             data = ols_data %>% # stargazer automatically applies heteroskedasticity-robust SE
               filter(year %in% c(2007)))

beet_3 <- lm(beet_mrio ~ gdp_pc, 
             data = ols_data %>% # stargazer automatically applies heteroskedasticity-robust SE
               filter(year %in% c(2018)))

stargazer(NONS_1, NONS_2, NONS_3, beet_1, beet_2, beet_3,
          type = "text",
          title = "Table xx. Cross-country regressions of net onshoring on income per capita",
          column.labels=c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)"),
          covariate.labels = c("Income"),
          omit.stat=c("LL","ser","f", "adj.rsq"),
          column.sep.width = "1pt",
          model.numbers=FALSE)



dt <- read.table("../data/timescaleseparation.csv", sep = " ", header = F)
dtplot = dt %>% pivot_longer(cols = c(x1, x1p), values_to = "valuesx", names_to = "scalex") %>% 
  pivot_longer(cols = c(x2, x2p), values_to = "valuesy", names_to = "scaley")
colnames(dt) = c("t", "x1", "x2", "x1p", "x2p")

x1dat = dt %>% select(t, x1)
x2dat = dt %>% select(t, x2)
x1polydat = dt %>% select(t, x1p)
x2polydat = dt %>% select(t, x2p)

ggplot(dtplot, aes(x= valuesx, y = valuesy))+geom_point(aes(color = scaley))

write.table(x1dat, '../data/x1dat.dat', sep = " ", row.names = F,
            quote = F)
write.table(x2dat, '../data/x2dat.dat', sep = " ", row.names = F,
            quote = F)
write.table(x1polydat, '../data/x1polydat.dat', sep = " ", row.names = F,
            quote = F)
write.table(x2polydat, '../data/x2polydat.dat', sep = " ", row.names = F,
            quote = F)

# Making the Upset plot comparing annotated genes shared/unique among MAGs, generated from OrthoVenn data.

library(UpSetR)
library(ggplot2)


# Inputting data. Data generated from OrthoVenn - reformatting here for UpSet plotting.

expressionInput <- c(Burkholderiaceae = 21, 
                     Limnohabitans_sp1 = 62, 
                     Limnohabitans_sp2 = 28,
                     Pedobacter = 117,
                     Polaromonas = 59,
                     `Burkholderiaceae&Limnohabitans_sp1` = 41,
                     `Burkholderiaceae&Limnohabitans_sp2` = 15,
                     `Burkholderiaceae&Pedobacter` = 59,
                     `Burkholderiaceae&Polaromonas` = 62,
                     `Limnohabitans_sp1&Limnohabitans_sp2` = 641,
                     `Limnohabitans_sp1&Pedobacter` = 35,
                     `Limnohabitans_sp1&Polaromonas` = 133,
                     `Limnohabitans_sp2&Pedobacter` = 12,
                     `Limnohabitans_sp2&Polaromonas` = 69,
                     `Pedobacter&Polaromonas` = 25,
                     `Burkholderiaceae&Limnohabitans_sp1&Limnohabitans_sp2` = 36,
                     `Burkholderiaceae&Limnohabitans_sp1&Pedobacter` = 14,
                     `Burkholderiaceae&Limnohabitans_sp1&Polaromonas` = 56,
                     `Burkholderiaceae&Limnohabitans_sp2&Pedobacter` = 6,
                     `Burkholderiaceae&Limnohabitans_sp2&Polaromonas` = 23,
                     `Burkholderiaceae&Pedobacter&Polaromonas` = 36,
                     `Limnohabitans_sp1&Limnohabitans_sp2&Pedobacter` = 32,
                     `Limnohabitans_sp1&Limnohabitans_sp2&Polaromonas` = 538,
                     `Limnohabitans_sp1&Pedobacter&Polaromonas` = 13,
                     `Limnohabitans_sp2&Pedobacter&Polaromonas` = 5,
                     `Burkholderiaceae&Limnohabitans_sp1&Limnohabitans_sp2&Pedobacter` = 13,
                     `Burkholderiaceae&Limnohabitans_sp1&Limnohabitans_sp2&Polaromonas` = 372,
                     `Burkholderiaceae&Limnohabitans_sp1&Pedobacter&Polaromonas` = 67,
                     `Burkholderiaceae&Limnohabitans_sp2&Pedobacter&Polaromonas` = 16,
                     `Limnohabitans_sp1&Limnohabitans_sp2&Pedobacter&Polaromonas` = 77,
                     `Burkholderiaceae&Limnohabitans_sp1&Limnohabitans_sp2&Pedobacter&Polaromonas` = 325)

# Making figure.
# Text scale order: intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars
upset(fromExpression(expressionInput), order.by = "freq", main.bar.color = "orangered3", matrix.color = "lightsteelblue3", point.size = 2.5, line.size = 1, 
      mainbar.y.label = "MAG Gene Intersections", sets.x.label = "MAG Set Size", 
      text.scale = c(1,1,1,1,1), set_size.angles = 45)
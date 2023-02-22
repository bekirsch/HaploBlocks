args <- commandArgs(trailingOnly = TRUE)

file=args[1]

HaploBlocks = read.csv(file)

png(filename = paste(file, "png", sep = "."))

plot(NA, xlim = c(min(HaploBlocks$bp.start), max(HaploBlocks$bp.end)), ylim = c(0, max(HaploBlocks$sHat)), xlab = "position [bp]", ylab = "inferred s-hat")
for (block in 1:nrow(HaploBlocks)) {
  rect(xleft = HaploBlocks$bp.start[block], ybottom = 0, xright = HaploBlocks$bp.end[block], ytop = HaploBlocks$sHat[block], col = 'black', border = NA)
}

dev.off()

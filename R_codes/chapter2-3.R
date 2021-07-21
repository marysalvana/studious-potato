


source("./pkg-config.R")



N <- 10
n <- N^2
TT <- 10
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

xlim <- c(-1, 1.5)
ylim <- c(-1, 1.5)
xseq <- seq(xlim[1],xlim[2],.05)
yseq <- seq(ylim[1],ylim[2],.05)
nx <- length(xseq)
ny <- length(yseq)
nx1 <- nx-1
ny1 <- ny-1
xygrid <- as.matrix(expand.grid(xseq,yseq))

Ygrid <- c(0.5, 0.5) + (xygrid - c(0.5, 0.5)) * as.numeric(distR_C(xygrid, matrix(c(0.5, 0.5), nrow = 1)))

w <- c(0.5, 0.5)

grid_x <- seq(from = -3, to = 3, length.out = 30)
newxygrid <- expand.grid(grid_x, grid_x) %>% as.matrix()

pdf(file = paste(root, 'Figures/estimation_landmarks.pdf', sep = ''), width = 15, height = 6)
split.screen( rbind(c(0.05,0.99,0.1,0.95), c(0.99,0.99,0.1,0.95)))
split.screen( figs = c( 1, 3 ), screen = 1 ) 

screen(3)
par(pty="s")
par(mai=c(0.3, 0.3, 0.3, 0.3))
eqscplot(range(Ygrid[,1]),range(Ygrid[,2]),type="n",xlab="",ylab="")
mtext(expression(s[x]), line = 3, cex = 2, side = 1)
mtext(expression(s[y]), line = 3, cex = 2, side = 2)
mtext("t = 0", line = 1, cex = 2, side = 3)
for (i in 1:ny) {
	lines(Ygrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx1),], lwd = 0.05)
}
for (j in 1:nx) {
	lines(Ygrid[seq((1+(j-1)),nx*ny,nx),], lwd = 0.05)
}
points(sim_grid_locations,pch=16,cex=0.5,col='red')
points(newxygrid,pch=3,cex=0.5,col=3)

screen(4)
par(pty="s")
par(mai=c(0.3, 0.3, 0.3, 0.3))
eqscplot(range(Ygrid[,1]),range(Ygrid[,2]),type="n",xlab="",ylab="", yaxt = 'n')
mtext(expression(s[x]), line = 3, cex = 2, side = 1)
mtext("t = 1", line = 1, cex = 2, side = 3)
Ygrid_temp <- cbind(Ygrid[, 1] + w[1], Ygrid[, 2] + w[2])
  
for (i in 1:ny) {
	lines(Ygrid_temp[(1+(i-1)*nx):((1+(i-1)*nx)+nx1),], lwd = 0.05)
}
for (j in 1:nx) {
	lines(Ygrid_temp[seq((1+(j-1)),nx*ny,nx),], lwd = 0.05)
}
points(sim_grid_locations,pch=16,cex=0.5,col='red')
points(newxygrid,pch=3,cex=0.5,col=3)

screen(5)
par(pty="s")
par(mai=c(0.3, 0.3, 0.3, 0.3))
eqscplot(range(Ygrid[,1]),range(Ygrid[,2]),type="n",xlab="",ylab="", yaxt = 'n')
mtext("t = 2", line = 1, cex = 2, side = 3)
mtext(expression(s[x]), line = 3, cex = 2, side = 1)
  
Ygrid_temp <- cbind(Ygrid[, 1] + w[1] * 2, Ygrid[, 2] + w[2] * 2)

for (i in 1:ny) {
	lines(Ygrid_temp[(1+(i-1)*nx):((1+(i-1)*nx)+nx1),], lwd = 0.05)
}
for (j in 1:nx) {
	lines(Ygrid_temp[seq((1+(j-1)),nx*ny,nx),], lwd = 0.05)
}
points(sim_grid_locations,pch=16,cex=0.5,col='red')
points(newxygrid,pch=3,cex=0.5,col=3)
close.screen( all=TRUE)
dev.off()


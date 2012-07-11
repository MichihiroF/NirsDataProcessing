#
#@auther = Michihiro Fukuhara
#@date = 2012/07/02
#@version:0.1.2

##Parameter Setting
ch24 <- c("Count","CH1","CH2","CH3","CH4","CH5","CH6","CH7","CH8","CH9",
		"CH10","CH11","CH12","CH13","CH14","CH15","CH16","CH17","CH18",
		"CH19","CH20","CH21","CH22","CH23","CH24","Mark","Time")
ch22 <- c("Count","CH1","CH2","CH3","CH4","CH5","CH6","CH7","CH8","CH9",
		"CH10","CH11","CH12","CH13","CH14","CH15","CH16","CH17","CH18",
		"CH19","CH20","CH21","CH22","Mark","Time")
ch22_con <- c("Count","CH1","CH2","CH3","CH4","CH5","CH6","CH7","CH8","CH9",
		"CH10","CH11","CH12","CH13","CH14","CH15","CH16","CH17","CH18",
		"CH19","CH20","CH21","CH22","Mark","Time","BodyMovement","RemoveMark","PreScan")
ch24_con <- c("Count","CH1","CH2","CH3","CH4","CH5","CH6","CH7","CH8","CH9",
		"CH10","CH11","CH12","CH13","CH14","CH15","CH16","CH17","CH18",
		"CH19","CH20","CH21","CH22","CH23","CH24","Mark","Time","BodyMovement","RemoveMark","PreScan")

#24ch*integral ver.
nirs_dataset1 <- function(filename,dirname){
	file1 <- paste(dirname,filename,sep="/")
	nirsdata <- read.csv(file1,col.names=ch24,skip=40)
}

#22ch*integral ver.
nirs_dataset2 <- function(filename,dirname){
	file1 <- paste(dirname,filename,sep="/")
	nirsdata <- read.csv(file1,col.names=ch22,skip=40)
}

#22ch*continuous ver.
nirs_dataset3 <- function(filename,dirname){
	file1 <- paste(dirname,filename,sep="/")
	nirsdata <- read.csv(file1,col.names=ch22_con,skip=40)
}
#24ch*continuous ver.
nirs_dataset4 <- function(filename,dirname){
	file1 <- paste(dirname,filename,sep="/")
	nirsdata <- read.csv(file1,col.names=ch24_con,skip=40)
}

#
#nirs_analysis
nirs_dataset <- function(
		filename,
		dirname,
		measurement=c("continuous","integral"),
		ch_size=c("22","24"))
{
	if(measurement == "integral"){
		switch(ch_size,
				"24" = nirs_dataset1(filename,dirname),
				"22" = nirs_dataset2(filename,dirname)
				)
	}else{
		switch(ch_size,
				"22" = nirs_dataset3(filename,dirname),
				"24" = nirs_dataset4(filename,dirname)
			)
	}
}

# data is numeric
nirs_freq = 10

nirs_fft <- function(data,visible=TRUE){
	sampling = (length(data))
	n = 0:(sampling-1)
	samplefreq = nirs_freq
	t = n/samplefreq
	f=(n*samplefreq)/sampling
	wave = data
	spec = abs(fft(wave))^2
	if(visible==TRUE){
		par(mfrow = c(2,1))
		plot(t,wave,type="l")
		xmax = samplefreq/2
		plot(f,spec,type="l",col = "navy",xlim=c(0,xmax))
	}
	s_number = round(sampling/2)
	f <- f[1:s_number]
	spec = spec[1:s_number]
	return(list(freq = f,s=s_number,spect = spec))
}

select_data <- function(data,s){
	result <- NULL
	for(i in 1:s){
		result <- append(result,data[length(data)/s*i])	
	}
	return(list(result = result,s_number = s))
}

eeps <- function(){
	dev.copy2eps(file=paste(as.POSIXlt(Sys.time()),"eps",sep="."))
	dev.off()
}
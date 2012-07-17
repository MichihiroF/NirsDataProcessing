#@date = 2012/07/16
#ETG-7100対応
#####実装関数
# nirs_dataset：nirsデータの読み込み
# nirs_fft：フーリエ変換

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

#NIRSデータの読み込み
nirs_dataset <- function(
		filename,
		dirname,
		measurement=c("continuous","integral"),
		ch_size=c("22","24"),
		info = FALSE){
	#ヘッダーの情報を表示
	if(info == TRUE){
		nirs_info(filename,dirname)
	}
	#データの種類選択(22or24ch、integral or continuous)
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

#nirsのヘッダー情報見る
nirs_info <- function(filename,dirname){
	file1 <- paste(dirname,filename,sep="/")
	f <- file(file1,"r")
	for(i in 1:31){
		res <- readLines(con=f,1)
		print(res)
	}
	close(f)
}

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

# フーリエ変換
nirs_freq = 10
nirs_fft <- function(data,visible=FALSE){
	sampling = (length(data))
	n = 0:(sampling-1)
	samplefreq = nirs_freq
	t = n/samplefreq
	f=(n*samplefreq)/sampling
	wave = data
	spec = abs(fft(wave))^2
	#fft結果の表示
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

#データを間引く
select_data <- function(data,s){
	result <- NULL
	for(i in 1:s){
		result <- append(result,data[length(data)/s*i])	
	}
	return(list(result = result,s_number = s))
}

#簡易図保存
eeps <- function(){
	dev.copy2eps(file=paste(as.POSIXlt(Sys.time()),"eps",sep="."))
	dev.off()
}

#@date = 2012/07/23
#ETG-7100対応
#####実装関数
# nirs_dataset：nirsデータの読み込み
# nirs_fft：フーリエ変換
# hamming：ハミング関数を作成
# baseline：ベースライン処理
# mSMA：移動平均処理
######################################################
#パッケージ
library(TTR) #移動平均用

#パラメータ
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
	if (!is.numeric(data)) {
		stop("not numeric !")
	}
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

#ベースライン処理
mbaseline <-function(
		data, #時系列データ
		Pretime, #pretimeの範囲を指定（例：Pretime = 200:300）
		Posttime, #posttimeの範囲を指定
		visible = TRUE #実行結果をプロットするかどうか
		){
	#PretimeとPosttimeのデータをBaselineに格納
	PreData  <- data.frame(Pretime,data[Pretime])
	PostData <- data.frame(Posttime,data[Posttime])
	colnames(PreData) <- c("Time","Data")
	colnames(PostData) <- c("Time","Data")
	Baseline <- rbind(PreData,PostData)
	
	##最小二乗法の実行 (y=a*x+b)
	res = lm(Baseline[,"Data"]~Baseline[,"Time"])
	a <- res$coefficients[2]
	b <- res$coefficients[1]
	x <- Baseline[,"Time"]
	y <- a*x+b
	
	#ベースライン後のデータ作成(yy = a*xx+b)
	x_min <- PreData[1,"Time"]
	y_min <- PostData[length(PostData[,"Time"]),"Time"]
	xx <- c(1:length(data[x_min:y_min]))
	yy <-data[x_min:y_min]- (a*xx+b)
	
	#結果のグラフ化
	if(visible==TRUE){
		par(mfrow=c(3,1))
		#data全体の作図と、Fitting直線の表示（処理する部分のみ）
		plot(data,type="l",ylab="data",xlab="Alldata")
		lines(x,y,col="red")
		#ベースラインの範囲のデータとFitting直線のみ表示
		plot(data[x_min:y_min],type="l",ylab="data",xlab="Baselinedata (Before)")
		lines(x,y,col="red")
		#ベースライン後のデータの出力
		plot(xx,yy,type="l",ylab="data",xlab="Baselinedata (After)")
	}
	return(list(x=xx,y=yy,a=a,b=b))
}

#移動平均処理（NAを省く）
mSMA <- function(data,l){
	sma_data = SMA(data,l)
	res = sma_data[l,length(sma_data)]
}

#タスク期間を色付けして描画
task_plot <- function(data,start,end,y_lim,x_lab,y_lab){
	x0 = start
	x1 = end
	y0 = -1.0
	y1 = 1.0
	bgname = "red"
	plot(data,xlab=x_lab,ylab=y_lab,type="l",ylim=y_lim,xaxt="n",yaxt="n");
	rect(x0,y0,x1,y1,col = bgname,density=10);axis(side=2,at=y_lim-0.1);axis(side=1,at=c(x0,x1));
}

##############窓関数###############
#短径窓を掛ける
rectangularW <- function(data){
	data_length <- length(data)
	result <- NULL
	result <- append(result,0)
	for(i in 2:data_length){
		if(i != data_length){
			result <- append(result,data[i])
		}else{
			result <- append(result,0)
		}
	}
	return(result)
}

#ハミング窓を掛ける
hammingW <- function(data){
	pi = 3.141593
	data_length <- length(data)
	result <- NULL
	for(i in 1:data_length){
		ham <- 0.54 - 0.46*cos((2*pi*i)/(data_length-1))
		result <- append(result,(data[i]*ham))
	}
	return(result)
}

#ハニング窓を掛ける
hanningW <- function(data){
	pi = 3.141593
	data_length <- length(data)
	result <- NULL
	for(i in 1:data_length){
			han <- 0.5 - 0.5*cos((2*pi*i)/(data_length-1))
			result <- append(result,(data[i]*han))
	}
	return(result)
}

#ガウス窓を掛ける
gaussianW <- function(data,m){
	data_length <- length(data)
	result <- NULL
	for(i in 1:data_length){
			gau <- exp(((-2*(m^2))/(data_length-1)^2)*(i-(data_length-1)/2)^2)
			result <- append(result,(data[i]*gau))
	}
	return(result)
}

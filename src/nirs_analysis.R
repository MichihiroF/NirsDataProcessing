#####実装関数
# nirs_dataset()：nirsデータの読み込み(ETG-7100)
# nirs_hot()：データ読み込み(nirs_hot)
# nirs_fft()：フーリエ変換
# mbaseline()：ベースライン処理
# mSMA()：移動平均処理
# task_plot()：作図
# ftest()：nirs用にt検定(select_dataと併せて使うほうが良い)
# select_data()：時系列データをn点ごとに抜き出す
# eeps(),eeps2()：簡易図保存
# overrap_fft()：FFTのオーバーラップ加算平均
# overrap_fft2()：FFTのオーバーラップ処理。加算はせず、全データをはきだす
# mfreq_plot：周波数の時間変化可視化関数

#窓関数色々：短絡、ハミング、ハニング、ガウス、ブラックマンハリス
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

#hot
ch_2 <-c("Time","Bloodflow(Left)","Bloodflow(Left1cm)","Bloodflow(Left13cm)","Bloodflow(Right)","Bloodflow(Right1cm)","Bloodflow(Right3cm)",
		"Pulse(Left)","Pulse(Right)","LF/HF(left)","LF/HF(Right)",
		"temparature","x-axis","y-axis","z-axis","DeepBreathDegree","Chart_x","Chart_y","Chart_syogen","Chart_radius",
		"","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","")

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

#hot
nirs_hot <- function(dirname,filename){
	file1 <- paste(dirname,filename,sep="/")
	read.csv(file1,col.names = ch_2,skip=24)
}

# フーリエ変換
nirs_freq = 10
nirs_fft <- function(data,visible=FALSE,window = FALSE){
	if (!is.numeric(data)) {
		stop("not numeric !")
	}
	if(window == TRUE){
		data = hanningW(data)
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

#FFTのオーバーラップ加算平均
overrap_fft <- function(data,window_size,overrap_rate){
	data_size = length(data)
	shift = floor(window_size*overrap_rate/100)
	#cat("data_size:");print(data_size)
	#cat("Shift:");print(shift)
	count = 0
	top = 0
	while(top <= data_size){
		if(count >= 1){
			bottom = (top+1)-shift
			top = top+window_size-shift
		}
		else{
			bottom = 1
			top = window_size
		}
		#cat("Bottom:");print(bottom)
		#cat("Top   :");print(top)
		count = count+1
	}
	max_count = count-1
	total_spec_data = 0
	#cat("MAXCOUNT:");print(max_count)
	for(i in 1:max_count){
		if(i > 1){
			bottom = (top+1)-shift
			top = top+window_size-shift
		}
		else{
			bottom = 1
			top = window_size
		}
		spec_data = nirs_fft(data[bottom:top],window=TRUE)$spect
		total_spec_data = total_spec_data + spec_data
	}
	res = total_spec_data/max_count
	return(res)
}

#FFTのオーバーラップ処理
overrap_fft2 <- function(data,window_size,overrap_rate){
	data_size = length(data)
	shift = floor(window_size*overrap_rate/100)
	count = 0
	top = 0
	while(top <= data_size){
		if(count >= 1){
			bottom = (top+1)-shift
			top = top+window_size-shift
		}
		else{
			bottom = 1
			top = window_size
		}
		count = count+1
	}
	max_count = count-1
	total_spec_data = 0
	for(i in 1:max_count){
		if(i > 1){
			bottom = (top+1)-shift
			top = top+window_size-shift
		}
		else{
			bottom = 1
			top = window_size
		}
		spec_data = nirs_fft(data[bottom:top],window=TRUE)$spect
		total_spec_data = append(total_spec_data,spec_data)
	}
	return(total_spec_data)
}


#データを間引く
select_data <- function(data,s){
	result <- NULL
	for(i in 1:s){
		result <- append(result,data[length(data)/s*i])	
	}
	return(list(result = result,s_number = s))
}

#簡易図保存1
eeps <- function(){
  ym = paste(as.POSIXlt(Sys.time())$year+1900,as.POSIXlt(Sys.time())$mon,sep="")
  ymd = paste(ym,as.POSIXlt(Sys.time())$mday,sep="")
  ymdh = paste(ymd,as.POSIXlt(Sys.time())$hour,sep="")
  ymdhm = paste(ymdh,as.POSIXlt(Sys.time())$min,sep="")
  ymdhms = paste(ymdhm,round(as.POSIXlt(Sys.time())$sec),sep="")
	dev.copy2eps(file=paste(ymdhms,"eps",sep="."))
	dev.off()
}

#簡易図保存2
eeps2 <- function(filename){
	dev.copy2eps(file=paste(filename,"eps",sep="."))
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
	res = sma_data[l:length(sma_data)]
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

#回帰直線つき散布図
liner_plot <- function(dt){
	a = lm(dt~c(1:length(dt)))$coefficients[2]
	b = lm(dt~c(1:length(dt)))$coefficients[1]
	plot(dt);abline(b,a,col="red")
	return(a)
}

#レストとタスクでt検定
ftest <- function(data,resttime,tasktime,alt){
	rest=data[resttime]
	task=data[tasktime]
	res=t.test(rest,task,alternative =alt)
	return(list(pval = res$p.value,tval = res$statistic))
}

#周波数の時間変化可視化関数
#datは時系列データ、densityは階調値（8あたりが標準）
mfreq_plot <- function(dat,density){
	res = overrap_fft2(dat,256,99)
	x <-1:(length(res)/128)
	y <- nirs_fft(dat[1:256])$freq
	z <- matrix(res[2:length(res)],length(res)/length(x),length(x))
	filled.contour(x,y,t(z),nlevels=density,col=gray((density:0)/density),ylab="Frequency",xlab="FFT count")
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

#ブラックマンハリス窓
blackman <- function(data){
	data_length <- length(data)
	result <- NULL
	for(i in 1:data_length){
		blk <- 0.35875-0.48829*cos(2*pi*i)+0.14128*cos(4*pi*i)-0.01168*cos(6*pi*i)
		result <- append(result,(data[i]*blk))
	}
	return(result)
}
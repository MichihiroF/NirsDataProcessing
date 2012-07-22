#@date = 2012/07/16
#未完成(とりあえず動く)
#nirs_analysis.Rの関数も使用
######################################################

#Install Package
library(psych)
library(MCMCpack)

#引数の説明
#data：NIRSデータ(行：時系列　列：被験者)
#std：標準化の有無
#fft：フーリエ変換の有無
#select：データ削減の有無
#max_Hz：使用する周波数帯域の最大値(最大5Hz)
#s_number：因子分析に使用する行のデータ数
#fmethod：パラメータ推定方法 (ml,pa,gls,wls,minres)
#fa_rotate：回転方法 visible：結果の表示
#visible：結果の表示

mfactanal <- function(
		data,
		std =TRUE,
		fft=TRUE,
		select = TRUE,
		max_Hz = 5,
		s_number,
		fmethod,
		fa_rotate,
		visible=TRUE,
		cluster=TRUE){	
	#被験者数とデータ長の取得
	nnumber <- length(data[1,])
	data_length <- length(data[,1])

	###########標準化###########
	if(std==TRUE){
		print("STD")
		s_data <- NULL
		for(i in 1:nnumber){
			tmp <- scale(data[,i])
			s_data <- append(s_data,tmp)
		}
		s_res <- matrix(s_data,data_length,nnumber)
	}else{
		s_res <- data
	}
	###########フーリエ変換###########	
	if(fft == TRUE){
		print("FFT")
		fft_data <- NULL
		for(i in 1:nnumber){
			tmp <- nirs_fft(s_res[,i],F)$spect
			fft_data <- append(fft_data,tmp)
		}
		fft_res <- matrix(fft_data,data_length/2,nnumber)
	}else{
		fft_res <- s_res
	}

	###########データの削減、間引き###########
	if(select == TRUE){
		print("Select")
		sel_data <- NULL
		for(i in 1:nnumber){
			tmp <-select_data(fft_res[,i][1:round((max_Hz*data_length)/10)],s_number)$result
			sel_data <- append(sel_data,tmp)
		}
		sel_res <- matrix(sel_data,s_number,nnumber)
	}else{
		sel_res <- fft_res
	}
	#行列ラベルの設定
	###########Parallel Analysis
	
	fnumber <- fa.parallel(sel_res,fa="n")$nfact
	
	###########Factor analysis###########
	#fm = ml,pa,gls,wls,minres
	
	factres = fa(sel_res, nfactors=fnumber ,rotate=fa_rotate, fm=fmethod, scores=T)
	#結果Logの出力
	###########結果の表示###########
	if((visible==TRUE) && (fnumber >1)){
		par(ps=24);axis(side=4,labels=F);
		biplot(factres$loading,factres$scores,var.axes = FALSE,expand=1.0,arrow.len = 0.1,xlab="Factor1",ylab="Factor2")
	}
	if((visible==TRUE)&&(fnumber==1)){
		barplot(factres$loading[,1],ylim=c(-0.2,1.0));
		#barplot(factres$loading[,1],ylim=c(-0.2,1.0),yaxt="n");
		#axis(side=2,at=c(-0.2,0.0,1.0));
	}
	if((visible==TRUE) && (fnumber >2)){
		par(mfrow=c(3,1),mex=0.8,ps=25,at=c(-0.2,0.0,1.0))
		barplot(factres$loading[,1],ylim=c(-0.2,1.0),yaxt="n");
		axis(side=2,at=c(-0.2,0.0,1.0));
		barplot(factres$loading[,2],ylim=c(-0.2,1.0),yaxt="n");
		axis(side=2,at=c(-0.2,0.0,1.0));
		barplot(factres$loading[,3],ylim=c(-0.2,1.0),yaxt="n");
		axis(side=2,at=c(-0.2,0.0,1.0));
	}
	if(cluster == TRUE){
		if(fnumber!=1){
			clst = kmeans(factres$loading[,],fnumber)
		}else{
			clst = kmeans(factres$loading[,],2)	
		}
		print("Cluster(kmeans):")
		print(clst)
	}
	return(factres)
}

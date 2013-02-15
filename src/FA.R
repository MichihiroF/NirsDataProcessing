######################################################
#必要なパッケージ
library(psych)
library(MCMCpack)

##引数の説明
#data：NIRSデータ(行：時系列　列：被験者)
#fmethod：パラメータ推定方法 (ml,pa,gls,wls,minres)
#fa_rotate：回転方法 visible：結果の表示
#std：標準化の有無
#dif：時系列データの差分をとる(diff関数を適用する)
#fft：フーリエ変換の有無
#max_Hz：使用する周波数帯域の最大値(最大5Hz)
#window：フーリエ変換の際に窓関数を適用するかどうか(hamming窓)
#select：データ削減の有無
#s_number：因子分析に使用する行のデータ数
#visible：結果の表示(今は因子1と因子2を軸に表示するのみ)
#cluster：因子分析結果のクラスタリング(k-means)

mfactanal <- function(
		data,
		fmethod,
		fa_rotate,
		std =FALSE,
		difference = FALSE,
		fft=FALSE,
		max_Hz = 5,
		window = FALSE,
		select = FALSE,
		s_number,
		visible=TRUE,
		cluster=TRUE){
	#被験者数とデータ長の取得
	nnumber <- length(data[1,])
	data_length <- length(data[,1])
	
	#####ラグ処理#####
	if(difference == TRUE){
		print("DIFF")
		dif_data <- NULL
		for(i in 1:nnumber){
			tmp <- diff(data[,i])
			dif_data <- append(dif_data,tmp)
		}
		dif_res <- matrix(dif_data,data_length,nnumber)
	}else{
		dif_res <- data
	}
	
	###########標準化###########
	if(std==TRUE){
		print("STD")
		s_data <- NULL
		for(i in 1:nnumber){
			tmp <- scale(dif_res[,i])
			s_data <- append(s_data,tmp)
		}
		s_res <- matrix(s_data,data_length,nnumber)
	}else{
		s_res <- dif_res
	}
	
	###########フーリエ変換###########	
	#窓関数の畳込み
	if(window == TRUE){
		win_data <- NULL
		for(i in 1:nnumber){
			tmp <- hanningW(s_res[,i])
			win_data <- append(win_data,tmp)
		}
		win_res <- matrix(win_data,data_length,nnumber)
	}else{
		win_res <- s_res
	}
	
	if(fft == TRUE){
		print("FFT")
		window_size <- 256
		fft_data <- NULL
		for(i in 1:nnumber){
			tmp <- overrap_fft(win_res[,i],window_size,90)
			fft_data <- append(fft_data,tmp)
		}
		fft_res <- matrix(fft_data,window_size/2,nnumber)
	}else{
		fft_res <- win_res
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
	pa_res <- fa.parallel(sel_res,fa="n")
	fnumber <- pa_res$nfact
	
	###########Factor analysis###########
	#fm = ml,pa,gls,wls,minres
	factres = fa(sel_res, nfactors=fnumber ,rotate=fa_rotate, fm=fmethod, scores=T)

	#結果Logの出力(微妙)
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
	if((visible==TRUE) && (fnumber >100000)){
		par(mfrow=c(2,1),mex=0.8,ps=25,at=c(-0.2,0.0,1.0))
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
	return(list(factres=factres,cluster=clst,dif_res=dif_res,pa = pa_res,s_res=s_res,win_res=win_res,fft_res=fft_res))
}

# KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy 
kmo = function( data ){ 
	
	library(MASS) 
	X <- cor(as.matrix(data)) 
	iX <- ginv(X) 
	S2 <- diag(diag((iX^-1))) 
	AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix 
	IS <- X+AIS-2*S2                         # image covariance matrix 
	Dai <- sqrt(diag(diag(AIS))) 
	IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix 
	AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix 
	a <- apply((AIR - diag(diag(AIR)))^2, 2, sum) 
	AA <- sum(a) 
	b <- apply((X - diag(nrow(X)))^2, 2, sum) 
	BB <- sum(b) 
	MSA <- b/(b+a)                        # indiv. measures of sampling adequacy 
	
	AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the 
	# correlation matrix. That is the 
	# negative of the partial correlations, 
	# partialling out all other variables. 
	
	kmo <- BB/(AA+BB)                     # overall KMO statistic 
	
	# Reporting the conclusion 
	if (kmo >= 0.00 && kmo < 0.50){ 
		test <- 'The KMO test yields a degree of common variance 
				unacceptable for FA.' 
	} else if (kmo >= 0.50 && kmo < 0.60){ 
		test <- 'The KMO test yields a degree of common variance miserable.' 
	} else if (kmo >= 0.60 && kmo < 0.70){ 
		test <- 'The KMO test yields a degree of common variance mediocre.' 
	} else if (kmo >= 0.70 && kmo < 0.80){ 
		test <- 'The KMO test yields a degree of common variance middling.' 
	} else if (kmo >= 0.80 && kmo < 0.90){ 
		test <- 'The KMO test yields a degree of common variance meritorious.' 
	} else { 
		test <- 'The KMO test yields a degree of common variance marvelous.' 
	} 
	
	ans <- list(  overall = kmo, 
			report = test, 
			individual = MSA, 
			AIS = AIS, 
			AIR = AIR ) 
	return(ans) 
	
}    # end of kmo() 

comp_est_Method_FA<- function(dt,rotate){
	fa_n= fa.parallel(dt,fa="fa")$nfact
	fa_rotate = rotate
	#最尤法(psych)(ML)
	factMLres = fa(dt, nfactors=fa_n ,rotate=fa_rotate, fm="ml", scores=T)
	#主因子法(PFA)
	factPAres = fa(dt, nfactors=fa_n ,rotate=fa_rotate, fm="pa", scores=T)
	#一般化最小二乗法
	factGLSres <- fa(dt, nfactors=fa_n ,rotate=fa_rotate, fm="gls", scores=T)
	#重みつき最小二乗法
	factWLSres <- fa(dt, nfactors=fa_n ,rotate=fa_rotate, fm="wls", scores=T)
	#最小残差法（重みづけない最小二乗法）
	factOLSres <- fa(dt, nfactors=fa_n ,rotate=fa_rotate, fm="minres", scores=T)
	if(fa_n != 1){
		par(mfrow=c(3,2))
		biplot(factMLres$scores,factMLres$loading,var.axes = TRUE,expand=1.0,arrow.len = 0.1,xlab="ML")
		biplot(factPAres$scores,factPAres$loading,var.axes = TRUE,expand=1.0,arrow.len = 0.1,xlab="PFA")
		biplot(factGLSres$scores,factGLSres$loading,var.axes = TRUE,expand=1.0,arrow.len = 0.1,xlab="GLS")
		biplot(factWLSres$scores,factWLSres$loading,var.axes = TRUE,expand=1.0,arrow.len = 0.1,xlab="WLS")
		biplot(factOLSres$scores,factOLSres$loading,var.axes = TRUE,expand=1.0,arrow.len = 0.1,xlab="MinRes")
	}else{
		par(mfrow=c(3,2))
		barplot(factMLres$loading[,1],ylim=c(-0.2,1.0),yaxt="n");axis(side=2,at=c(-0.4,0.0,1.0));
		barplot(factPAres$loading[,1],ylim=c(-0.2,1.0),yaxt="n");axis(side=2,at=c(-0.4,0.0,1.0));
		barplot(factGLSres$loading[,1],ylim=c(-0.2,1.0),yaxt="n");axis(side=2,at=c(-0.4,0.0,1.0));
		barplot(factWLSres$loading[,1],ylim=c(-0.2,1.0),yaxt="n");axis(side=2,at=c(-0.4,0.0,1.0));
		barplot(factOLSres$loading[,1],ylim=c(-0.2,1.0),yaxt="n");axis(side=2,at=c(-0.4,0.0,1.0));
	}
}
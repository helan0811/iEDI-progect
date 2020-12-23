################## Clusting for image features
library(mclust)

data=read.csv("filename.csv",header=F)
Mdata=matrix(data=0,nrow=1,ncol=9)
data[is.na(data)]=0
Mdata[i,1]=as.numeric(name1[i])
	Mdata[i,2]=dim(data)[2]
	if (Mdata[i,2]<2)
	{
		Mdata[i,3:11]=1
	}
	else{
	for (j in 1:9)
	{
		if (((Mdata[i,2]==2)&((rowMeans(data[j,]))==data[j,1]))|((rowMeans(data[j,]))==data[j,1]))
		{
			Mdata[i,j+2]=1
		}
		else {
		M=Mclust(data[j,], G=1:5)
		Mdata[i,j+2]=M$G
		}
	}
	}

write.csv(Mdata,'EDI_feature.csv')


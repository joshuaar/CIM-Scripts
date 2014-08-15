require(e1071)
require(sva)

trainTestSplit = function(data,classes){
	classIndexes = list()
	# gather classes into a dictionary-like datastructure
	for(i in 1:length(classes)){
		currentClass = classes[i]
		currentClassList = classIndexes[[currentClass]]
		if(is.null(currentClassList))
			currentClassList = c()
		currentClassList = c(i,currentClassList)
		classIndexes[[currentClass]] = currentClassList
	}
	# select one from each class
	testSetIndexes = sapply(classIndexes, function(i)sample(i,1))
	XTest = data[,testSetIndexes]
	XTrain = data[,-testSetIndexes]
	classesTest = classes[testSetIndexes]
	classesTrain = classes[-testSetIndexes]
	out = list()
	out$train_data = XTrain
	out$test_data = XTest
	out$train_classes = classesTrain
	out$test_classes = classesTest
	out
}

crossValIteration = function(data, classes, p, p.adjust=F, shuffle=F){
	# Perform "leave one from each class out" cross validation
	# Does SVA adjustment
	# data: an n by m matrix where columns are samples and rows are features
	# classes: a length m vector of class labels
	# p: P value cutoff
	# p.adjust: adjust p value method. 
	#	FALSE: Do not do p adjustment
	#	"BH": Benjamini and Hochberg adjustment
	#	check documentation for p.adjust for more methods
	# shuffle: Shuffle the class labels randomly without replation
	splittedData = trainTestSplit(data,classes)# list with train_data, test_data, train_classes, test_classes fields
	
	XTrain = splittedData$train_data	
	XTest = splittedData$test_data	
	YTrain = as.factor(splittedData$train_classes)
	YTest = as.factor(splittedData$test_classes)	
	
	YTrain.model = model.matrix(~YTrain)
	YTrain.model0 = model.matrix(~1,data=YTrain)
	trainSv = sva(XTrain,YTrain.model,YTrain.model0)
	
	fsvaobj = fsva(dbdat=XTrain, mod=YTrain.model, sv=trainSv, newdat=XTest)
	XTrain_adjusted = fsvaobj$db
}	

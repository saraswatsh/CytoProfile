#########################################################################################################
# Function to generate a pdf file of plots of pls-da, components, and VIP score of components.
# Author: Shubh Saraswat
# Arguments:
#   x.df: a matrix or data frame with groups, stimulation, and continuous variables
#   colors:  option to have specific colors representing different groups
#   ellipse: option to have ellipse drawn in the classification figure. Default set to false.
#   bg: option to have a prediction background drawn in the classification figure. Default set to false.
#   conf.mat: option to print confusion matrix to show classification from model. Default set to false.
#########################################################################################################

#' Analyze data with PLS-DA method.
#'
#' This function takes a data frame assuming that the first two columns consists
#' of groups and stimulation of patients and using the rest of the columns as
#' components for the PLS-DA analysis.
#'
#' @param x.df A matrix or data frame of variables.
#' @param colors Vector of colors to be set, list of colors to be set equal to the number of groups or treatments. Default set to NULL and will result in random colors.
#' @param title Title of the PDF file to be saved.
#' @param ellipse If ellipse should be drawn in the figures or not. Default set to FALSE.
#' @param bg If prediction background should be drawn in the figures or not. Default set to FALSE.
#' @param conf.mat If confusion matrix of the classifications be printed at the end or not. Default set to FALSE.
#' @description
#' This function takes a matrix or data frame of variables where the first two columns must consist of "Group"
#' and "Stimulation" or "Group" and "Treatment". The rest of columns will be cytokines that will be analyzed.
#' The function also uses log2 transformation of the data, so the data that is inputed, must be raw data without
#' any transformation.
#'
#' @return A PDF file consisting of classification figures, component figures including Variable of importance in projection (VIP) scores, and
#' classification of groups with VIP scores greater than 1.
#' @examples
#' cyt.plsda(cytdata.df[,-c(1,3)], colors = c("black", "purple"), title = "Example Analysis.pdf", bg = TRUE, conf.mat = TRUE)
#' cyt.plsda(cytdata.df[,-c(1,3)], colors = c("black", "purple"), title = "Example Analysis.pdf", ellipse = TRUE, conf.mat = TRUE)
#' cyt.plsda(cytdata.df[,-c(1,4)], ellipse = TRUE, title = "Example Analysis.pdf", conf.mat = TRUE)
#' @export

cyt.plsda = function(x.df, colors = NULL, title, ellipse = FALSE, bg = FALSE, conf.mat = FALSE){
  # Log2 transforming cytokines
  x.df = data.frame(x.df[,c(1:2)], log2(x.df[, -c(1:2)]))

  # Making the first two columns to be lowercase
  names(x.df)[1:2] <- tolower(names(x.df)[1:2])
  # Creating a table to have two separate vectors for group and stimulation
  a = table( x.df[, c(1,2)] ); a

  if("treatment" %in% names(x.df)[1:2]){
    Treatment.vec = dimnames(a)$treatment; Treatment.vec
  }else{
    Stimulation.vec = dimnames(a)$stimulation; Stimulation.vec
  }
  Group.vec = dimnames(a)$group; Group.vec

  # Generate a color palette based on the number of groups
  if (is.null(colors)) {
    num_groups <- length(unique(x.df$group))
    colors <- rainbow(num_groups)
  }

  if("treatment" %in% names(x.df)[1:2]){
    pdf( file= title)
    for(i in 1:length(Treatment.vec) ) {
      theTrt = Treatment.vec[i]
      condt  = x.df[, "treatment"] == theTrt
      Title  = theTrt
      theData.df = x.df[condt,-c(1:2)]
      theGroups  = x.df[condt, "group"]

      cytokine.splsda =
        splsda(theData.df, theGroups, scale=TRUE, ncomp=2, keepX = c(ncol(x.df)-2, ncol(x.df)-2))

      splsda.predict = predict( cytokine.splsda, theData.df, dist="max.dist")
      Prediction1 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy1 = sum(Prediction1[,1] == Prediction1[,2])/length(Prediction1[,1]); accuracy1
      acc1 = 100*signif(accuracy1, digits = 2)

      bg.maxdist <-  background.predict(cytokine.splsda,
                                        comp.predicted = 2,
                                        dist = 'max.dist')
      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }

      plotLoadings(cytokine.splsda, comp=1, contrib = 'max', method = 'mean',
                   legend.color=colors, title=paste("Component 1: ", Title), size.title=1 )
      plotLoadings(cytokine.splsda, comp=2, contrib = 'max', method = 'mean',
                   legend.color=colors, title=paste("Component 2: ", Title), size.title=1 )
      loading.vip.df = data.frame(cytokine.splsda$loadings$X, vip(cytokine.splsda))
      colnames(loading.vip.df)= c("loading.comp1", "loading.comp2", "vip.comp1", "vip.comp2")
      loading.vip.df

      #####VIPscore
      # component 1
      VIPscore=vip(cytokine.splsda)
      Vscore=as.data.frame(VIPscore)
      Vscore$metabo= rownames(Vscore)
      Vscore$comp = Vscore$comp1
      bar=Vscore[,c('metabo','comp')]
      rownames(bar)=c(1:nrow(bar))
      bar=bar[order(bar$comp, decreasing=T),]

      a1 = ggplot(bar, aes(x=metabo, y=comp)) +
        geom_bar(stat="identity", position = "dodge") +scale_y_continuous(limits=c(0,max(bar[,2])))+
        geom_hline(yintercept=1, color="grey") +     #draw a horizontal line
        scale_x_discrete(limits=factor(bar[,1])) +   #set the first column to be a factor
        theme(axis.text.x=element_text(angle=45, hjust=1)) +    # angle
        labs(x="", y="VIP score")+
        ggtitle("Component 1")+
        theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
      print(a1)

      # component 2
      VIPscore=vip(cytokine.splsda)
      Vscore=as.data.frame(VIPscore)
      Vscore$metabo= rownames(Vscore)
      Vscore$comp = Vscore$comp2
      bar=Vscore[,c('metabo','comp')]
      rownames(bar)=c(1:nrow(bar))
      bar=bar[order(bar$comp, decreasing=T),]


      a2=ggplot(bar, aes(x=metabo, y=comp)) +
        geom_bar(stat="identity", position = "dodge") +scale_y_continuous(limits=c(0,max(bar[,2])))+
        geom_hline(yintercept=1, color="grey") +     #draw a horizontal line
        scale_x_discrete(limits=factor(bar[,1])) +   #set the first column to be a factor
        theme(axis.text.x=element_text(angle=45, hjust=1)) +    # angle
        labs(x="", y="VIP score")+
        ggtitle("Component 2")+
        theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
      print(a2)

      # conduct PLSDA on selected variable based on the first component
      condtVariable = Vscore$comp1 > 1
      KeepX = sum( condtVariable )
      KeepX
      theData.mat = theData.df[, condtVariable]
      cytokine.splsda =
        splsda(theData.mat, theGroups, scale=TRUE, ncomp=2, keepX = c(KeepX, KeepX))


      splsda.predict = predict( cytokine.splsda, theData.mat, dist="max.dist")
      Prediction2 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy2 = sum(Prediction2[,1] == Prediction2[,2])/length(Prediction2[,1]); accuracy2
      acc2 = 100*signif(accuracy2, digits = 2)

      bg.maxdist <-  background.predict(cytokine.splsda,
                                        comp.predicted = 2,
                                        dist = 'max.dist')
      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE) {
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title=paste(Title, "(VIP>1)","With Accuracy", acc2, "%" ))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title=paste(Title, "(VIP>1)","With Accuracy", acc2, "%" ),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"))
      }

      plotLoadings(cytokine.splsda, comp=1, contrib = 'max', method = 'mean',
                   legend.color=colors,
                   title=paste("Component 1: ", Title, "(VIP>1)"), size.title=1 )
    }
    dev.off()
    # Prints confusion matrix
    if(conf.mat == TRUE){
      cat("Confusion Matrix for PLS-DA Comparison \n")
      print(get.confusion_matrix(truth = Prediction1[,1], predicted = Prediction1[,2])) # Confusion matrix for all variables in model
      cat("Confusion Matrix for PLS-DA Comparison with VIP Score > 1 \n")
      print(get.confusion_matrix(truth = Prediction2[,1], predicted = Prediction2[,2])) # Confusion matrix for all variables with VIP Score > 1
    }
  }else{
    pdf( file= title)
    for(i in 1:length(Stimulation.vec) ) {
      theTrts = Stimulation.vec[i]
      condt  = x.df[, "stimulation"] == theTrts
      Title  = theTrts
      theData.df = x.df[condt,-c(1:2)]
      theGroups  = x.df[condt, "group"]

      cytokine.splsda =
        splsda(theData.df, theGroups, scale=TRUE, ncomp=2, keepX = c(ncol(x.df)-2, ncol(x.df)-2))

      splsda.predict = predict( cytokine.splsda, theData.df, dist="max.dist")
      Prediction1 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy1 = sum(Prediction1[,1] == Prediction1[,2])/length(Prediction1[,1]); accuracy1
      acc1 = 100*signif(accuracy1, digits = 2)

      bg.maxdist <-  background.predict(cytokine.splsda,
                                        comp.predicted = 2,
                                        dist = 'max.dist')
      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title= paste (Title, "With Accuracy:",acc1,"%"))
      }

      plotLoadings(cytokine.splsda, comp=1, contrib = 'max', method = 'mean',
                   legend.color=colors, title=paste("Component 1: ", Title), size.title=1 )
      plotLoadings(cytokine.splsda, comp=2, contrib = 'max', method = 'mean',
                   legend.color=colors, title=paste("Component 2: ", Title), size.title=1 )
      loading.vip.df = data.frame(cytokine.splsda$loadings$X, vip(cytokine.splsda))
      colnames(loading.vip.df)= c("loading.comp1", "loading.comp2", "vip.comp1", "vip.comp2")
      loading.vip.df

      #####VIPscore
      # component 1
      VIPscore=vip(cytokine.splsda)
      Vscore=as.data.frame(VIPscore)
      Vscore$metabo= rownames(Vscore)
      Vscore$comp = Vscore$comp1
      bar=Vscore[,c('metabo','comp')]
      rownames(bar)=c(1:nrow(bar))
      bar=bar[order(bar$comp, decreasing=T),]

      a1 = ggplot(bar, aes(x=metabo, y=comp)) +
        geom_bar(stat="identity", position = "dodge") +scale_y_continuous(limits=c(0,max(bar[,2])))+
        geom_hline(yintercept=1, color="grey") +     #draw a horizontal line
        scale_x_discrete(limits=factor(bar[,1])) +   #set the first column to be a factor
        theme(axis.text.x=element_text(angle=45, hjust=1)) +    # angle
        labs(x="", y="VIP score")+
        ggtitle("Component 1")+
        theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
      print(a1)

      # component 2
      VIPscore=vip(cytokine.splsda)
      Vscore=as.data.frame(VIPscore)
      Vscore$metabo= rownames(Vscore)
      Vscore$comp = Vscore$comp2
      bar=Vscore[,c('metabo','comp')]
      rownames(bar)=c(1:nrow(bar))
      bar=bar[order(bar$comp, decreasing=T),]


      a2=ggplot(bar, aes(x=metabo, y=comp)) +
        geom_bar(stat="identity", position = "dodge") +scale_y_continuous(limits=c(0,max(bar[,2])))+
        geom_hline(yintercept=1, color="grey") +     #draw a horizontal line
        scale_x_discrete(limits=factor(bar[,1])) +   #set the first column to be a factor
        theme(axis.text.x=element_text(angle=45, hjust=1)) +    # angle
        labs(x="", y="VIP score")+
        ggtitle("Component 2")+
        theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
      print(a2)

      # conduct PLSDA on selected variable based on the first component
      condtVariable = Vscore$comp1 > 1
      KeepX = sum( condtVariable )
      KeepX
      theData.mat = theData.df[, condtVariable]
      cytokine.splsda =
        splsda(theData.mat, theGroups, scale=TRUE, ncomp=2, keepX = c(KeepX, KeepX))


      splsda.predict = predict( cytokine.splsda, theData.mat, dist="max.dist")
      Prediction2 = cbind( original = theGroups, splsda.predict$class$max.dist )
      accuracy2 = sum(Prediction2[,1] == Prediction2[,2])/length(Prediction2[,1]); accuracy2
      acc2 = 100*signif(accuracy2, digits = 2)

      bg.maxdist <-  background.predict(cytokine.splsda,
                                        comp.predicted = 2,
                                        dist = 'max.dist')
      # Conditions to have either ellipses or prediction background or both or neither on graphs.
      if(ellipse == TRUE & bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == TRUE) {
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=TRUE, title=paste(Title, "(VIP>1)", "With Accuracy:",acc2,"%"))
      }
      else if(bg == TRUE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title=paste(Title, "(VIP>1)", "With Accuracy:",acc2,"%"),
                  background = bg.maxdist)
      }
      else if(ellipse == FALSE & bg == FALSE){
        plotIndiv(cytokine.splsda, ind.names=NA, legend=TRUE, col=colors, pch=16,
                  ellipse=FALSE, title= paste (Title, "(VIP>1)", "With Accuracy:",acc2,"%"))
      }

      plotLoadings(cytokine.splsda, comp=1, contrib = 'max', method = 'mean',
                   legend.color=colors,
                   title=paste("Component 1: ", Title, "(VIP>1)"), size.title=1 )
    }
    dev.off()
    # Prints confusion matrix
    if(conf.mat == TRUE){
      cat("Confusion Matrix for PLS-DA Comparison \n")
      print(get.confusion_matrix(truth = Prediction1[,1], predicted = Prediction1[,2])) # Confusion matrix for all variables in model
      cat("Confusion Matrix for PLS-DA Comparison with VIP Score > 1 \n")
      print(get.confusion_matrix(truth = Prediction2[,1], predicted = Prediction2[,2])) # Confusion matrix for all variables with VIP Score > 1
    }
  }
}

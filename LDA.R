# Latend Dirichlet Allocation using (collapsed) Gibbs Sampling.
#
# All you need is provide the input file, and set few parameters.
#
# Please visit https://umdrive.memphis.edu/rbanjade/public/ for updates (if any) and further documentation
# about the code.
#
# Rajendra Banjade
# Dec 2015
#
# Main reference:
#  https://www.youtube.com/watch?v=DDq3OVp9dNA
#  http://u.cs.biu.ac.il/~89-680/darling-lda.pdf 
##############################

# using hash package (https://cran.r-project.org/web/packages/hash/hash.pdf)
library(hash)

#dataFile <- "c:/Users/Rajendra/workspace/RCodes/LDAUsingGibbsSampling/testdata.txt"
dataFile <- "c:/Users/Rajendra/workspace/RCodes/LDAUsingGibbsSampling/msr_train_lda_lemmatized.txt"


####### Hyper parameters. User can set different values.#############
# number of topics to infer. We don't know how many topics are there. 
# We just estimate and try using different number of topics.
K = 5 

#
# Parameter for Dirichlet prior to be used to find document distribution over topics.
# Alpha is a K dimensional vector. However, we set uniform values as we don't know the density.
# < 1 is set, so that the distribution becomes sparse. Smaller value means more sparse.
alpha = 0.001
#alphas[1:K] = rep(alpha, times = K) 
# Parameter for Dirichlet prior used to find topic distribution over words.
beta = 0.001
#betas[1:V] = rep(beta, times = V) we don't know V yet.

#
# Number of iterations to run.
# In one iteration, we go through all documents and all words.
niters = 100


##### Read documents, and initialze various data structures ########## 
#read data <one document in a line> and <tokens separated by one or more spaces or tabs>
conn <- file(dataFile, open="r")
lines <-readLines(conn)
close(conn)

##
# D - total number of documents.
# N - total number of words (in all documents). same word appearing in different places
#     within the document or in another documents are counted separately.
# V - total number of unique words (i.e. size of vocabulary). Valid word or invalid word
#      space separated token is considered as word. 
# NOTE: words are case sensitive. If you want, preprocess your input file separately.
#
D <- length(lines)
N = 0
V = 0
# actual words and their topic assignments. 
documents <- list(D)
wordTopics <- list(D) # size of D x size of documents (variable size). The matrix is zagged.

# vocabulary hash, each entry is: <word, numeric id>
dict <- hash()  # forward dictionary <word, id> 
rdict <- hash() # reverse dictionary <id, word> 

# read tokens (using tokens and words interchangeably).
# d - document id. 
for (d in 1:length(lines)){
  doc <- lines[d]
  # remove punctuations before tokenizing the texts.
  doc <- gsub('[[:punct:]]','',doc)
  #print(line)
  docWords <- strsplit(doc, "\\s+")[[1]]
  N <- N + length(docWords)
  docWordIds <- list()
  #update vocabulary hash table <word, id>
  for(wIdx in 1: length(docWords)) {
    word <- docWords[[wIdx]]
    sprintf("word: %s", word)
    wordId <- dict[[word]]
    #if is not present in the hashtable before.. add it.
    if (is.null(wordId)) {
      # increase vocabulary id.
      V <- V + 1
      .set(dict, keys=word, values=V)
      wordId <- V
      # set word in reverse dictionary as well
      .set(rdict, keys=V, values=word)
    }
    # replacing actual word by its numeric id. 
    # We can recover the actual word by using reverse dictionary <id, word>
    docWordIds[[wIdx]] <- wordId
  }
  documents[[d]] <- docWordIds
}
sprintf("Total number of words: %i, Vocabulary size (unique words): %i", N, V)


# release memory used to store documents as lines
lines <- ""; 

#
# number of words in each document assigned to topic k = 1:K
# its a fixed size matrix D x K
ndkMatrix <- matrix(0,nrow=D, ncol=K)
#number of times word w (unique) is assigned to topic k=1:K
nkwMatrix <- matrix(0,nrow=V, ncol=K)
# number of words assigned to topic k=1:K, one dimensional array.
nk <- rep(0,K)



### randomly initialize topic for each word in each doument; [1,K].
# d - document index
# w - word index in the document d
for (d in 1:length(documents)){
  docWords <- documents[[d]]
  docWordTopics <- list(length(docWords)) 
  for (wd in 1: length(docWords)) {
    # randomly set a topic (1:K)
    wordId <- docWords[[wd]]
    topic <- sample(K, 1)
    docWordTopics[[wd]] <- topic
    ndkMatrix[d,topic] <- ndkMatrix[d,topic] + 1
    nkwMatrix[wordId, topic] <- nkwMatrix[wordId, topic] + 1
    nk[topic] <- nk[topic] + 1
  }
  wordTopics[[d]] <- docWordTopics
}


#
# iteratations of gibbs sampling
for(it in 1:niters) {
  sprintf("Iteration: %i, out of %i", it, niters)
  # iterate through all documents and all words in them.
  for (d in 1:length(documents)){
    # get list of words and their current topic assignments in document d.
    docWords <- documents[[d]]
    docWordTopics <- wordTopics[[d]] 
    
    # iterate through all words in the document
    for (wdIdx in 1: length(docWords)) {
      wordId <- docWords[[wdIdx]]
      wordTopic <- docWordTopics[[wdIdx]]
      
      # reduce the count as we are going to assin the topic 
      # and only depend on assignment of topic to all other words.
      ndkMatrix[d, wordTopic] <- ndkMatrix[d, wordTopic] - 1
      nkwMatrix[wordId, wordTopic] <- nkwMatrix[wordId, wordTopic] - 1
      nk[wordTopic] <- nk[wordTopic] - 1
      
      # find the probability of topic k generating the word w, it will be multi nomial
      multkw <- rep(0.0, K)
      for (k in 1:K) {
        multkw[k] <- (ndkMatrix[d,k] + alpha)*(nkwMatrix[wordId,k]+ beta)/(nk[k] + beta*V) 
        #multkw[k] <- 0.3
      }
      # sample new topic [1,K] from multinomial distribution and update the topic assignment for current word.
      wordNewTopic <- sample(K,size=1, prob = multkw)
      docWordTopics[[wdIdx]] <- wordNewTopic
      
      # increment count based on newly assigned topic.
      ndkMatrix[d, wordNewTopic] <- ndkMatrix[d, wordNewTopic] + 1
      nkwMatrix[wordId, wordNewTopic] <- nkwMatrix[wordId, wordNewTopic] + 1
      nk[wordNewTopic] <- nk[wordNewTopic] + 1
    }
    wordTopics[[d]] <- docWordTopics
  } 
} #gibbs sampling iteration loop ends here.


#
############## calculate thetas and phis ##################################
phi <- matrix(nrow=K, ncol = V)
#thetas <- matrix(nrow=D, ncol=K)

# calculate phi
for(k in 1:K) {
  # calculate normalization factor.
  nzw <- 0
  for (v in 1:V) {
    nzw <- nzw + nkwMatrix[v, k] # + beta; we can add beta at once outside loop as all beta are same 
  }
  nzw <- nzw + V*beta
  #calculate phi for each word for the given topic.
  for(v in 1:V) {
    phi[k,v] <- (nkwMatrix[v, k] + beta)/nzw
  }
}

############# print topics words (top x words out of V) ##################################.
# Each topic is a distribution over words. We can sort them by their probability
# and print only top words. eg. top 20 words.
# get the words in the descending order of their probabilities of being drawn from that topic.
#order(dt[,2, drop=FALSE], decreasing=FALSE)
displayWordsPerTopic <- 10
finalWordTopic <- list()
finalTopicWords <- list()
for (k in 5:5) {
  # get the index of words in V in the decreasing order of their probabilities of being in topic k
  wordIndices <- order(phi[k,, drop=FALSE], decreasing=TRUE)
  #sprintf("Topic id: %i",k)
  paste("Topic..")
  for(idx in 1:displayWordsPerTopic) {
    wordId <- wordIndices[[idx]]
    prob <- phi[k, wordId]
    #sprintf("word and probability: %i %.4f", wordId, prob)
    finalWordTopic[[idx]]<- prob
    finalTopicWords[[idx]] <- rdict[[as.character(wordId)]]
  }
}

#fileConn<-file("output.txt")
#writeLines(c("Hello","World"), fileConn)
#close(fileConn)



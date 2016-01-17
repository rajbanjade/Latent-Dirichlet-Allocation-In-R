# Latent-Dirichlet-Allocation-In-R
 This is an R implementation of 
 Latent Dirichlet Allocation (LDA) using (Collapsed) Gibbs Sampling.

Please go through the code. I hope that the comments in the code itself give you a walk through of it.
Additional details maybe available, if not here, at https://umdrive.memphis.edu/rbanjade/public/

To run it, all you need is to:   
    1. provide an input file, each line containing a document and tokens separated by one or more whitespaces or tabs. A sample input file can be found along with this file. The test file contains texts from Microsoft Research Paraphrase (MSRP) corpus. I don't think it's a best corpus for topic modeling but I am using it just for demonstration purpose. You may test with some different input, preferably a bigger collection of documents.    
    2. set few parameters (optional). The default parameter values will be used otherwise.   

 I would suggest you review following topics before reading LDA.  
   (a) Distributions - Multinomial distribution, Dirichlet distribution  
   (b) Sampling - Gibbs Sampling  
   (c) Conjugate distributions  

 Main references:  
     https://www.youtube.com/watch?v=DDq3OVp9dNA (David Blei's talk)  
     http://u.cs.biu.ac.il/~89-680/darling-lda.pdf  


  Thanks!
 
  Rajendra Banjade  
   https://umdrive.memphis.edu/rbanjade/public/  
